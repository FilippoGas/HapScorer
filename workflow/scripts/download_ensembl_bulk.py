# workflow/scripts/download_ensembl_bulk.py

"""
Securely fetches and verifies a target reference dataset from the Ensembl server.

Automates the download of genomic annotations and FASTA/TSV files with network retry resilience.
It implements a rigorous two-stage data integrity pipeline:
1. Network transfer validation against upstream Ensembl CHECKSUMS utilizing the legacy UNIX `sum` algorithm.
2. Cryptographic SHA-256 digest generation for local storage provenance.

Supports dual execution modes:
1. Native Snakemake workflow integration via the `script:` directive (receives target URL via `params.url`
   and output destination via `output["dataset"]`).
2. Standalone Command Line Interface (CLI) execution via `argparse` for manual debugging.

Expected Inputs/Outputs:
    * **url** (`params.url` or `--url`): Direct Ensembl FTP/HTTP URL.
    * **out** (`output["dataset"]` or `--out`): Destination path for the downloaded dataset.

Optional Parameters:
    * **max_retries** (`params.max_retries` or `--max-retries`): Maximum download retry attempts (default: 3).
    * **backoff_factor** (`params.backoff_factor` or `--backoff-factor`): Exponential backoff delay factor in seconds (default: 5.0).
"""

import argparse
import hashlib
import logging
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Union

from utils import init_snakemake_log, setup_logger, setup_warnings

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake = Snakemake()

logger = setup_logger(__name__)
setup_warnings()

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def parse_args() -> argparse.Namespace:
    """
    Parses command-line arguments for single-dataset retrieval and verification.

    Returns:
        args (:obj:`argparse.Namespace`): An object containing parsed execution parameters.
    """
    parser = argparse.ArgumentParser(description="Securely fetches and verifies a single dataset from Ensembl.")
    parser.add_argument("--url", required=True, help="Direct FTP/HTTP URL of the Ensembl dataset.")
    parser.add_argument("--out", required=True, help="Output destination path for the dataset file.")
    parser.add_argument("--max-retries", type=int, default=3, help="Maximum download retry attempts (default: 3).")
    parser.add_argument(
        "--backoff-factor", type=float, default=5.0, help="Exponential backoff factor in seconds (default: 5.0)."
    )
    parser.add_argument("--verbose", action="store_true", help="Enable detailed DEBUG logging output.")
    return parser.parse_args()


def generate_sha256(file_path: Union[str, Path]) -> str:
    """
    Computes the SHA-256 cryptographic hash of a given file.

    This function reads the file in 4KB chunks, optimizing memory usage
    and allowing the hashing of large genomic datasets that exceed available RAM
    without causing memory overflow.

    Args:
        file_path (:obj:`str` or :obj:`pathlib.Path`): The local file system path to the target file.

    Returns:
        hash_digest (:obj:`str`): The resulting SHA-256 hash formatted as a hexadecimal string.
    """
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        # Read and update hash string value in blocks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def make_progress_hook(file_name: str) -> Callable[[int, int, int], None]:
    """
    Generates a callback function to monitor and report download progress.

    Utilizes a closure to track the last reported state, ensuring that standard
    output is not overwhelmed by logging progress iteratively. It logs updates
    in 10 MB increments.

    Args:
        file_name (:obj:`str`): The name of the file being downloaded, used for logging context.

    Returns:
        progress_hook (:obj:`typing.Callable`[[int, int, int], None]): A function matching the signature
            required by `urllib.request.urlretrieve`'s `reporthook` parameter.
    """
    last_reported_mb = -1

    def progress_hook(block_num: int, block_size: int, total_size: int):
        nonlocal last_reported_mb
        downloaded_bytes = block_num * block_size
        downloaded_mb = downloaded_bytes / (1024 * 1024)

        # Log progress every 10 MB or when complete to avoid spamming the log
        if int(downloaded_mb) % 10 == 0 and int(downloaded_mb) != last_reported_mb:
            last_reported_mb = int(downloaded_mb)
            if total_size > 0:
                total_mb = total_size / (1024 * 1024)
                percent = (downloaded_bytes / total_size) * 100
                logger.info(f"[{file_name}] Progress: {downloaded_mb:.1f} MB / {total_mb:.1f} MB ({percent:.1f}%)")
            else:
                logger.info(f"[{file_name}] Progress: {downloaded_mb:.1f} MB downloaded")

    return progress_hook


def fetch_and_verify(
    url: str,
    output_path: Union[str, Path],
    max_retries: int = 3,
    backoff_factor: float = 5.0,
) -> None:
    """
    Executes the secure download and dual-stage verification of a target dataset with network retry resilience.

    Orchestrates the retrieval of both the target dataset and its corresponding CHECKSUMS
    registry from the Ensembl FTP/HTTP service. Retries transient network failures using
    an exponential backoff strategy. Validates network transfer integrity via the system's
    legacy `sum` algorithm against upstream checksums, and computes a cryptographic SHA-256
    digest saved to an adjacent `.sha256` file upon successful verification.

    Args:
        url (:obj:`str`): Direct FTP or HTTP URL of the target Ensembl dataset.
        output_path (:obj:`str` or :obj:`pathlib.Path`): Destination path for the downloaded dataset file.
        max_retries (:obj:`int`, optional): Maximum number of retry attempts upon encountering network errors (default: 3).
        backoff_factor (:obj:`float`, optional): Multiplicative delay factor in seconds applied between retries (default: 5.0).

    Raises:
        ValueError: If the target file is omitted from upstream CHECKSUMS registry or checksum verification fails.
        Exception: Propagates underlying network or filesystem errors after exhausting retries and executing cleanup.
    """
    output_file = Path(output_path)
    file_name = output_file.name
    base_dir_url = url.rsplit("/", 1)[0]
    checksum_url = f"{base_dir_url}/CHECKSUMS"
    tmp_checksum_file = Path(f"{output_path}.checksums.tmp")

    logger.info(f"[{file_name}] Initiating dataset acquisition from: {url}")

    # Retry loop for primary file download
    for attempt in range(1, max_retries + 1):
        try:
            # Download the target data file
            hook = make_progress_hook(file_name)
            urllib.request.urlretrieve(url, output_file, reporthook=hook)
            logger.info(f"[{file_name}] Download completed successfully on attempt {attempt}/{max_retries}.")
            break
        except Exception as err:
            if output_file.exists():
                output_file.unlink()

            if attempt < max_retries:
                sleep_time = backoff_factor * (2 ** (attempt - 1))
                logger.warning(
                    f"[{file_name}] Transfer attempt {attempt}/{max_retries} failed: {err}. "
                    f"Retrying in {sleep_time:.1f} seconds..."
                )
                time.sleep(sleep_time)
            else:
                logger.error(f"[{file_name}] CRITICAL ERROR: All {max_retries} download attempts failed.")
                raise err

    # Retry loop for CHECKSUMS registry fetch
    for attempt in range(1, max_retries + 1):
        try:
            # Download the corresponding Ensembl CHECKSUMS file
            logger.info(f"[{file_name}] Fetching upstream CHECKSUMS registry...")
            urllib.request.urlretrieve(checksum_url, tmp_checksum_file)
            break
        except Exception as err:
            if tmp_checksum_file.exists():
                tmp_checksum_file.unlink()

            if attempt < max_retries:
                sleep_time = backoff_factor * (2 ** (attempt - 1))
                logger.warning(
                    f"[{file_name}] CHECKSUMS fetch attempt {attempt}/{max_retries} failed: {err}. "
                    f"Retrying in {sleep_time:.1f} seconds..."
                )
                time.sleep(sleep_time)
            else:
                logger.error(
                    f"[{file_name}] CRITICAL ERROR: Failed to retrieve CHECKSUMS registry after {max_retries} attempts."
                )
                if output_file.exists():
                    output_file.unlink()
                raise err

    try:
        logger.debug(f"[{file_name}] Parsing upstream CHECKSUMS registry entries...")

        # Parse the expected checksum (UNIX sum) from the downloaded text file
        expected_sum_line = ""
        with open(tmp_checksum_file, "r") as chk_f:
            for line in chk_f:
                if file_name in line:
                    expected_sum_line = line.strip()
                    break

        if not expected_sum_line:
            raise ValueError(f"Target file {file_name} not found in the upstream CHECKSUMS registry.")

        logger.debug(f"[{file_name}] Matched upstream CHECKSUMS line: '{expected_sum_line}'")
        expected_checksum_value, expected_blocks = expected_sum_line.split()[:2]

        # Compute the local UNIX sum using the system subsystem
        logger.info(f"[{file_name}] Verifying network transfer integrity (Legacy 'sum')...")
        logger.debug(f"[{file_name}] Executing subprocess command: 'sum {output_file}'")
        sum_process = subprocess.run(["sum", str(output_file)], capture_output=True, text=True, check=True)
        actual_checksum_value, actual_blocks = sum_process.stdout.strip().split()[:2]

        if expected_checksum_value != actual_checksum_value or expected_blocks != actual_blocks:
            raise ValueError(
                f"Checksum mismatch for {file_name}.\n"
                f"Expected: {expected_checksum_value} (Blocks: {expected_blocks})\n"
                f"Actual:   {actual_checksum_value} (Blocks: {actual_blocks})"
            )

        logger.info(f"[{file_name}] Transfer integrity verified successfully.")

        # Generate a robust cryptographic hash for local provenance
        logger.info(f"[{file_name}] Generating secure SHA-256 digest for local data provenance...")
        sha256_digest = generate_sha256(output_file)
        sha256_file = Path(f"{output_path}.sha256")
        with open(sha256_file, "w") as sha_f:
            sha_f.write(f"{sha256_digest}  {file_name}\n")

        logger.info(f"[{file_name}] SHA-256 digest generated: {sha256_digest}\n")

    except Exception as e:
        logger.error(f"[{file_name}] CRITICAL ERROR: {str(e)}")
        # Clean up potentially corrupted files upon failure
        if output_file.exists():
            output_file.unlink()
        raise e
    finally:
        # Guarantee cleanup of temporary checksum file
        if tmp_checksum_file.exists():
            tmp_checksum_file.unlink()


# =============================================================================
# MAIN LOGIC
# =============================================================================


def run_pipeline(
    url: str,
    out_path: Union[str, Path],
    max_retries: int = 3,
    backoff_factor: float = 5.0,
) -> None:
    """
    Executes the secure acquisition, retry-backed download, and verification of a target dataset.

    Args:
        url (:obj:`str`): Direct FTP or HTTP URL of the target Ensembl dataset.
        out_path (:obj:`str` or :obj:`pathlib.Path`): Destination file path for the downloaded dataset.
        max_retries (:obj:`int`, optional): Maximum retry attempts for transient network errors (default: 3).
        backoff_factor (:obj:`float`, optional): Exponential backoff delay factor in seconds (default: 5.0).
    """

    target_name = Path(out_path).name
    logger.info(f"Initiating target dataset acquisition and verification for '{target_name}' (source: {url})...")

    fetch_and_verify(
        url=url,
        output_path=out_path,
        max_retries=max_retries,
        backoff_factor=backoff_factor,
    )

    logger.info(f"Target dataset '{target_name}' successfully retrieved and verified.")


if __name__ == "__main__":
    is_snakemake = False

    if "snakemake" in globals():
        is_snakemake = True
        init_snakemake_log(snakemake, base_dir=Path.cwd())
    else:
        args = parse_args()

    url = snakemake.params["url"] if is_snakemake else args.url
    out_path = snakemake.output["dataset"] if is_snakemake else args.out
    max_retries = snakemake.params.get("max_retries", 3) if is_snakemake else args.max_retries
    backoff_factor = snakemake.params.get("backoff_factor", 5.0) if is_snakemake else args.backoff_factor

    if not is_snakemake and args.verbose:
        logger.setLevel(logging.DEBUG)

    run_pipeline(
        url=url,
        out_path=out_path,
        max_retries=max_retries,
        backoff_factor=backoff_factor,
    )
