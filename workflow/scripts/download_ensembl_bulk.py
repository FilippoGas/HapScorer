# workflow/scripts/download_ensembl_bulk.py

"""
Securely fetches and verifies bulk reference datasets from the Ensembl server.

Designed for execution within a Snakemake workflow, this script automates the
download of genomic annotations and FASTA files. It implements a rigorous two-stage
data integrity pipeline:
1. Network transfer validation against upstream Ensembl CHECKSUMS utilizing the
   legacy UNIX `sum` algorithm.
2. Cryptographic SHA-256 digest generation for local storage provenance, mitigating
   downstream risks associated with bit rot or unintended modifications.

Expected Snakemake Parameters:
    - release (str or int): The Ensembl release version (e.g., 113).

Expected Snakemake Outputs:
    - gtf (str or Path): Output path for the GTF annotation file.
    - cds (str or Path): Output path for the CDS FASTA file.
    - pep (str or Path): Output path for the peptide FASTA file.
    - cdna (str or Path): Output path for the cDNA FASTA file.
    - uniprot (str or Path): Output path for the UniProt TSV mapping file.
"""

import hashlib
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Union

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake = Snakemake()

# Redirect standard output and standard error to the Snakemake log file
log_file = open(snakemake.log[0], "w", buffering=1)
sys.stdout = log_file
sys.stderr = log_file


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def generate_sha256(file_path: Union[str, Path]) -> str:
    """
    Computes the SHA-256 cryptographic hash of a given file.

    This function reads the file in 4KB chunks, optimizing memory usage
    and allowing the hashing of large genomic datasets that exceed available RAM
    without causing memory overflow.

    Args:
        file_path (str or pathlib.Path): The local file system path to the target file.

    Returns:
        str: The resulting SHA-256 hash formatted as a hexadecimal string.
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
        file_name (str): The name of the file being downloaded, used for logging context.

    Returns:
        Callable[[int, int, int], None]: A function matching the signature required by
        `urllib.request.urlretrieve`'s `reporthook` parameter.
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
                print(
                    f"[{file_name}] Progress: {downloaded_mb:.1f} MB / {total_mb:.1f} MB ({percent:.1f}%)", flush=True
                )
            else:
                print(f"[{file_name}] Progress: {downloaded_mb:.1f} MB downloaded", flush=True)

    return progress_hook


def fetch_and_verify(url: str, output_path: Union[str, Path]) -> None:
    """
    Executes the secure download and dual-stage verification of a target dataset.

    The function orchestrates the download of both the target data and its
    corresponding CHECKSUMS registry. It validates the network transfer using
    the system's native `sum` command, and upon success, computes a robust
    SHA-256 digest saved to an adjacent `.sha256` file. Any integrity failure
    triggers an automatic cleanup of the corrupted data.

    Args:
        url (str): The direct FTP/HTTP URL of the Ensembl dataset.
        output_path (str or pathlib.Path): The local file system path where the
            dataset should be stored.

    Raises:
        ValueError: If the file is missing from the CHECKSUMS registry or if
            the calculated UNIX sum mismatches the upstream expectation.
        Exception: Propagates any underlying network or filesystem errors after
            performing cleanup of temporary or corrupted files.
    """
    output_file = Path(output_path)
    file_name = output_file.name
    base_dir_url = url.rsplit("/", 1)[0]
    checksum_url = f"{base_dir_url}/CHECKSUMS"
    tmp_checksum_file = Path(f"{output_path}.checksums.tmp")

    print(f"[{file_name}] Initiating download from: {url}", flush=True)

    try:
        # 1. Download the target data file
        hook = make_progress_hook(file_name)
        urllib.request.urlretrieve(url, output_file, reporthook=hook)
        print(f"[{file_name}] Download completed.", flush=True)

        # 2. Download the corresponding Ensembl CHECKSUMS file
        print(f"[{file_name}] Fetching upstream CHECKSUMS file...", flush=True)
        urllib.request.urlretrieve(checksum_url, tmp_checksum_file)

        # 3. Parse the expected checksum (UNIX sum) from the downloaded text file
        expected_sum_line = ""
        with open(tmp_checksum_file, "r") as chk_f:
            for line in chk_f:
                if file_name in line:
                    expected_sum_line = line.strip()
                    break

        if not expected_sum_line:
            raise ValueError(f"Target file {file_name} not found in the upstream CHECKSUMS registry.")

        expected_checksum_value = expected_sum_line.split()[0]
        expected_blocks = expected_sum_line.split()[1]

        # 4. Compute the local UNIX sum using the system subsystem
        print(f"[{file_name}] Verifying network transfer integrity (Legacy 'sum')...", flush=True)
        sum_process = subprocess.run(["sum", str(output_file)], capture_output=True, text=True, check=True)
        actual_checksum_value, actual_blocks = sum_process.stdout.strip().split()[:2]

        if expected_checksum_value != actual_checksum_value or expected_blocks != actual_blocks:
            raise ValueError(
                f"Checksum mismatch for {file_name}.\n"
                f"Expected: {expected_checksum_value} (Blocks: {expected_blocks})\n"
                f"Actual:   {actual_checksum_value} (Blocks: {actual_blocks})"
            )

        print(f"[{file_name}] Transfer integrity verified successfully.", flush=True)

        # 5. Generate a robust cryptographic hash for local provenance
        print(f"[{file_name}] Generating secure SHA-256 digest for local data provenance...", flush=True)
        sha256_digest = generate_sha256(output_file)
        sha256_file = Path(f"{output_path}.sha256")
        with open(sha256_file, "w") as sha_f:
            sha_f.write(f"{sha256_digest}  {file_name}\n")

        print(f"[{file_name}] SHA-256 digest generated: {sha256_digest}\n", flush=True)

    except Exception as e:
        print(f"[{file_name}] CRITICAL ERROR: {str(e)}", file=sys.stderr, flush=True)
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


def main():
    """
    Primary execution function for Snakemake integration.
    """
    # Extract execution parameters dynamically from the Snakemake object
    release = snakemake.params.release

    # Map output variables to their exact Ensembl URL pathways
    # These URLs reflect the canonical directory structure of the Ensembl FTP service
    datasets = {
        snakemake.output.gtf: f"https://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{release}.gtf.gz",
        snakemake.output.cds: f"https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz",
        snakemake.output.pep: f"https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz",
        snakemake.output.cdna: f"https://ftp.ensembl.org/pub/release-{release}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        snakemake.output.uniprot: f"https://ftp.ensembl.org/pub/release-{release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{release}.uniprot.tsv.gz",
    }

    print("===============================================================", flush=True)
    print(f" Commencing data acquisition for Ensembl Release {release}", flush=True)
    print("===============================================================\n", flush=True)

    for output_path, origin_url in datasets.items():
        fetch_and_verify(origin_url, output_path)

    print("======================================================================", flush=True)
    print(" All requested datasets have been successfully retrieved and verified.", flush=True)
    print("======================================================================\n", flush=True)


if __name__ == "__main__":
    main()
