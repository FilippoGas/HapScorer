# workflow/scripts/generate_mutated_sequences.py

"""
Generates mutated coding nucleotide sequences (CDS) for each haplotype.

Parses haplotype variant signatures (SNPs, insertions, deletions) and applies
them to canonical transcript reference sequences (5' UTR, CDS, 3' UTR). Accounts
for start-codon gains in the 5' UTR and stop-codon losses extending into the 3' UTR.

Supports dual execution modes:
1. Native Snakemake workflow integration via the `script:` directive.
2. Standalone Command Line Interface (CLI) execution via `argparse` for manual debugging.

Expected Snakemake Inputs (or CLI Arguments):
    - haplotypes (--haplotypes): Path to the extracted haplotypes CSV file.
    - wt_cds (--wt-cds): Path to the canonical reference CDS sequences CSV file.
    - annotations (--annotations): Path to transcript structural annotations CSV file.

Expected Snakemake Outputs (or CLI Arguments):
    - out_mutated (--out-mutated): Path for the generated mutated nucleotide sequences CSV.

Optional Parameters:
    - threads (--threads): Number of CPU cores to utilize for multiprocessing (default: 1).
"""

import argparse
import multiprocessing as mp
import re
import sys
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import pandas as pd
from Bio.Seq import translate

if TYPE_CHECKING:
    from snakemake.script import Snakemake

    snakemake: Snakemake = Snakemake()

# =============================================================================
# HELPER FUNCTIONS & VARIANT PARSING
# =============================================================================


def parse_args() -> argparse.Namespace:
    """
    Parses command-line arguments for the mutated sequence generation workflow.

    Returns:
        argparse.Namespace: An object containing the parsed execution parameters.
    """
    parser = argparse.ArgumentParser(description="Generates mutated nucleotide sequences for haplotypes.")
    parser.add_argument("--haplotypes", required=True, help="Path to input haplotypes CSV.")
    parser.add_argument("--wt-cds", required=True, help="Path to input wild-type CDS CSV.")
    parser.add_argument("--annotations", required=True, help="Path to transcript annotations CSV.")
    parser.add_argument("--out-mutated", required=True, help="Output path for mutated sequences CSV.")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads to utilize.")
    return parser.parse_args()


def parse_variant(variant_type: str, nuc_change: str) -> Optional[Dict[str, Any]]:
    """
    Parses an HGVS-like nucleotide change string into a structured dictionary.

    Deconstructs HGVS notation for single nucleotide polymorphisms (SNPs),
    deletions, insertions, and duplications across 5' UTR, CDS, and 3' UTR
    regions while filtering out intronic variants containing '+'.

    Args:
        variant_type (str): The variant consequence type (e.g., 'missense_variant').
        nuc_change (str): HGVS nucleotide change string (e.g., 'c.123A>G', 'c.-15G>A', 'c.*10del').

    Returns:
        Optional[Dict[str, Any]]: A dictionary detailing parsed variant attributes:
            - variant_type (str): Original variant consequence type.
            - nucleotide_change (str): Original HGVS string.
            - pos (int): 1-based relative nucleotide position.
            - ref (str): Reference nucleotide allele.
            - alt (str): Alternative nucleotide allele.
            - desc (str): Variant class ('snp', 'del', or 'ins').
            Returns None if the variant is intronic or unparseable.
    """
    if not nuc_change or "+" in nuc_change:
        return None  # Exclude intronic variants

    is_5utr = "-" in nuc_change
    is_3utr = "*" in nuc_change

    if "del" in nuc_change:
        desc = "del"
        ref = nuc_change.split("del")[-1]
        alt = ""
        sep = "_-" if is_5utr else ("_*" if is_3utr else "_")
        parts = nuc_change.split("del")[0].split(sep)
        pos_str = (
            parts[-1]
            if len(parts) > 1
            else nuc_change.split("del")[0].replace("-", "").replace("*", "").replace("c.", "")
        )

        try:
            raw_pos = int(re.sub(r"[^\d]", "", pos_str))
            pos = (raw_pos + len(ref) - 1) if is_5utr else (raw_pos - len(ref) + 1)
        except ValueError:
            return None

    elif "ins" in nuc_change:
        desc = "ins"
        ref = ""
        alt = nuc_change.split("ins")[-1]
        sep = "_-" if is_5utr else ("_*" if is_3utr else "_")
        parts = nuc_change.split("ins")[0].split(sep)
        pos_str = (
            parts[-1]
            if len(parts) > 1
            else nuc_change.split("ins")[0].replace("-", "").replace("*", "").replace("c.", "")
        )

        try:
            raw_pos = int(re.sub(r"[^\d]", "", pos_str))
            pos = (raw_pos + 1) if is_5utr else raw_pos
        except ValueError:
            return None

    elif "dup" in nuc_change:
        desc = "ins"
        ref = ""
        alt = nuc_change.split("dup")[-1]
        sep = "_-" if is_5utr else ("_*" if is_3utr else "_")
        parts = nuc_change.split("dup")[0].split(sep)
        pos_str = (
            parts[-1]
            if len(parts) > 1
            else nuc_change.split("dup")[0].replace("-", "").replace("*", "").replace("c.", "")
        )

        try:
            raw_pos = int(re.sub(r"[^\d]", "", pos_str))
            pos = (raw_pos + len(alt) - 1) if is_5utr else (raw_pos - len(alt) + 1)
        except ValueError:
            return None

    elif ">" in nuc_change:
        desc = "snp"
        parts = nuc_change.split(">")
        alt = parts[1]
        ref = parts[0][-1]
        pos_str = parts[0][:-1]
        try:
            pos = int(re.sub(r"[^\d]", "", pos_str))
        except ValueError:
            return None
    else:
        return None

    return {
        "variant_type": variant_type,
        "nucleotide_change": nuc_change,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "desc": desc,
    }


def apply_snp(sequence: str, snps: List[Dict[str, Any]], is_5utr: bool) -> str:
    """
    Applies single nucleotide substitutions (SNPs) to a target sequence string.

    Modifies characters in place according to 1-based positions, correctly handling
    directional indexing differences between 5' UTR (reverse 3'-to-5' offset)
    and CDS/3' UTR regions.

    Args:
        sequence (str): The target nucleotide sequence string.
        snps (List[Dict[str, Any]]): List of parsed SNP variant dictionaries.
        is_5utr (bool): Flag indicating whether the sequence belongs to the 5' UTR region.

    Returns:
        str: The updated nucleotide sequence string containing applied substitutions.
    """
    seq_list = list(sequence)
    seq_len = len(sequence)
    for snp in snps:
        pos = snp["pos"]
        idx = (seq_len - pos) if is_5utr else (pos - 1)
        if 0 <= idx < seq_len:
            seq_list[idx] = snp["alt"]
    return "".join(seq_list)


def apply_deletion(sequence: str, dels: List[Dict[str, Any]], is_5utr: bool) -> str:
    """
    Applies deletion events to a target sequence by masking deleted positions.

    Replaces deleted nucleotide bases with placeholder characters ('D') to preserve
    relative positional indexing prior to downstream string cleanup.

    Args:
        sequence (str): The target nucleotide sequence string.
        dels (List[Dict[str, Any]]): List of parsed deletion variant dictionaries.
        is_5utr (bool): Flag indicating whether the sequence belongs to the 5' UTR region.

    Returns:
        str: The sequence string containing 'D' placeholder masks at deleted sites.
    """
    seq_list = list(sequence)
    seq_len = len(sequence)
    for d in dels:
        pos = d["pos"]
        ref_len = len(d["ref"])
        start = (seq_len - pos) if is_5utr else (pos - 1)
        end = start + ref_len
        for i in range(max(0, start), min(seq_len, end)):
            seq_list[i] = "D"
    return "".join(seq_list)


def apply_insertion(sequence: str, inss: List[Dict[str, Any]], is_5utr: bool) -> str:
    """
    Inserts alternative nucleotide sequences into a target sequence string.

    Sorts insertion events by position and inserts alternative alleles into
    the target sequence string, accounting for 5' UTR coordinate mapping.

    Args:
        sequence (str): The target nucleotide sequence string.
        inss (List[Dict[str, Any]]): List of parsed insertion or duplication variant dictionaries.
        is_5utr (bool): Flag indicating whether the sequence belongs to the 5' UTR region.

    Returns:
        str: The updated sequence string containing inserted nucleotide segments.
    """
    sorted_inss = sorted(inss, key=lambda x: x["pos"])
    seq = sequence
    for ins in sorted_inss:
        pos = ins["pos"]
        idx = (len(seq) - pos + 1) if is_5utr else pos
        if 0 <= idx <= len(seq):
            seq = seq[:idx] + ins["alt"] + seq[idx:]
    return seq


def compute_shift(shift: List[int], inss: List[Dict[str, Any]], dels: List[Dict[str, Any]]) -> List[int]:
    """
    Calculates coordinate positional shift vectors across 5' UTR insertion and deletion events.

    Tracks accumulated length adjustments induced by indels within 5' UTR sequences
    to evaluate positional shifts of downstream canonical start codons (ATG).

    Args:
        shift (List[int]): Initial coordinate shift vector initialized with zeros.
        inss (List[Dict[str, Any]]): List of parsed insertion variants in the 5' UTR.
        dels (List[Dict[str, Any]]): List of parsed deletion variants in the 5' UTR.

    Returns:
        List[int]: Updated shift vector mapping original positions to length offsets,
            with deleted positions represented as None.
    """
    shift = list(shift)
    for ins in inss:
        start = len(shift) - ins["pos"]
        alt_len = len(ins["alt"])
        if 0 <= start < len(shift):
            shift = shift[:start] + [shift[start] + alt_len] * (len(shift) - start)
    for d in dels:
        start = len(shift) - d["pos"]
        ref_len = len(d["ref"])
        end = start + ref_len + 1
        if 0 <= start < len(shift):
            shift = shift[:start] + [None] * ref_len + ([s - ref_len for s in shift[end:]] if end < len(shift) else [])
    return shift


# =============================================================================
# CORE SEQUENCE MUTATION LOGIC
# =============================================================================


def build_mutated_sequence(wt_sequence: Dict[str, str], variant_types_str: str, nuc_change_str: str) -> Dict[str, Any]:
    """
    Constructs the full mutated nucleotide sequence for a specific haplotype signature.

    Evaluates and applies variants across 5' UTR, CDS, and 3' UTR regions:
    1. Modifies 5' UTR and detects newly gained start codons (ATG) that initiate translation.
    2. Applies non-synonymous CDS variants (SNPs, insertions, deletions).
    3. Translates the resulting sequence and extends into the 3' UTR if the canonical
       stop codon is lost due to mutation.

    Args:
        wt_sequence (Dict[str, str]): Dictionary containing reference sequence domains:
            '5utr', 'coding', and '3utr'.
        variant_types_str (str): Comma-separated string of variant consequence types.
        nuc_change_str (str): Comma-separated string of HGVS nucleotide changes.

    Returns:
        Dict[str, Any]: A dictionary containing mutated sequence evaluation metrics:
            - nn_sequence (str): Reconstructed mutated nucleotide sequence.
            - start_gain (bool): Indicates if a novel start codon was gained in the 5' UTR.
            - multiple_starts_in_5 (bool): Indicates presence of multiple ATG sites in 5' UTR.
            - stop_loss_3_translated (bool): Indicates translation extension into 3' UTR due to stop-loss
    """
    v_types = [v.strip() for v in str(variant_types_str).split(",") if v.strip()]
    n_changes = [n.strip() for n in str(nuc_change_str).split(",") if n.strip()]

    parsed_variants = []
    for vt, nc in zip(v_types, n_changes):
        v_info = parse_variant(vt, nc)
        if v_info:
            parsed_variants.append(v_info)

    utr5_seq = str(wt_sequence.get("5utr", "")) if pd.notna(wt_sequence.get("5utr")) else ""
    coding_seq = str(wt_sequence.get("coding", "")) if pd.notna(wt_sequence.get("coding")) else ""
    utr3_seq = str(wt_sequence.get("3utr", "")) if pd.notna(wt_sequence.get("3utr")) else ""

    sequence = ""
    start_gain = False
    multiple_starts_in_5 = False
    stop_loss_3_translated = False

    # Process 5' UTR
    vars_5utr = [v for v in parsed_variants if "-" in v["nucleotide_change"] and v["pos"] < len(utr5_seq)]
    if vars_5utr and utr5_seq:
        utr5_mod = utr5_seq
        snps = [v for v in vars_5utr if v["desc"] == "snp"]
        dels = [v for v in vars_5utr if v["desc"] == "del"]
        inss = [v for v in vars_5utr if v["desc"] == "ins"]

        shift = [0] * len(utr5_seq)
        if dels or inss:
            shift = compute_shift(shift, inss, dels)

        if snps:
            utr5_mod = apply_snp(utr5_mod, snps, is_5utr=True)
        if dels:
            utr5_mod = apply_deletion(utr5_mod, dels, is_5utr=True)
        if inss:
            utr5_mod = apply_insertion(utr5_mod, inss, is_5utr=True)

        utr5_mod = utr5_mod.replace("D", "")

        if "ATG" in utr5_mod:
            wt_starts = [m.start() for m in re.finditer("ATG", utr5_seq)]
            wt_shifted_starts = {(st + shift[st]) for st in wt_starts if st < len(shift) and shift[st] is not None}
            mod_starts = [m.start() for m in re.finditer("ATG", utr5_mod)]
            multiple_starts_in_5 = True

            new_starts = [st for st in mod_starts if st not in wt_shifted_starts]
            if new_starts:
                first_start = min(new_starts)
                sequence = utr5_mod[first_start:]
                start_gain = True

    # Process CDS
    cds_vars = [v for v in parsed_variants if "-" not in v["nucleotide_change"] and "*" not in v["nucleotide_change"]]
    if cds_vars:
        cds_mod = coding_seq
        snps = [v for v in cds_vars if v["desc"] == "snp"]
        dels = [v for v in cds_vars if v["desc"] == "del"]
        inss = [v for v in cds_vars if v["desc"] == "ins"]

        if snps:
            cds_mod = apply_snp(cds_mod, snps, is_5utr=False)
        if dels:
            cds_mod = apply_deletion(cds_mod, dels, is_5utr=False)
        if inss:
            cds_mod = apply_insertion(cds_mod, inss, is_5utr=False)

        cds_mod = cds_mod.replace("D", "")

        if sequence:
            sequence += cds_mod
        else:
            match = re.search("ATG", cds_mod)
            if match:
                sequence = cds_mod[match.start() :]
    else:
        sequence += coding_seq

    # Process 3' UTR
    try:
        translated = str(translate(sequence))
    except Exception:
        translated = ""

    if "*" not in translated:
        stop_loss_3_translated = True
        vars_3utr = [v for v in parsed_variants if "*" in v["nucleotide_change"]]
        if vars_3utr and utr3_seq:
            utr3_mod = utr3_seq
            snps = [v for v in vars_3utr if v["desc"] == "snp"]
            dels = [v for v in vars_3utr if v["desc"] == "del"]
            inss = [v for v in vars_3utr if v["desc"] == "ins"]

            if snps:
                utr3_mod = apply_snp(utr3_mod, snps, is_5utr=False)
            if dels:
                utr3_mod = apply_deletion(utr3_mod, dels, is_5utr=False)
            if inss:
                utr3_mod = apply_insertion(utr3_mod, inss, is_5utr=False)

            utr3_mod = utr3_mod.replace("D", "")
            sequence += utr3_mod
        else:
            sequence += utr3_seq

    return {
        "nn_sequence": sequence,
        "start_gain": start_gain,
        "multiple_starts_in_5": multiple_starts_in_5,
        "stop_loss_3_translated": stop_loss_3_translated,
    }


def process_row(row: Tuple[Any, pd.Series], sequences_dict: Dict[str, Dict[str, str]]) -> Dict[str, Any]:
    """
    Worker function for parallelized row-by-row haplotype sequence reconstruction.

    Distinguishes wild-type haplotypes (.0) from mutated variants and delegates
    mutated sequence construction to `build_mutated_sequence`.

    Args:
        row (Tuple[Any, pandas.Series]): Tuple containing row index and DataFrame Series
            representing a single haplotype record.
        sequences_dict (Dict[str, Dict[str, str]]): Lookup dictionary mapping transcript IDs
            to wild-type domain sequences ('coding', '5utr', '3utr').

    Returns:
        Dict[str, Any]: Dictionary containing the haplotype transcript ID and its
            reconstructed mutated nucleotide sequence metadata.
    """
    _, data = row
    hap_id = str(data["haplotype_transcript_id"])
    tx_id = str(data["ensembl_transcript_id"])
    wt_seq = sequences_dict.get(tx_id, {"coding": "", "5utr": "", "3utr": ""})

    if str(data["haplotype_gene_id"]).endswith(".0"):
        res = {
            "haplotype_transcript_id": hap_id,
            "start_gain": False,
            "multiple_starts_in_5": False,
            "stop_loss_3_translated": False,
            "nn_sequence": wt_seq.get("coding", ""),
        }
    else:
        res_mut = build_mutated_sequence(
            wt_sequence=wt_seq,
            variant_types_str=data.get("variant_types", ""),
            nuc_change_str=data.get("nucleotide_change", ""),
        )
        res = {"haplotype_transcript_id": hap_id, **res_mut}

    return res


def run_pipeline(
    haplotypes_path: Union[str, Path],
    wt_cds_path: Union[str, Path],
    annotations_path: Union[str, Path],
    out_mutated_path: Union[str, Path],
    threads: int,
) -> None:
    """
    Executes the complete mutated sequence generation pipeline.

    Loads input datasets, filters reference sequences for the targeted chromosome,
    parallelizes haplotype processing across worker threads using `multiprocessing.Pool`,
    merges sequence outputs back into the haplotype metadata, and exports the final CSV dataset.

    Args:
        haplotypes_path (str or pathlib.Path): Path to the input haplotypes metadata CSV file.
        wt_cds_path (str or pathlib.Path): Path to the reference wild-type CDS sequences CSV file.
        annotations_path (str or pathlib.Path): Path to the transcript structural annotations CSV file.
        out_mutated_path (str or pathlib.Path): Output file path for the mutated sequences CSV file.
        threads (int): Number of CPU worker threads to spawn for multiprocessing.

    Returns:
        None
    """
    haplotypes_path = Path(haplotypes_path)
    wt_cds_path = Path(wt_cds_path)
    annotations_path = Path(annotations_path)
    out_mutated_path = Path(out_mutated_path)

    df_haplo = pd.read_csv(haplotypes_path)
    df_wt = pd.read_csv(wt_cds_path)
    df_anno = pd.read_csv(annotations_path)

    chrom_from_file = haplotypes_path.name.split("chr")[-1].split("_")[0]
    anno_txs = set(df_anno[df_anno["chromosome_name"].astype(str) == chrom_from_file]["ensembl_transcript_id"].unique())

    df_wt = df_wt[df_wt["ensembl_transcript_id"].isin(anno_txs)]
    sequences_dict = df_wt.set_index("ensembl_transcript_id")[["coding", "5utr", "3utr"]].to_dict("index")

    pool = mp.Pool(threads)
    worker_func = partial(process_row, sequences_dict=sequences_dict)
    results = pool.map(worker_func, df_haplo.iterrows())
    pool.close()
    pool.join()

    df_res = pd.DataFrame(results)
    df_out = pd.merge(df_haplo, df_res, on="haplotype_transcript_id", how="left")

    # Uppercase Boolean formatting for parity with R
    bool_cols = ["start_gain", "multiple_starts_in_5", "stop_loss_3_translated"]
    for col in bool_cols:
        df_out[col] = df_out[col].map({True: "TRUE", False: "FALSE"})

    # Explicit reordering of columns
    cols_order = [
        "haplotype_transcript_id",
        "ensembl_gene_id",
        "ensembl_transcript_id",
        "haplotype_gene_id",
        "variants",
        "rsid",
        "variant_types",
        "hgnc_symbol",
        "strand",
        "start_gain",
        "multiple_starts_in_5",
        "stop_loss_3_translated",
        "nn_sequence",
    ]
    df_out = df_out[cols_order]

    out_mutated_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_mutated_path, index=False, na_rep="NA")


# =============================================================================
# MAIN LOGIC
# =============================================================================

if __name__ == "__main__":
    is_snakemake = False

    if "snakemake" in globals():
        is_snakemake = True
        if hasattr(snakemake, "log") and snakemake.log:
            log_f = open(snakemake.log[0], "w", buffering=1)
            sys.stdout = log_f
            sys.stderr = log_f
    else:
        args = parse_args()

    haplotypes_path = snakemake.input["haplotypes_csv"] if is_snakemake else args.haplotypes
    wt_cds_path = snakemake.input["wt_cds"] if is_snakemake else args.wt_cds
    annotations_path = snakemake.input["annotations"] if is_snakemake else args.annotations
    out_mutated_path = snakemake.output["mutated_cds_csv"] if is_snakemake else args.out_mutated
    threads = snakemake.threads if is_snakemake else args.threads

    run_pipeline(
        haplotypes_path=haplotypes_path,
        wt_cds_path=wt_cds_path,
        annotations_path=annotations_path,
        out_mutated_path=out_mutated_path,
        threads=threads,
    )
