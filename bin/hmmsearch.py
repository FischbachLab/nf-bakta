# Description: This script is used to search for Pfam domains in the predicted proteins of the CDSs.
import argparse
import logging
from collections import OrderedDict
from typing import Optional, Sequence

import polars as pl
import pyhmmer
from pydantic import AliasChoices, BaseModel, ConfigDict, Field

log = logging.getLogger(__name__)


class Pfam(BaseModel):
    model_config = ConfigDict(coerce_numbers_to_str=True)
    pfam_id: str
    feature_id: str
    pfam_start: int = Field(
        serialization_alias="pfam_start",
        validation_alias=AliasChoices("start", "pfam_start"),
    )
    pfam_end: int = Field(
        serialization_alias="pfam_end",
        validation_alias=AliasChoices("end", "pfam_end", "pfam_stop", "stop"),
    )
    pfam_len: Optional[int | None] = None
    aln_len: Optional[int | None] = None
    evalue: Optional[float | None] = None
    score: Optional[float | None] = None
    aa_cov: Optional[float | None] = None
    domain_cov: Optional[float | None] = None

    def to_dict(self) -> dict:
        return self.model_dump(by_alias=True, exclude_none=True)


# function to read a protein fasta file using biopython with minimum sequence length greater than 10 amino acids
def read_protein_fasta(fasta_file: str) -> Sequence[pyhmmer.easel.DigitalSequence]:
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(
        fasta_file, digital=True, alphabet=alphabet
    ) as seq_file:
        sequences = seq_file.read_block()

    return sequences


def run_hmmsearch(
    proteins_fasta: str, pfam_db_path: str, output_file: str, threads: int
) -> Sequence[dict]:
    """Detect Pfam-A entries"""
    pfam_hits = []

    proteins = read_protein_fasta(proteins_fasta)

    with pyhmmer.plan7.HMMFile(pfam_db_path) as hmm, open(output_file, "wb") as fh:
        for top_hits in pyhmmer.hmmsearch(
            hmm, proteins, bit_cutoffs="gathering", cpus=threads, Z=1000000
        ):
            for hit in top_hits:
                cds_len = hit.best_domain.alignment.target_length
                domain_cov = (
                    hit.best_domain.alignment.hmm_to
                    - hit.best_domain.alignment.hmm_from
                    + 1
                ) / len(hit.best_domain.alignment.hmm_sequence)
                aa_cov = (
                    hit.best_domain.alignment.target_to
                    - hit.best_domain.alignment.target_from
                    + 1
                ) / cds_len

                pfam = OrderedDict()
                pfam["pfam_id"] = hit.best_domain.alignment.hmm_accession.decode()
                pfam["feature_id"] = hit.name.decode()
                pfam["start"] = hit.best_domain.alignment.target_from
                pfam["stop"] = hit.best_domain.alignment.target_to
                pfam["aln_len"] = len(hit.best_domain.alignment.hmm_sequence)
                pfam["pfam_len"] = hit.best_domain.alignment.hmm_length
                pfam["evalue"] = hit.evalue
                pfam["score"] = hit.score
                pfam["name"] = hit.best_domain.alignment.hmm_name.decode()
                pfam["aa_cov"] = aa_cov
                pfam["domain_cov"] = domain_cov

                pfam_hits.append(Pfam(**pfam).to_dict())

            top_hits.write(fh, format="domains", header=True)
    return pfam_hits


def compute_overlap(segment1, segment2):
    """
    # Example usage:
    segment1 = (10, 40)
    segment2 = (20, 30)
    print(calculate_percentage_overlap(segment1, segment2))
    """
    start1, end1 = segment1
    start2, end2 = segment2

    # Calculate overlap
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_end > overlap_start:
        # There's an overlap
        overlap_length = overlap_end - overlap_start
    else:
        # No overlap
        return 0

    # Length of the shorter segment
    length_shorter_segment = min(end1 - start1, end2 - start2)

    # Calculate percentage overlap
    percentage_overlap = (overlap_length / length_shorter_segment) * 100

    return percentage_overlap


def usage():
    parser = argparse.ArgumentParser(description="Pfam domain search")
    parser.add_argument(
        "-i",
        "--input",
        help="Input protein file in fasta format",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file prefix",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--pfam_db",
        help="Pfam database file",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads",
        default=12,
        type=int,
    )
    parser.add_argument(
        "--max_evalue",
        help="Maximum evalue threshold",
        default=1e-5,
        type=float,
    )
    parser.add_argument(
        "--min_domain_coverage",
        help="Minimum domain coverage threshold [0-1]",
        default=0.4,
        type=float,
    )
    parser.add_argument(
        "--min_overlap",
        help="Minimum overlap threshold [0-100]",
        default=50,
        type=float,
    )
    return parser.parse_args()


def main():
    args = usage()
    protein_file = args.input
    output_prefix = args.output
    pfam_db = args.pfam_db

    threads = args.threads
    max_evalue = args.max_evalue
    min_domain_coverage = args.min_domain_coverage
    min_overlap = args.min_overlap

    ## TEST DATA
    # protein_file = "test_data/proteins/Clostridium-sporogenes-ATCC-15579-MAF-2.faa"
    # output_prefix = "test_data/proteins/Clostridium-sporogenes-ATCC-15579-MAF-2"
    # pfam_db = "/mnt/efs/databases/HMMs/pfam/v36.0/Pfam-A.hmm"

    # threads = 12
    # max_evalue = 1e-5
    # min_domain_coverage = 0.4
    # min_overlap = 50

    domtbl_output = f"{output_prefix}.hmmsearch_domtblout.tsv"
    processed_output = f"{output_prefix}.hmmsearch_all_hits.csv"
    drop_candidates_output = f"{output_prefix}.hmmsearch_dropped_hits.csv"
    final_output = f"{output_prefix}.hmmsearch_tbl.csv"

    hits = (
        pl.DataFrame(run_hmmsearch(protein_file, pfam_db, domtbl_output, threads))
        .sort(["feature_id", "pfam_start"], descending=[False, False])
        .with_row_index("row_index")
        .with_columns(pl.col("row_index").cast(pl.Int32))
    )

    log.info(f"hits={hits.height}")

    """
    group_by feature_id:
    for each row in the group
        compute_overlap with the remaining rows in the group
        if overlap > 0.5
            keep the row with the lower evalue
        else
            keep the row
    """
    drop = []

    dropped_by_thresholds = (
        hits.filter(
            (pl.col("evalue") > max_evalue)
            & (pl.col("domain_cov") < min_domain_coverage)
        )
        .with_columns(
            pl.lit(f"'evalue > {max_evalue} or aa_cov < {min_domain_coverage}'").alias(
                "reason"
            )
        )
        .to_dicts()
    )
    drop.extend(dropped_by_thresholds)

    filtered_hits = hits.filter(
        (pl.col("evalue") <= max_evalue) & (pl.col("domain_cov") >= min_domain_coverage)
    )
    for protein_hit in filtered_hits.partition_by("feature_id"):
        hit = protein_hit.to_dicts()
        for i, row in enumerate(hit):
            for j, row2 in enumerate(hit):
                if i == j:
                    continue
                overlap = compute_overlap(
                    (row["pfam_start"], row["pfam_end"]),
                    (row2["pfam_start"], row2["pfam_end"]),
                )
                if overlap > min_overlap:
                    # pick the one with the lower evalue
                    if row["evalue"] < row2["evalue"]:
                        row["reason"] = f"'Overlap > {min_overlap}%; lower evalue (R2)'"
                        drop.append(row2)
                    else:
                        row["reason"] = f"'Overlap > {min_overlap}%; lower evalue (R1)'"
                        drop.append(row)

    drop_candidates = (
        pl.DataFrame(drop)
        .with_columns(pl.col("row_index").cast(pl.Int32))
        .unique()
        .sort(["feature_id", "pfam_start"], descending=[False, False])
    )
    log.info(f"drop_candidates={drop_candidates.height}")

    hits.write_csv(processed_output)
    drop_candidates.write_csv(drop_candidates_output)

    final_df = (
        hits.join(drop_candidates, on="row_index", how="anti")
        .rename({"pfam_id": "pfam_id_version"})
        .with_columns(
            pl.col("pfam_id_version")
            .str.split(".")
            .list.first()
            .str.replace("PF", "pfam")
            .alias("pfam_id"),
        )
    )
    log.info(f"final_df={final_df.height}")

    log.info(
        f"Orignal hmmsearch output annotated {hits.get_column('feature_id').n_unique()} features"
    )
    log.info(
        f"Filtering marked {drop_candidates.get_column('feature_id').n_unique()} features"
    )
    log.info(
        f"Finally {final_df.get_column('feature_id').n_unique()} features were sucessfully annotated"
    )
    all_features = hits.get_column("feature_id").unique()
    final_features = final_df.get_column("feature_id").unique()
    missing_features = all_features.filter(~all_features.is_in(final_features))
    log.info(f"Missing features: {missing_features.to_list()}")

    final_df.write_csv(final_output)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    main()
