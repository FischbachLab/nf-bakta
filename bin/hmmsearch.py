#!/usr/bin/env python3
# Description: This script is used to search for Pfam domains in the predicted proteins of the CDSs
# and filter them such that only the best hit is retained for each region on a protein.
import argparse
import logging
from collections import OrderedDict, defaultdict
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
    compute_overlap(segment1, segment2))
    >>> 100.0
    """
    # 0  10 20 30 40
    # |==|==|==|==| segment1
    #       |==|    segment2

    # => 10/10 * 100 = 100.0

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
    percentage_overlap = overlap_length / length_shorter_segment

    return percentage_overlap


def usage():
    parser = argparse.ArgumentParser(description="Pfam domain search")
    parser.add_argument(
        "--faa",
        help="Input protein file in fasta format",
        required=True,
    )
    parser.add_argument(
        "--genome_id",
        help="genome_id, to be used as output file prefix",
        required=True,
    )
    parser.add_argument(
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
    # parser.add_argument(
    #     "--max_evalue",
    #     help="Maximum evalue threshold",
    #     default=1e-5,
    #     type=float,
    # )
    parser.add_argument(
        "--min_domain_coverage",
        help="Minimum domain coverage threshold [0-1]",
        default=0.4,
        type=float,
    )
    parser.add_argument(
        "--min_overlap",
        help="Minimum overlap threshold [0-1]",
        default=0.5,
        type=float,
    )
    return parser.parse_args()


def main():
    args = usage()
    protein_file = args.input
    output_prefix = args.output
    pfam_db = args.pfam_db

    threads = args.threads
    # max_evalue = args.max_evalue
    min_domain_coverage = args.min_domain_coverage
    min_overlap = args.min_overlap

    ## TEST DATA
    # protein_file = "test_data/proteins/Clostridium-sporogenes-ATCC-15579-MAF-2.faa"
    # output_prefix = "test_data/proteins/Clostridium-sporogenes-ATCC-15579-MAF-2"
    # pfam_db = "/mnt/efs/databases/HMMs/pfam/v36.0/Pfam-A.hmm"

    # threads = 12
    # max_evalue = 1e-5 ## no longer supported
    # min_domain_coverage = 0.4
    # min_overlap = 50

    assert 0 <= min_domain_coverage <= 1, "min_domain_coverage must be between 0 and 1"
    assert 0 <= min_overlap <= 1, "min_overlap must be between 0 and 1"

    domtbl_output = f"{output_prefix}.hmmsearch_domtblout.tsv"
    processed_output = f"{output_prefix}.hmmsearch_all_hits.csv"
    drop_candidates_output = f"{output_prefix}.hmmsearch_dropped_hits.csv"
    stats_output = f"{output_prefix}.hmmsearch_stats.tsv"
    final_output = f"{output_prefix}.hmmsearch_tbl.csv"

    stats = defaultdict(str)

    hits = pl.DataFrame(
        run_hmmsearch(protein_file, pfam_db, domtbl_output, threads)
    ).sort(["feature_id", "pfam_start"], descending=[False, False])
    hits.write_csv(processed_output)
    log.debug(f"hits={hits.height}")
    stats["num_rows_original_output"] = hits.height

    """
    group_by feature_id:
    for each row in the group
        compute_overlap with the remaining rows in the group
        if overlap > 0.5
            keep the row with the lower evalue
    """
    drop = []
    hits = (
        hits.filter(pl.col("domain_cov") >= min_domain_coverage)
        .with_row_index("row_index")
        .sort(["feature_id", "row_index"], descending=[False, False])
    )
    stats["num_rows_after_dom_cov_filter"] = hits.height

    for protein_hit in hits.partition_by("feature_id"):
        hit = protein_hit.to_dicts()
        for i, row in enumerate(hit):
            for j, row2 in enumerate(hit):
                if j <= i:
                    continue
                overlap = compute_overlap(
                    (row["pfam_start"], row["pfam_end"]),
                    (row2["pfam_start"], row2["pfam_end"]),
                )
                if overlap > min_overlap:
                    # pick the one with the lower evalue
                    if row["evalue"] < row2["evalue"]:
                        row2["reason"] = f"Overlaps: {row['pfam_id']}"
                        drop.append(row2)
                    else:
                        row["reason"] = f"Overlaps: {row2['pfam_id']}"
                        drop.append(row)

    drop_candidates = pl.DataFrame(drop).unique()
    drop_rows = drop_candidates.get_column("row_index").unique()
    log.debug(f"drop_candidates={len(drop_rows)}")
    stats["num_rows_dropped"] = len(drop_rows)
    drop_candidates.write_csv(drop_candidates_output)

    final_df = (
        hits.filter(~pl.col("row_index").is_in(drop_rows))
        .rename({"pfam_id": "pfam_id_version"})
        .with_columns(
            pl.col("pfam_id_version")
            .str.split(".")
            .list.first()
            .str.replace("PF", "pfam")
            .alias("pfam_id"),
        )
    )
    final_df.write_csv(final_output)
    log.debug(f"final_df={final_df.height}")
    stats["num_rows_final"] = final_df.height

    stats["original_features_w_pfams"] = hits.get_column("feature_id").n_unique()
    stats["dropped_features_w_pfams"] = drop_candidates.get_column(
        "feature_id"
    ).n_unique()
    stats["final_features_w_pfams"] = final_df.get_column("feature_id").n_unique()

    ## Did we completely miss/remove any features?
    dropped_features = drop_candidates.get_column("feature_id").unique()
    all_viable_features = (
        hits.filter(~pl.col("feature_id").is_in(dropped_features))
        .get_column("feature_id")
        .unique()
    )
    final_features = final_df.get_column("feature_id").unique()
    missing_features = all_viable_features.filter(
        ~all_viable_features.is_in(final_features)
    )

    stats["all_viable_features"] = all_viable_features.n_unique()
    stats["missing_features_num"] = missing_features.n_unique()
    stats["missing_features"] = ",".join(missing_features.to_list())

    log.info(f"Stats: {stats}")
    # save stats to file
    with open(stats_output, "w") as f:
        for key, value in stats.items():
            f.write(f"{key}\t{value}\n")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    main()
