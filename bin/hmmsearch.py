# Description: This script is used to search for Pfam domains in the predicted proteins of the CDSs.
import logging
from collections import OrderedDict
from typing import Optional, Sequence

import polars as pl
import pyhmmer
from Bio import SeqIO
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
    pfam_length: Optional[int | None] = None
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


def get_protein_metadata(proteins_fasta: str) -> dict:
    """
    Use BioPython to read the protein fasta file and return the metadata
    """
    return {
        record.id: {
            "description": record.description,
            "name": record.name,
            "length": len(record.seq),
        }
        for record in SeqIO.parse(proteins_fasta, "fasta")
    }


def run_hmmsearch(
    proteins_fasta: str, pfam_db_path: str, output_file: str, threads: int
) -> Sequence[dict]:
    """Detect Pfam-A entries"""
    pfam_hits = []
    protein_metadata = get_protein_metadata(proteins_fasta)
    proteins = read_protein_fasta(proteins_fasta)
    # cds = OrderedDict()
    with pyhmmer.plan7.HMMFile(pfam_db_path) as hmm, open(output_file, "wb") as fh:
        for top_hits in pyhmmer.hmmsearch(
            hmm, proteins, bit_cutoffs="gathering", cpus=threads
        ):
            for hit in top_hits:
                aa_identifier = hit.name.decode()
                cds_len = protein_metadata[aa_identifier]["length"]
                domain_cov = (
                    hit.best_domain.alignment.hmm_to
                    - hit.best_domain.alignment.hmm_from
                    + 1
                ) / len(hit.best_domain.alignment.hmm_sequence)
                aa_cov = (
                    hit.best_domain.alignment.target_to
                    - hit.best_domain.alignment.target_from
                    + 1
                    # ) / len(cds["aa"])
                ) / cds_len

                pfam = OrderedDict()
                pfam["pfam_id"] = hit.best_domain.alignment.hmm_accession.decode()
                pfam["feature_id"] = hit.name.decode()
                pfam["start"] = hit.best_domain.alignment.target_from
                pfam["stop"] = hit.best_domain.alignment.target_to
                pfam["pfam_length"] = len(hit.best_domain.alignment.hmm_sequence)
                pfam["evalue"] = hit.evalue
                pfam["score"] = hit.score
                pfam["name"] = hit.best_domain.alignment.hmm_name.decode()
                pfam["aa_cov"] = aa_cov
                pfam["domain_cov"] = domain_cov

                pfam_hits.append(Pfam(**pfam).to_dict())

            top_hits.write(fh, format="targets", header=True)
    return pfam_hits


def main():
    protein_file = "test_data/proteins/Clostridium-sporogenes-ATCC-15579-MAF-2.faa"
    output_file = (
        "test_data/proteins/Clostridium-sporogenes-ATCC-15579-MAF-2.tblout.tsv"
    )
    pfam_db = "/mnt/efs/databases/HMMs/pfam/v36.0/Pfam-A.hmm"
    threads = 12
    hits = pl.DataFrame(
        run_hmmsearch(protein_file, pfam_db, output_file, threads)
    ).sort(["feature_id", "pfam_start"], descending=[False, False])

    log.info(f"hits={hits.height}")
    hits.write_csv(
        "test_data/proteins/Clostridium-sporogenes-ATCC-15579-MAF-2.filtered_hits.tsv"
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    main()
