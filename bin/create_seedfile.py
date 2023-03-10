#!/usr/bin/env python3
"""Generate seedfile for nf-bakta pipeline
USAGE: python create_seedfile.py -g S3path_to_genome_folder -e extension_of_files -project Project_Name -prefix Prefix_Name -s seedfile_name
EXAMPLE: python create_seedfile.py -g s3://nextflow-pipelines/nf-bakta/test/data -e .fna -project 00_TEST -prefix 20221201 -s test_20221201.seedfile.csv

# pip install -U cloudpathlib[s3] pandas
"""
import argparse
import logging
from datetime import date

import pandas as pd
from cloudpathlib import AnyPath, CloudPath

PIPELINE_RESULTS_BASE_S3PATH = "s3://genomics-workflow-core/Results/Bakta"


def usage():
    parser = argparse.ArgumentParser(
        description="Create submission seedfile for the Bakta pipeline"
    )
    parser.add_argument(
        "-g",
        "--genome-dir",
        type=str,
        required=True,
        help="Directory with all the genomes with extension '.fna'",
    )
    parser.add_argument(
        "-project",
        type=str,
        required=False,
        default="00_TEST",
        help=(
            """
    Name of the project that this analysis belongs to.
    This will become part of your output path.
    No spaces, please.
    """
        ),
    )
    parser.add_argument(
        "-prefix",
        type=str,
        required=False,
        default="00_test_prefix",
        help=(
            """
        Name of the subset of the data within this Project that this analysis belongs to.
        This will become part of your output path.
        No spaces, please.
        """
        ),
    )
    parser.add_argument(
        "-e",
        "--extension",
        type=str,
        default=".fna",
        help="extension of the genome files. Please include the '.'",
    )
    parser.add_argument(
        "-s",
        "--seedfile",
        type=str,
        required=False,
        default="",
        help="Name of the seedfile for the pipeline. By default, it is a combination of <PROJECT> and <PREFIX>",
    )
    parser.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        required=False,
        default=False,
        help="Search the S3 location recursively. [This might take more time]",
    )
    return parser.parse_args()


def generate_seedfile(
    genomes_dir: str, extension: str, recursive: bool
) -> pd.DataFrame:
    genomes = AnyPath(genomes_dir)
    if recursive:
        return pd.DataFrame(
            [
                {
                    "genome_id": str(genome.name)
                    .replace(f"{extension}", "")
                    .rstrip("."),
                    "genome_path": str(genome),
                }
                for genome in genomes.rglob(f"*{extension}")
            ]
        )
    else:
        return pd.DataFrame(
            [
                {
                    "genome_id": str(genome.name)
                    .replace(f"{extension}", "")
                    .rstrip("."),
                    "genome_path": str(genome),
                }
                for genome in genomes.glob(f"*{extension}")
            ]
        )


def upload_to_s3(
    df: pd.DataFrame, project: str, prefix: str, seedfile: str = "", s3base: str = ""
) -> CloudPath:
    if not s3base:
        s3base = PIPELINE_RESULTS_BASE_S3PATH

    if not seedfile:
        timestamp = date.today().isoformat().replace("-", "")
        seedfile = f"{timestamp}-{project}.seedfile.csv"

    s3base_obj = CloudPath(s3base)

    seedfile_s3path = s3base_obj / project / prefix / "00_seedfile" / seedfile  # type: ignore
    df.to_csv(seedfile_s3path.as_uri(), index=False)

    return seedfile_s3path


def main():
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        encoding="utf-8",
        level=logging.INFO,
    )
    args = usage()

    seedfile_df = generate_seedfile(args.genome_dir, args.extension, args.recursive)
    num_rows, _ = seedfile_df.shape
    logging.info(f"Seedfile contains paths to {num_rows} genomes")

    seedfile_s3path = upload_to_s3(
        seedfile_df, args.project, args.prefix, args.seedfile
    )
    logging.info(f"Seedfile has been uploaded here: {seedfile_s3path}")

    logging.info(
        f"""Submission command:\n
aws batch submit-job \\
    --job-name nf-bakta-{args.prefix} \\
    --job-queue priority-maf-pipelines \\
    --job-definition nextflow-production \\
    --container-overrides command=FischbachLab/nf-bakta,\\
"--seedfile","{seedfile_s3path}",\\
"--project","{args.project}",\\
"--prefix","{args.prefix}"

*NOTE*\n\tThe above command assumes that you'll be running the code from the `main` branch.
    If that is not the case, please add `"-r","branch-name"` to the submission command as well.
    If you don't know what that means, this probably doesn't apply to you, please feel free to 
    ignore this message.
    """
    )


if __name__ == "__main__":
    main()
