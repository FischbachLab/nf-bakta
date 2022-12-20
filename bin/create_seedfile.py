#!/usr/bin/env python3
"""Generate seedfile for nf-binqc pipeline
USAGE: python create_seedfile.py S3path_to_genome_folder extension_of_files local_seedfile_output
EXAMPLE: python create_seedfile.py s3://nextflow-pipelines/nf-binqc/test/data fna ../test/test_20221201.seedfile.csv

# pip install -U cloudpathlib[s3] pandas
"""
from cloudpathlib import S3Path
import pandas as pd
import sys


def main():
    s3path = S3Path(sys.argv[1])
    extension = sys.argv[2]
    seedfile = sys.argv[3]

    pd.DataFrame(
        [
            {
                "genome_id": str(genome.name).replace(f".{extension}", ""),
                "genome_path": str(genome),
            }
            for genome in s3path.glob(f"*.{extension}")
        ]
    ).to_csv(seedfile, index=False)


if __name__ == "__main__":
    main()
