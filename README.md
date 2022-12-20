# NF-BINQC

Requires that all the genomes in the `--fastas` location be in fasta format, with the same extension. By default, the pipeline assumes that they all have the extension `fna`. If this is not the case, please use the flag `--ext` to reflect the appropriate extension.

## Seedfile format

- comma delimited
- two columns with headers: `genome_id`, `genome_path`
- Example: [test_20221201.seedfile.csv](test/test_20221201.seedfile.csv)

The helper script, [create_seedfile.py](bin/create_seedfile.py), will create the properly formatted seedfile for you if you can point it to an S3 path.

```bash
cd nf-binqc/bin
python create_seedfile.py s3://nextflow-pipelines/nf-binqc/test/data fna ../test/test_20221201.seedfile.csv
```

## Test

```bash
aws batch submit-job \
    --job-name nf-binqc-test \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-binqc,\
"-r","sj-sdlc",\
"--seedfile","s3://nextflow-pipelines/nf-binqc/test/test_20221201.seedfile.csv",\
"--project","00_Test",\
"--prefix","20221205"
```

## Actual sample command

```bash
aws batch submit-job \
    --job-name nf-binqc-SCv2_4_20210212 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-binqc,\
"--fastas","s3://maf-versioned/ninjamap/Index/SCv2_4_20210212/fasta",\
"--project","SCv2_4_20210212"
```

or locally

```bash
nextflow run . --fastas test/data --project 00_TEST --outdir test/result/
```

## Update

Note that updating the code here will *not* update the pipeline automatically.

```{bash}
cd nf-binqc
aws s3 sync . s3://nextflow-pipelines/nf-binqc --exclude ".git/*" --exclude "test/result/*" --delete --profile maf
```
