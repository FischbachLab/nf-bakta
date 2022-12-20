# NF-BAKTA

## Seedfile format

- comma delimited
- two columns with headers: `genome_id`, `genome_path`
- Example: [test_20221201.seedfile.csv](test/test_20221201.seedfile.csv)

The helper script, [create_seedfile.py](bin/create_seedfile.py), will create the properly formatted seedfile for you if you can point it to an S3 path.

```bash
cd nf-binqc
python bin/create_seedfile.py s3://nextflow-pipelines/nf-binqc/test/data fna test/test_20221220_0.seedfile.csv
```

## Test

```bash
aws batch submit-job \
    --job-name nf-bakta-test-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-bakta,\
"--seedfile","s3://nextflow-pipelines/nf-bakta/00_Test/seedfiles/test_20221220_1.seedfile.csv",\
"--project","00_Test",\
"--prefix","20221220-1"
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
