# NF-BAKTA

## Seedfile format

- comma delimited
- two columns with headers: `genome_id`, `genome_path`
- Example: [test_20221220_0.seedfile.csv](test/test_20221220_0.seedfile.csv)

The helper script, [create_seedfile.py](bin/create_seedfile.py), will create the properly formatted seedfile for you if you can point it to an S3 path.

```bash
cd nf-bakta
python bin/create_seedfile.py s3://nextflow-pipelines/nf-bakta/test/data fna test/test_20221220_0.seedfile.csv
```

### OR

```bash
cd nf-bakta
python bin/create_seedfile.py s3://maf-users/Nathan_Johns/Scratch/ fna test/test_20221220.seedfile.csv
```

## Upload the seedfile to S3

```bash
aws s3 cp test/test_20221220_1.seedfile.csv s3://genomics-workflow-core/Results/Bakta/00_Test/seedfiles/test_20221220_1.seedfile.csv
```

## Test

```bash
aws batch submit-job \
    --job-name nf-bakta-test-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-bakta,\
"--seedfile","s3://genomics-workflow-core/Results/Bakta/00_Test/seedfiles/test_20221220_1.seedfile.csv",\
"--project","00_Test",\
"--prefix","20221220-1"
```
