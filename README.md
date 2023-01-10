# NF-BAKTA

## Seedfile format

- comma delimited
- two columns with headers: `genome_id`, `genome_path`
- Example: [test_20221220_0.seedfile.csv](test/test_20221220_0.seedfile.csv)

The helper script, [create_seedfile.py](bin/create_seedfile.py), will create the properly formatted seedfile for you if you can point it to an S3 path.

```bash
cd nf-bakta
python bin/create_seedfile.py \
    -g s3://maf-users/Nathan_Johns/DBs/Segata_Genomes/Fastas/ \
    -project UHGG_Annotation \
    -prefix 20221221 \
    --extension .fasta
```

This helper script will also recommend a job submission command that you can use to launch your job using the seedfile that was just created.

## Test

```bash
aws batch submit-job \
    --job-name nf-bakta-test-1 \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-bakta,\
"--seedfile","s3://genomics-workflow-core/Results/Bakta/00_Test/seedfiles/test_20221220_1.seedfile.csv",\
"--project","00_Test",\
"--prefix","20221220-2"
```

## Bakta Database

The database for this pipeline is stored on our EFS at `/mnt/efs/databases/Bakta/db/v4.0`. This path is provided as the `bakta_db` parameter. Note that this path should not be staged within the pipleine, but just passed as a value. This is done because all containers have access to that path, i.e. it's already available/staged/mounted for the container to use.
