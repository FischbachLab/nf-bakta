# NF-BAKTA

## Python Environment

```bash
conda create -n bakta -c conda-forge python=3.11 cloudpathlib-s3 pandas notebook fsspec s3fs=2023.3.0
```

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

### Local

```bash
nextflow run main.nf \
    --seedfile test/local.seedfile.csv \
    --project 00_Test \
    --prefix 20240306-pfam
```

### Remote

```bash
aws batch submit-job \
    --job-name nf-bakta-pfam-test \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command=FischbachLab/nf-bakta,\
"-r","pfam-annotation",\
"--seedfile","s3://genomics-workflow-core/Results/Bakta/00_Test/seedfiles/test_20221220_1.seedfile.csv",\
"--project","00_Test",\
"--prefix","20240401-pfam"
```

## Bakta Database

```text
v4.0 = ???
v5.0, type=full, 2023-02-20, DOI: 10.5281/zenodo.7669534
v5.1, type=full, 2024-01-19, DOI: 10.5281/zenodo.10522951
```

The database for this pipeline is stored on our EFS at `/mnt/efs/databases/Bakta/db/v5.0`. This path is provided as the `bakta_db` parameter. Note that this path should not be staged within the pipleine, but just passed as a value. This is done because all containers have access to that path, i.e. it's already available/staged/mounted for the container to use.

### Download new database

This was needed when Bakta moved from db schema v4.0 to v5.0.

```bash
cd /mnt/efs/databases/Bakta/db
mkdir v5
docker container run \
    --rm \
    -u $(id -u):$(id -g) \
    -v /mnt/efs/databases/Bakta/db/v5:/db \
    458432034220.dkr.ecr.us-west-2.amazonaws.com/bakta:1.9.3 \
    bakta_db download --output /db --type full
```

### Update existing database

```bash
mkdir -p /mnt/efs/databases/Bakta/db/tmp
cd /mnt/efs/databases/Bakta/db
docker container run \
-it \
    --rm \
    -v /mnt/efs/databases/Bakta/db/v5:/db \
    -v /mnt/efs/databases/Bakta/db/db_tmp:/bakta_tmp \
    458432034220.dkr.ecr.us-west-2.amazonaws.com/bakta:1.9.3 \
    bakta_db update --db /db --tmp-dir /tmp
```

## PFam Annotation Parsing

References:

- `hmmscan` vs `hmmsearch`: [Use hmmsearch](http://cryptogenomicon.org/hmmscan-vs-hmmsearch-speed-the-numerology.html) to annotate the proteins from the bakta pipeline. 

- `pyhmmer`:
    - [repo](https://github.com/althonos/pyhmmer)
    - [Read the docs](https://pyhmmer.readthedocs.io/en/stable/)
    - [Example of usage in Bakta](https://github.com/oschwengers/bakta/blob/d6443639958750c3bece5822e84978271d1a4dc7/bakta/features/cds.py#L106)

- [Parsing `HMMsearch` output](https://dark.metagenomics.eu/pfam-annotation)