# Some useful oneliners

## List all jobs in the `RUNNING` state

```bash
aws batch list-jobs --no-paginate --job-queue default-maf-pipelines --job-status "RUNNING" --output text > jobs.tsv
```

Other options for `--job-status` are: `SUBMITTED`, `RUNNABLE` and `FINISHED`.

## Terminate all jobs from the `list-jobs` ouptut

```bash
cut -f 4 jobs.tsv | parallel -j 4 "aws batch terminate-job  --job-id {} --reason 'Terminated by maf-user'"
```
