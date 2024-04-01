#!/bin/bash -x

set -euo pipefail

num_jobs=10000

while (( num_jobs > 10 )); do
    aws batch list-jobs --no-paginate --job-queue default-maf-pipelines --job-status "RUNNABLE" --output text > jobs.tsv
    num_jobs=$(wc -l jobs.tsv | awk '{print $1}')
    tail -n +2 jobs.tsv | cut -f 4 | parallel -j 3 "aws batch terminate-job  --job-id {} --reason 'Terminated by SJ'"
done

aws batch list-jobs --no-paginate --job-queue default-maf-pipelines --job-status "RUNNING" --output text > running_jobs.tsv
cut -f 4 running_jobs.tsv | parallel -j 3 "aws batch terminate-job  --job-id {} --reason 'Terminated by SJ'"