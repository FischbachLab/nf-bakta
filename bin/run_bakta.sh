#!/bin/bash -x

set -euo pipefail

assembly_id=${1}
assembly=${2}
cpus=${3}
bakta_db=${4}

min_seqs=1
num_seqs=$(grep -c '>' "${assembly}")

if ((num_seqs >= min_seqs)); then
    bakta \
    --db "${bakta_db}" \
    --verbose \
    --prefix "${assembly_id}" \
    --compliant \
    --threads "${cpus}" \
    --skip-plot \
    "${assembly}"

    sha256sum "${assembly}" &> "${assembly_id}.sha256"
else
    echo "${assembly} contains ${num_seqs} sequences. Cannot process further"
fi

exit 0