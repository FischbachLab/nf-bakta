#!/bin/bash -x

set -euo pipefail

assembly_id=${1}
assembly=${2}
cpus=${3}
bakta_db=${4}

# check if the baktadb directory is accessible
if [[ ! -d "${bakta_db}" ]]; then
    echo "Error: bakta_db directory not found: ${bakta_db}"
    exit 1
fi

min_seqs=1
num_seqs=$(grep -c '>' "${assembly}" || true)

if (( num_seqs >= min_seqs )); then
    bakta \
    --db "${bakta_db}" \
    --verbose \
    --prefix "${assembly_id}" \
    --compliant \
    --threads "${cpus}" \
    "${assembly}"
    # --skip-plot \

    sha256sum "${assembly}" &> "${assembly_id}.sha256"
else
    echo "${assembly} contains ${num_seqs} sequences. Cannot process further"
fi

exit 0