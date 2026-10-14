#!/usr/bin/env bash
# nf-core/mag, cobinning mode. Build the samplesheets first (see README.md).
set -euo pipefail
cd "$(dirname "$0")"
export NXF_SYNTAX_PARSER=v1   # nf-core/mag 3.3.1 config predates the Nextflow >= 25.10 strict parser

nextflow run nf-core/mag -r 3.3.1 \
    -profile apptainer \
    -params-file params/cobinning.json \
    -resume "$@"
