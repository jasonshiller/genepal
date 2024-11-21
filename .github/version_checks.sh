#!/usr/bin/env bash

set -euo pipefail

config_version=$(sed -n "/^\s*version\s*=\s*'/s/version//p" nextflow.config | tr -d "=[:space:]'")
cff_version=$(sed -n '/^version: /s/version: //p' CITATION.cff | tr -d '[:space:]')

if [[ $config_version != $cff_version ]]; then
    echo 'config_version != cff_version'
    exit 1
fi

# Check CHANGELOG version
head -10 CHANGELOG.md | grep "## v$config_version - " >/dev/null \
    || (echo 'Failed to match CHANGELOG version'; exit 1)

# Check .nf-core.yml version
tail -5 .nf-core.yml | grep "version: $config_version" >/dev/null \
    || (echo 'Failed to match .nf-core.yml version'; exit 1)

# Check multiqc_config.yml version
config_version_plain=$(echo $config_version | sed 's/.*\dev$/dev/; t; s/\+.*//')
head -n 5 assets/multiqc_config.yml | grep "https://github.com/plant-food-research-open/genepal/blob/$config_version_plain/docs/usage.md" >/dev/null \
    || (echo 'Failed to match multiqc_config.yml version'; exit 1)
