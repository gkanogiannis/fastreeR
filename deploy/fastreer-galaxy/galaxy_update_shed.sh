#!/usr/bin/env bash
set -euo pipefail

tools_root="$1"
SHED="$2"
KEY="$3"

planemo \
    shed_update \
    --force_repository_creation\
    --recursive\
    --shed_target $SHED \
    --shed_key $KEY \
    ${tools_root}

#tools='DIST2TREE FASTA2DIST VCF2DIST VCF2TREE'
#for t in $tools; do
#    t=$(echo $t | tr '[:upper:]' '[:lower:]')
#    planemo shed_update \
#            --force_repository_creation \
#            --shed_target $SHED \
#            --shed_key $KEY \
#            ${tools_root}/$t
#done
