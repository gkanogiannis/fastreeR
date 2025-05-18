#!/usr/bin/env bash
set -euo pipefail

KEY="$1"
SHED="$2"
dir=`dirname $0`

planemo shed_update \
            --force_repository_creation \
            --shed_target $SHED \
            --shed_key $KEY \
            ${dir}
