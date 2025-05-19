#!/usr/bin/env bash
set -euo pipefail

tools_root="$1"

xml_files=$(find ${tools_root} -name '*.xml' -print0 | tr '\0' ' ')
planemo \
    test --skip_venv \
    --test_output ${tools_root}/galaxy-test-report.html \
    --test_output_json /dev/null \
    ${xml_files}
