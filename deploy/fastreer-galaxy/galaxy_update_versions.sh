#!/usr/bin/env bash
set -euo pipefail

tools_root="$1"
VERSION="$2"

echo "Updating Galaxy tool XMLs to fastreer version: ${VERSION}"

# Update version in <tool ... version="...">
find ${tools_root}/tools -name ".xml" -print0 | while IFS= read -r -d '' file; do
    echo "Updating tool version in: $file"
    sed -i -E "s|(<token name=\"@TOOL_VERSION@\">)[^<]+(</token>)|\1${VERSION}\2|" "$file"
done

echo "All Galaxy tool versions updated to ${VERSION}"
