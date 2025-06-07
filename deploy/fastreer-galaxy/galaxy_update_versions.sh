#!/usr/bin/env bash
set -euo pipefail

tools_root="$1"
VERSION="$2"

echo "Updating Galaxy tool XMLs to fastreer version: ${VERSION}"

# Update version in <tool ... version="...">
find ${tools_root}/tools -name "*.xml" -print0 | while IFS= read -r -d '' file; do
    echo "Updating tool version in: $file"
    sed -i -E "s/(<tool[^>]*version=\")([^\"]+)(\"[^>]*>)/\1${VERSION}\3/" "$file"
done

# Update fastreer package requirement version
find ${tools_root}/tools -name "*.xml" -print0 | while IFS= read -r -d '' file; do
    echo "Updating Conda requirement in: $file"
    sed -i -E "s|(<requirement[^>]*>fastreer)([^<]*)(</requirement>)|\1=${VERSION}\3|g" "$file"
done

echo "All Galaxy tool versions updated to ${VERSION}"
