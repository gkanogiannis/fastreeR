#!/usr/bin/env bash
set -euo pipefail

VERSION="$1"
dir=`dirname $0`

echo "🔄 Updating Galaxy tool XMLs to fastreer version: ${VERSION}"

# Update version in <tool ... version="...">
find ${dir}/tools -name "*.xml" -print0 | while IFS= read -r -d '' file; do
    echo "📄 Updating tool version in: $file"
    sed -i -E "s/(<tool[^>]*version=\")([^\"]+)(\"[^>]*>)/\1${VERSION}\3/" "$file"
done

# Update fastreer package requirement version
find ${dir}/tools -name "*.xml" -print0 | while IFS= read -r -d '' file; do
    echo "📦 Updating Conda requirement in: $file"
    sed -i -E "s|(<requirement[^>]*>fastreer)([ <>=0-9\.,]*)(</requirement>)|\1=${VERSION}\3|g" "$file"
done

echo "✅ All Galaxy tool versions updated to ${VERSION}"
