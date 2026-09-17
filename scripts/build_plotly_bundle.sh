#!/bin/bash
# Rebuilds the vendored plotly.js bundle in variantgrid/static_files/default_static/js/lib/
#
# None of plotly's off-the-shelf partial bundles cover what we plot - sankey ships only in the
# full ~4.8MB build - so we bundle plotly's core plus the six extra traces we use, which comes
# to about a third of that. Run this to move to a new plotly version, then point the <script>
# tags at the new filename.
#
# Usage: scripts/build_plotly_bundle.sh [version]

set -euo pipefail

PLOTLY_VERSION=${1:-3.7.0}
VARIANTGRID_DIR=$(cd "$(dirname "$0")/.." && pwd)
OUT=${VARIANTGRID_DIR}/variantgrid/static_files/default_static/js/lib/plotly-custom-${PLOTLY_VERSION}.min.js

BUILD_DIR=$(mktemp -d)
trap 'rm -rf "${BUILD_DIR}"' EXIT
cd "${BUILD_DIR}"

npm init -y > /dev/null
npm install --no-audit --no-fund "plotly.js@${PLOTLY_VERSION}" esbuild > /dev/null

cat > entry.js <<'ENTRY'
const Plotly = require('plotly.js/lib/core');  // scatter is built in
Plotly.register([
    require('plotly.js/lib/bar'),
    require('plotly.js/lib/box'),
    require('plotly.js/lib/heatmap'),
    require('plotly.js/lib/histogram'),
    require('plotly.js/lib/pie'),
    require('plotly.js/lib/sankey'),
]);
module.exports = Plotly;
ENTRY

# define:global stands in for the shim browserify (plotly's own bundler) gives has-hover,
# and loader:.css=empty drops the maplibre stylesheet pulled in for map traces we don't register
npx esbuild entry.js --bundle --minify --format=iife --global-name=Plotly \
    --define:global=globalThis --loader:.css=empty --legal-comments=none --outfile=bundle.js

cat > "${OUT}" <<HEADER
/**
* plotly.js (VariantGrid custom bundle - minified) v${PLOTLY_VERSION}
* Traces: scatter, bar, box, heatmap, histogram, pie, sankey
* Copyright 2012-2026, Plotly, Inc.
* All rights reserved.
* Licensed under the MIT license
*/
HEADER
cat bundle.js >> "${OUT}"

echo "Wrote ${OUT} ($(du -h "${OUT}" | cut -f1))"
