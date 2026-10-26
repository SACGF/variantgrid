#!/bin/bash
# Set up the headless browser harness in $VG_BROWSER_HOME (default ~/.cache/claude-browser): a uv venv with
# playwright, the chromium headless shell, and - on a box without root - the shared libraries chromium needs
# that aren't installed, extracted from Ubuntu packages with apt-get download (no root needed).
# See scripts/browser/README.md.
set -euo pipefail

HOME_DIR="${VG_BROWSER_HOME:-$HOME/.cache/claude-browser}"
PLAYWRIGHT_VERSION="1.64.0"
# Ubuntu noble package names for every library chromium-headless-shell 156 wanted on vg-test2
SIDELOAD_PACKAGES="libatk1.0-0t64 libatk-bridge2.0-0t64 libatspi2.0-0t64 libxcomposite1 libxdamage1 libxfixes3 libxrandr2 libxi6 libgbm1 libasound2t64"

mkdir -p "$HOME_DIR"
cd "$HOME_DIR"

if [ ! -x pw/bin/python ]; then
    uv venv pw
fi
uv pip install --python pw/bin/python "playwright==${PLAYWRIGHT_VERSION}"
pw/bin/python -m playwright install chromium-headless-shell

shell=$(ls -d ~/.cache/ms-playwright/chromium_headless_shell-*/chrome-headless-shell-linux64/chrome-headless-shell | tail -1)
missing() { LD_LIBRARY_PATH="$HOME_DIR/libs/usr/lib/x86_64-linux-gnu" ldd "$shell" | grep "not found" || true; }

if [ -n "$(missing)" ]; then
    echo "chromium is missing shared libraries; side-loading them into $HOME_DIR/libs"
    mkdir -p debs libs
    (cd debs && apt-get download $SIDELOAD_PACKAGES)
    for deb in debs/*.deb; do dpkg -x "$deb" libs; done
    rm -rf debs
fi

if [ -n "$(missing)" ]; then
    echo "still missing after side-loading (add the package to SIDELOAD_PACKAGES):"
    missing
    exit 1
fi

LD_LIBRARY_PATH="$HOME_DIR/libs/usr/lib/x86_64-linux-gnu" pw/bin/python - <<'EOF'
from playwright.sync_api import sync_playwright
with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page()
    page.set_content("<html><body>ok</body></html>")
    print("browser ok:", browser.version, page.evaluate("document.compatMode"))
    browser.close()
EOF
echo "harness ready in $HOME_DIR - next: scripts/browser/README.md 'Each session'"
