"""
Evaluate one JS expression on a page in quirks and standards mode and print both results. An `async () =>`
is fine, so a tab or a node can be clicked and awaited before measuring. See scripts/browser/README.md.

    $VG_BROWSER_HOME/pw/bin/python scripts/browser/jsq.py [--wait MS] [--modes quirks,standards] URL 'JS'
"""
import argparse
import json

from harness import BASE, launch, new_context
from playwright.sync_api import sync_playwright


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("url")
    ap.add_argument("js")
    ap.add_argument("--wait", type=int, default=2500, help="ms after networkidle before evaluating")
    ap.add_argument("--modes", default="quirks,standards")
    args = ap.parse_args()
    with sync_playwright() as p:
        browser = launch(p)
        for mode in args.modes.split(","):
            ctx = new_context(browser, mode)
            page = ctx.new_page()
            page.goto(BASE + args.url, wait_until="networkidle", timeout=90000)
            page.wait_for_timeout(args.wait)
            result = page.evaluate(args.js)
            print(f"== {mode} ({page.evaluate('document.compatMode')})")
            print(result if isinstance(result, str) else json.dumps(result, indent=1))
            ctx.close()
        browser.close()


if __name__ == "__main__":
    main()
