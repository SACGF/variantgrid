"""
Load one page and report what its JavaScript complained about: uncaught exceptions (with stacks), console
errors and warnings, failed requests and 4xx/5xx responses. Prints one JSON object, including the served
document so the caller can map an inline script's line number back to its template. `vg page --js-writes-to-database` runs this;
it can also be run against a runserver as in scripts/browser/README.md.

    $VG_BROWSER_HOME/pw/bin/python scripts/browser/console.py [--wait MS] [--click SELECTOR ...] URL
"""
import argparse
import json

from harness import BASE, launch, new_context
from playwright.sync_api import Error, sync_playwright

CONSOLE_TYPES = {"error", "warning", "assert"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("url")
    ap.add_argument("--wait", type=int, default=1500, help="ms after networkidle, for grids and deferred ajax")
    ap.add_argument("--click", action="append", default=[], help="click this selector after load (a tab), in order")
    args = ap.parse_args()

    errors, console, requests, click_failures = [], [], [], []
    with sync_playwright() as p:
        browser = launch(p)
        ctx = new_context(browser, "quirks")
        page = ctx.new_page()
        page.on("pageerror", lambda e: errors.append({"message": e.message, "stack": e.stack or ""}))
        page.on("console", lambda m: console.append({"type": m.type, "text": m.text, **m.location})
                if m.type in CONSOLE_TYPES else None)
        page.on("requestfailed", lambda r: requests.append({"method": r.method, "url": r.url, "error": r.failure}))
        page.on("response", lambda r: requests.append({"method": r.request.method, "url": r.url, "status": r.status})
                if r.status >= 400 else None)

        response = page.goto(BASE + args.url, wait_until="networkidle", timeout=90000)
        page.wait_for_timeout(args.wait)
        for selector in args.click:
            try:
                page.click(selector, timeout=10000)
            except Error as e:
                click_failures.append(f"{selector}: {e.message.splitlines()[0]}")
                continue
            page.wait_for_load_state("networkidle")
            page.wait_for_timeout(args.wait)

        print(json.dumps({
            "url": page.url,
            "status": response.status if response else None,
            "document": response.text() if response else "",
            "errors": errors,
            "console": console,
            "requests": [r for r in requests if not r["url"].endswith("favicon.ico")],
            "click_failures": click_failures,
        }))
        ctx.close()
        browser.close()


if __name__ == "__main__":
    main()
