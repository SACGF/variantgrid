"""
Screenshot and measure pages served by a working-tree runserver, as claude_agent, in quirks mode (as served)
and/or standards mode (<!DOCTYPE html> prepended by a playwright route). Prints one JSON line per render with
status, compatMode, page height, the boxes of the --measure selectors and page errors; saves PNGs and
results.json under --out. See scripts/browser/README.md.

    $VG_BROWSER_HOME/pw/bin/python scripts/browser/shot.py [options] URL...
"""
import argparse
import json
import re
import sys
from pathlib import Path

from harness import BASE, launch, new_context
from playwright.sync_api import sync_playwright

# tag#id.class:nth for every element, keyed so the two modes can be matched by compare.py
DUMP_JS = """() => {
  const out = {};
  const walk = (el, prefix) => {
    let i = 0;
    for (const c of el.children) {
      i++;
      if (c.tagName === 'SCRIPT' || c.tagName === 'STYLE' || c.tagName === 'LINK') continue;
      const id = c.id ? '#' + c.id.replace(/\\d+/g, 'N') : '';
      const cls = typeof c.className === 'string' && c.className ? '.' + c.className.trim().split(/\\s+/).slice(0, 3).join('.') : '';
      const key = prefix + '/' + c.tagName.toLowerCase() + id + cls + ':' + i;
      const r = c.getBoundingClientRect();
      const cs = getComputedStyle(c);
      out[key] = [Math.round(r.width), Math.round(r.height), Math.round(r.top + window.scrollY), cs.display, cs.fontSize];
      walk(c, key);
    }
  };
  walk(document.body, '');
  return out;
}"""

MEASURE_JS = """(sel) => Array.from(document.querySelectorAll(sel)).slice(0, 5).map(e => {
    const r = e.getBoundingClientRect();
    return {w: Math.round(r.width), h: Math.round(r.height), top: Math.round(r.top + window.scrollY),
            visible: !!(e.offsetWidth || e.offsetHeight || e.getClientRects().length)};
})"""

CLIP_JS = """(sel) => { const e = document.querySelector(sel); if (!e) return null; const r = e.getBoundingClientRect();
    return {x: 0, y: Math.max(0, r.top + window.scrollY - 120), width: window.innerWidth, height: r.height + 240}; }"""

# the base layout scrolls .main-content, not the document: opened up for a whole-page screenshot only
UNROLL_CSS = "body, .main-content { overflow: visible !important; height: auto !important; max-height: none !important; }"


def slug(url: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", url.strip("/")).strip("_") or "root"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("urls", nargs="+", help="paths on the server, e.g. /variantopedia/dashboard")
    ap.add_argument("--modes", default="quirks,standards")
    ap.add_argument("--vw", type=int, default=1280)
    ap.add_argument("--vh", type=int, default=900)
    ap.add_argument("--wait", type=int, default=1500, help="ms after networkidle, for plotly and grids")
    ap.add_argument("--measure", default="", help="comma-separated selectors whose boxes to report")
    ap.add_argument("--click", default="", help="click this selector after load (open a tab) before measuring")
    ap.add_argument("--patch", action="append", default=[], help="OLD=>NEW replacement applied to served HTML")
    ap.add_argument("--dump", action="store_true", help="record every element's box for compare.py")
    ap.add_argument("--unroll", action="store_true", help="open the .main-content scroller so the screenshot is the whole page")
    ap.add_argument("--clip", default="", help="screenshot only the region around the first match of this selector")
    ap.add_argument("--no-full", action="store_true", help="viewport-only screenshot")
    ap.add_argument("--out", default="shots")
    ap.add_argument("--tag", default="", help="suffix for the png names")
    return ap.parse_args()


def render(page, url, args, selectors):
    errors = []
    page.on("pageerror", lambda e: errors.append(str(e)))
    page.on("console", lambda m: errors.append(f"console.error: {m.text}") if m.type == "error" else None)
    resp = page.goto(BASE + url, wait_until="networkidle", timeout=90000)
    if args.click:
        page.click(args.click)
        page.wait_for_load_state("networkidle")
    page.wait_for_timeout(args.wait)
    return {
        "url": url, "status": resp.status if resp else None,
        "compat": page.evaluate("document.compatMode"),
        "doc_height": page.evaluate("document.documentElement.scrollHeight"),
        "body_height": page.evaluate("document.body.scrollHeight"),
        "measure": {sel: page.evaluate(MEASURE_JS, sel) for sel in selectors},
        "rects": page.evaluate(DUMP_JS) if args.dump else None,
        "errors": [e for e in errors if "favicon" not in e][:5],
    }


def screenshot(page, png, args):
    if args.unroll:
        page.add_style_tag(content=UNROLL_CSS)
        page.wait_for_timeout(300)
    clip = page.evaluate(CLIP_JS, args.clip) if args.clip else None
    if clip:
        png = png.with_name(png.stem + "_clip.png")
    page.screenshot(path=str(png), full_page=not args.no_full, clip=clip)
    return png


def main():
    args = parse_args()
    selectors = [s for s in args.measure.split(",") if s]
    patches = [tuple(p.split("=>", 1)) for p in args.patch]
    out_dir = Path(args.out)
    out_dir.mkdir(exist_ok=True)
    results = []
    with sync_playwright() as p:
        browser = launch(p)
        for mode in args.modes.split(","):
            ctx = new_context(browser, mode, patches, {"width": args.vw, "height": args.vh})
            page = ctx.new_page()
            for url in args.urls:
                info = {"mode": mode, **render(page, url, args, selectors)}
                png = screenshot(page, out_dir / f"{slug(url)}__{mode}_{args.vh}{args.tag}.png", args)
                info["png"] = str(png)
                results.append(info)
                print(json.dumps({k: v for k, v in info.items() if k != "rects"}), flush=True)
                (out_dir / "results.json").write_text(json.dumps(results, indent=1))
            ctx.close()
        browser.close()


if __name__ == "__main__":
    sys.exit(main())
