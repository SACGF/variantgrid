"""
Diff the element boxes shot.py --dump recorded in quirks vs standards mode for each URL: elements that
collapsed (visible to zero height), appeared, resized by more than --min-delta, or changed font size.
Resized and font lines are innermost elements only, so one taller cell isn't reported as its whole table.

    $VG_BROWSER_HOME/pw/bin/python scripts/browser/compare.py shots/results.json [--min-delta 8] [--top 12]
"""
import argparse
import json
from collections import defaultdict


def short(key: str) -> str:
    return "/".join(key.split("/")[-4:])[-150:]


def leaves(items, key_of):
    keys = {key_of(i) for i in items}
    return [i for i in items if not any(k != key_of(i) and k.startswith(key_of(i) + "/") for k in keys)]


def report(url, quirks, standards, min_delta, top):
    print(f"\n=== {url}  status q={quirks['status']} s={standards['status']}")
    qr, sr = quirks.get("rects"), standards.get("rects")
    if not (qr and sr):
        print("   (run shot.py with --dump and both modes)")
        return
    print(f"   page height quirks {quirks['body_height']} standards {standards['body_height']}; "
          f"errors q={quirks['errors']} s={standards['errors']}")
    collapsed, appeared, resized, fonts = [], [], [], []
    for key in qr.keys() & sr.keys():
        qw, qh, _, _, qf = qr[key]
        sw, sh, _, _, sf = sr[key]
        if qh > 4 and sh <= 1:
            collapsed.append((key, qh))
        elif sh > 4 and qh <= 1:
            appeared.append((key, sh))
        elif abs(qh - sh) >= min_delta or abs(qw - sw) >= min_delta:
            resized.append((key, (qw, qh), (sw, sh)))
        if qf != sf:
            fonts.append((key, qf, sf))
    print(f"   elements {len(qr)}: collapsed {len(collapsed)}, appeared {len(appeared)}, "
          f"resized>={min_delta}px {len(resized)}, font-size changed {len(fonts)}")
    for name, items in (("COLLAPSED", collapsed), ("APPEARED", appeared)):
        for key, h in items[:top]:
            print(f"   {name} {h}px  {short(key)}")
    for key, q, s in sorted(leaves(resized, lambda i: i[0]), key=lambda i: -abs(i[1][1] - i[2][1]))[:top]:
        print(f"   RESIZED q={q} s={s}  {short(key)}")
    for key, qf, sf in leaves(fonts, lambda i: i[0])[:top]:
        print(f"   FONT q={qf} s={sf}  {short(key)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results")
    ap.add_argument("--min-delta", type=int, default=8)
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args()
    by_url = defaultdict(dict)
    for r in json.load(open(args.results)):
        by_url[r["url"]][r["mode"]] = r
    for url, modes in by_url.items():
        if "quirks" in modes and "standards" in modes:
            report(url, modes["quirks"], modes["standards"], args.min_delta, args.top)
        else:
            print(f"\n=== {url}: only {list(modes)} rendered")


if __name__ == "__main__":
    main()
