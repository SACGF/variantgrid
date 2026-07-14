"""
PoC: annotate SV conservation (phastCons/phyloP max) with pyBigWig instead of VEP --custom bigwig.

VEP's --custom ...,summary_stats=max,format=bigwig computes the max bigWig value over the
variant's reference span, 1-based inclusive [POS, END]  ==  0-based half-open [POS-1, END).

Usage: poc_bigwig.py <query.vcf.gz> [--exact 0|1]  ->  writes TSV of maxes to stdout, timing to stderr
"""
import sys, time, gzip
import pyBigWig

D = "/data/annotation/VEP/annotation_data/GRCh37"
TRACKS = [
    ("phastCons100way", f"{D}/hg19.100way.phastCons.bw"),
    ("phastCons46way",  f"{D}/hg19.phastCons46way.placental.bw"),
    ("phyloP100way",    f"{D}/hg19.100way.phyloP100way.bw"),
    ("phyloP46way",     f"{D}/hg19.phyloP46way.placental.bw"),
]

# Hybrid: exact on the edge windows (where zoom bins can straddle the query boundary and
# overestimate a MAX), zoom on the fully-interior middle (interior zoom bins are safe for max).
# EDGE must exceed the largest zoom-reduction bin the file uses; 128 kb is comfortably safe here.
EDGE = 128_000

# Aligned-block decomposition (your idea): cover [s,e) with the fewest aligned windows -
# exact per-base on the two ragged ends (< BLOCK wide each), and the interior read straight
# from the zoom summaries at an ALIGNED BLOCK-wide granularity (nBins), which is bit-exact
# because each requested bin lands on precomputed bin boundaries (no straddling => no overshoot).
BLOCK = 16_384

def span_max(bw, chrom, s, e, mode):
    if mode == "exact":
        return bw.stats(chrom, s, e, type="max", exact=True)[0]
    if mode == "zoom":
        return bw.stats(chrom, s, e, type="max", exact=False)[0]
    if mode == "hybrid":
        if e - s <= 2 * EDGE:
            return bw.stats(chrom, s, e, type="max", exact=True)[0]
        vals = [
            bw.stats(chrom, s, s + EDGE, type="max", exact=True)[0],
            bw.stats(chrom, e - EDGE, e, type="max", exact=True)[0],
            bw.stats(chrom, s + EDGE, e - EDGE, type="max", exact=False)[0],
        ]
        vals = [v for v in vals if v is not None]
        return max(vals) if vals else None
    # aligned: exact per-base on ragged ends of width >= EDGE, and the interior read from zoom
    # summaries on a FINE grid (bin width G). A fine bin's zoom reduction is <= G, so its worst-case
    # straddle reaches at most ~G outside the interior - which lands inside the exact-scanned ends
    # (EDGE >> G) and is capped there. Net: bit-exact, but the whole span is read from zoom + two
    # small exact ends instead of a full per-base scan.
    G = 16384       # interior zoom bin width -> zoom reduction used is <= G
    E = 2 * G       # exact-end width must exceed that reduction so any straddle is capped exactly
    if e - s <= 2 * E:
        return bw.stats(chrom, s, e, type="max", exact=True)[0]
    lo, hi = s + E, e - E
    vals = [
        bw.stats(chrom, s, lo, type="max", exact=True)[0],
        bw.stats(chrom, hi, e, type="max", exact=True)[0],
    ]
    nb = max(1, (hi - lo) // G)
    vals.extend(v for v in bw.stats(chrom, lo, hi, type="max", exact=False, nBins=nb) if v is not None)
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None

def main():
    query = sys.argv[1]
    mode = "exact"
    if "--mode" in sys.argv:
        mode = sys.argv[sys.argv.index("--mode") + 1]
    exact = mode == "exact"

    bws = [(name, pyBigWig.open(path)) for name, path in TRACKS]
    # chrom naming: bigwig uses "chr1"? decide per file
    prefixes = [("chr" if "chr1" in bw.chroms() else "") for _, bw in bws]
    chromlens = [bw.chroms() for _, bw in bws]

    # read query
    variants = []
    with gzip.open(query, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            chrom, pos, info = c[0], int(c[1]), c[7]
            kv = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
            end = int(kv.get("END", pos))
            svlen = abs(int(kv["SVLEN"])) if "SVLEN" in kv else 0
            # VEP's reference footprint: 1-based inclusive [POS, max(END, POS+|SVLEN|)]
            # (DEL/DUP -> END; INS with tiny END but large SVLEN -> POS+|SVLEN|)
            eff_end = max(end, pos + svlen)
            variants.append((chrom, pos, eff_end, kv.get("SVTYPE", ".")))

    t0 = time.time()
    out = []
    for chrom, pos, end in ((v[0], v[1], v[2]) for v in variants):
        vals = []
        for (name, bw), pre, clen in zip(bws, prefixes, chromlens):
            ck = pre + chrom
            s = pos - 1                      # 0-based inclusive start
            e = end                          # 0-based exclusive end == 1-based inclusive END
            L = clen.get(ck)
            if L is None:
                vals.append(None); continue
            if s < 0: s = 0
            if e > L: e = L
            if e <= s: e = s + 1
            try:
                v = span_max(bw, ck, s, e, mode)
            except Exception:
                v = None
            vals.append(v)
        out.append(vals)
    dt = time.time() - t0

    print(f"PoC pyBigWig mode={mode}: {len(variants)} SVs x {len(TRACKS)} tracks "
          f"in {dt*1000:.1f} ms ({dt/max(1,len(variants))*1000:.4f} ms/SV)", file=sys.stderr)

    sys.stdout.write("chrom\tpos\tend\tsvtype\t" + "\t".join(n for n, _ in TRACKS) + "\n")
    for v, vals in zip(variants, out):
        cells = ["" if x is None else repr(x) for x in vals]
        sys.stdout.write("\t".join([v[0], str(v[1]), str(v[2]), v[3], *cells]) + "\n")

if __name__ == "__main__":
    main()
