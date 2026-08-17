"""Generate sparse test reference fastas from recorded fetch regions.

Usage:
    python3 scripts/generate_sparse_test_fastas.py <regions.jsonl> <output_dir>

See variantgrid/data/reference/sparse_test_fastas/README.md - the regions file is
written by running the test suite under variantgrid.test_runner.FastaRecordingRunner
on a machine with the real reference fastas.

For each source fasta this writes a bgzipped sparse copy: identical contig
names/order/lengths, all-N sequence except the recorded regions (plus padding)
which carry the real sequence. Contigs under 100kb that were fetched at all
(eg mitochondria) are included in full.

Every recorded region is verified against the real fasta after generation.
"""
import json
import os
import sys
from collections import defaultdict

import pysam

REGIONS_FILE = sys.argv[1]
OUT_DIR = sys.argv[2]
PAD = 10_000
FULL_CONTIG_MAX = 100_000  # contigs this small with any fetch recorded are included whole
LINE = 60

regions_by_file = defaultdict(lambda: defaultdict(list))
with open(REGIONS_FILE) as f:
    for line in f:
        filename, contig, start, end = json.loads(line)
        regions_by_file[filename][contig].append((start, end))

os.makedirs(OUT_DIR, exist_ok=True)

for filename, contig_regions in sorted(regions_by_file.items()):
    real = pysam.FastaFile(filename)
    lengths = dict(zip(real.references, real.lengths))
    out_path = os.path.join(OUT_DIR, os.path.basename(filename))
    print(f"== {filename} -> {out_path}")

    # merge padded intervals per contig
    merged = {}
    for contig, regs in contig_regions.items():
        length = lengths[contig]
        if length <= FULL_CONTIG_MAX:
            merged[contig] = [(0, length)]
            continue
        ivs = []
        for start, end in regs:
            start = 0 if start is None else start
            end = length if end is None else end
            ivs.append((max(0, start - PAD), min(length, end + PAD)))
        ivs.sort()
        out = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[contig] = [tuple(x) for x in out]

    n_line = ("N" * LINE + "\n").encode()
    with pysam.BGZFile(out_path, "wb") as out:
        for contig in real.references:
            length = lengths[contig]
            out.write(f">{contig}\n".encode())
            ivs = merged.get(contig)
            if not ivs:
                full_lines, rem = divmod(length, LINE)
                block = n_line * 4096
                while full_lines >= 4096:
                    out.write(block)
                    full_lines -= 4096
                out.write(n_line * full_lines)
                if rem:
                    out.write(b"N" * rem + b"\n")
                continue
            seq = bytearray(b"N" * length)
            for s, e in ivs:
                seq[s:e] = real.fetch(contig, s, e).upper().encode()
            CHUNK_LINES = 16384
            buf = []
            for i in range(0, length, LINE):
                buf.append(seq[i:i + LINE])
                if len(buf) == CHUNK_LINES:
                    out.write(b"\n".join(buf) + b"\n")
                    buf = []
            if buf:
                out.write(b"\n".join(buf) + b"\n")
            del seq

    pysam.faidx(out_path)

    # verify every recorded region round-trips
    sparse = pysam.FastaFile(out_path)
    assert list(sparse.references) == list(real.references)
    assert list(sparse.lengths) == list(real.lengths)
    checked = 0
    for contig, regs in contig_regions.items():
        for start, end in regs:
            s = 0 if start is None else start
            e = lengths[contig] if end is None else end
            want = real.fetch(contig, s, e).upper()
            got = sparse.fetch(contig, s, e).upper()
            assert got == want, f"mismatch {contig}:{s}-{e}"
            checked += 1
    print(f"   contigs={len(real.references)} patched_contigs={len(merged)} regions_verified={checked}")
    for ext in (".fai", ".gzi"):
        sz = os.path.getsize(out_path + ext)
        print(f"   {os.path.basename(out_path)}{ext}: {sz:,} bytes")
    print(f"   {os.path.basename(out_path)}: {os.path.getsize(out_path):,} bytes")
