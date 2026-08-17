#!/usr/bin/env python3
"""De-identify a VCF by shifting variant positions and renaming samples.

Positions move by a random offset, and REF is re-read from the reference FASTA at the new
position so the record stays internally consistent. ALT is rebuilt to preserve the variant
type (SNV stays an SNV, a 3bp deletion stays a 3bp deletion).

Records that describe the same underlying event move together: rows sharing an MNVTAG, a
phase set (PS), or a position get one shift between them, so decomposed MNVs and phased
calls keep their relationships. Everything else gets an independent offset, which is what
destroys the inter-variant spacing that makes a genotype set matchable.

Shifts are bounded (--max-shift) so variants stay in roughly the same part of the same gene.
Pass --bed to constrain them to the feature they started in.

Usage:
    scripts/vcf_fuzz.py in.vcf --fasta ref.fa.gz --rename-sample OLD=NEW -o out.vcf
"""
import argparse
import random
import re
import sys

import pysam

BASES = "ACGT"
SYMBOLIC_ALT = re.compile(r"^[<.*]")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf")
    parser.add_argument("--fasta", required=True, help="Reference FASTA (bgzipped + faidx)")
    parser.add_argument("-o", "--out", help="Output VCF (default: overwrite input)")
    parser.add_argument("--max-shift", type=int, default=200,
                        help="Maximum bp to move a variant in either direction (default 200). "
                             "0 renames samples and leaves coordinates alone - right for files "
                             "whose positions are fixed panel intervals rather than patient data")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed, for reproducible output")
    parser.add_argument("--rename-sample", action="append", default=[], metavar="OLD=NEW",
                        help="Rename a sample everywhere it appears, including header text")
    parser.add_argument("--replace", action="append", default=[], metavar="OLD=NEW",
                        help="Replace arbitrary header text - run dates, timezones, site paths "
                             "left in caller command lines. Repeatable")
    parser.add_argument("--bed", help="Keep each variant inside the feature it started in")
    return parser.parse_args()


def load_bed(path):
    features = {}
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            features.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    return features


def containing_feature(features, contig, pos):
    for start, end in features.get(contig, []):
        if start < pos <= end:
            return start, end
    return None


class Reference:
    """FASTA lookup that copes with the VCF and the FASTA disagreeing about 'chr' prefixes."""

    def __init__(self, path):
        self.fasta = pysam.FastaFile(path)
        self.names = set(self.fasta.references)

    def contig(self, contig):
        if contig in self.names:
            return contig
        alternative = contig[3:] if contig.startswith("chr") else f"chr{contig}"
        return alternative if alternative in self.names else None

    def fetch(self, contig, pos, length):
        """1-based pos, returns `length` bases or None if unavailable/ambiguous."""
        name = self.contig(contig)
        if name is None:
            return None
        try:
            seq = self.fasta.fetch(name, pos - 1, pos - 1 + length).upper()
        except (ValueError, KeyError):
            return None
        if len(seq) != length or "N" in seq:
            return None
        return seq


def rebuild_alt(old_ref, old_alt, new_ref, rng):
    """Keep the variant type and length change, against the new reference base(s)."""
    if SYMBOLIC_ALT.match(old_alt):
        return old_alt
    if len(old_alt) < len(old_ref):            # deletion
        return new_ref[:len(old_alt)]
    if len(old_alt) > len(old_ref):            # insertion
        return new_ref + old_alt[len(old_ref):]
    if old_alt != new_ref:                     # substitution, still a substitution
        return old_alt
    # The new reference happens to equal the old alt - pick something else, same length
    replacement = rng.choice([b for b in BASES if b != new_ref[0]])
    return replacement + old_alt[1:]


def group_key(fields):
    """Records describing one event share a key, so they get one shift between them."""
    contig, pos, _, _, _, _, _, info, *rest = fields
    if match := re.search(r"MNVTAG=([^;\t]+)", info):
        return ("mnv", match.group(1))
    if len(rest) >= 2 and "PS" in rest[0].split(":"):
        values = dict(zip(rest[0].split(":"), rest[1].split(":")))
        if phase_set := values.get("PS"):
            return ("ps", contig, phase_set)
    return ("pos", contig, pos)


def shifted_info(info, delta, anchor=None):
    """Shift positional INFO fields. `anchor` is the group's (ref, alt) after shifting,
        which MNVTAG quotes alongside the position. """
    def shift_end(match):
        return f"END={int(match.group(1)) + delta}"

    info = re.sub(r"END=(\d+)", shift_end, info)

    def shift_mnvtag(match):
        contig, pos, allele = match.group(1), int(match.group(2)) + delta, match.group(3)
        if anchor:
            allele = f"{anchor[0]}->{anchor[1]}"
        return f"MNVTAG={contig}:{pos}_{allele}"

    return re.sub(r"MNVTAG=([^:;\t]+):(\d+)_([^;\t]+)", shift_mnvtag, info)


def shifted_sample_columns(columns, delta):
    fmt = columns[0].split(":")
    if "PS" not in fmt:
        return columns
    out = [columns[0]]
    index = fmt.index("PS")
    for sample in columns[1:]:
        values = sample.split(":")
        if len(values) > index and values[index] not in (".", ""):
            values[index] = str(int(values[index]) + delta)
        out.append(":".join(values))
    return out


def main():
    args = parse_args()
    rng = random.Random(args.seed)
    reference = Reference(args.fasta)
    features = load_bed(args.bed) if args.bed else None
    renames = dict(pair.split("=", 1) for pair in args.rename_sample + args.replace)

    header, records = [], []
    with open(args.vcf) as f:
        for line in f:
            line = line.rstrip("\n")
            (header if line.startswith("#") else records).append(line)

    for old, new in renames.items():
        header = [line.replace(old, new) for line in header]

    if args.max_shift == 0:
        with open(args.out or args.vcf, "w") as f:
            for line in header + records:
                f.write(line + "\n")
        print(f"{args.vcf}: samples renamed, {len(records)} records left in place")
        return

    groups = {}
    for record in records:
        groups.setdefault(group_key(record.split("\t")), []).append(record)

    contig_order = {line.split("ID=")[1].split(",")[0]: i
                    for i, line in enumerate(header) if line.startswith("##contig=")}

    out_records, skipped = [], 0
    for members in groups.values():
        parsed = [r.split("\t") for r in members]
        longest = max(len(f[3]) for f in parsed)
        anchor = min(int(f[1]) for f in parsed)
        contig = parsed[0][0]

        delta = None
        feature = containing_feature(features, contig, anchor) if features else None
        for _ in range(200):
            candidate = rng.randint(-args.max_shift, args.max_shift)
            if candidate == 0:
                continue
            if any(int(f[1]) + candidate < 1 for f in parsed):
                continue
            if feature and not all(feature[0] < int(f[1]) + candidate <= feature[1] for f in parsed):
                continue
            spread = max(int(f[1]) for f in parsed) - anchor
            if reference.fetch(contig, anchor + candidate, spread + longest) is None:
                continue
            delta = candidate
            break

        if delta is None:
            skipped += len(parsed)
            continue

        for fields in parsed:
            new_pos = int(fields[1]) + delta
            new_ref = reference.fetch(contig, new_pos, len(fields[3]))
            fields[1] = str(new_pos)
            fields[4] = rebuild_alt(fields[3], fields[4], new_ref, rng)
            fields[3] = new_ref

        # MNVTAG names the group's longest allele at its anchor - rewrite it from the
        # record that actually carries it, now that every REF has been re-read.
        anchor_fields = max((f for f in parsed if int(f[1]) == anchor + delta),
                            key=lambda f: len(f[3]), default=None)
        anchor_alleles = (anchor_fields[3], anchor_fields[4]) if anchor_fields else None

        for fields in parsed:
            fields[7] = shifted_info(fields[7], delta, anchor_alleles)
            if len(fields) > 8:
                fields[8:] = shifted_sample_columns(fields[8:], delta)
            out_records.append("\t".join(fields))

    out_records.sort(key=lambda r: (contig_order.get(r.split("\t")[0], 1 << 30),
                                    int(r.split("\t")[1])))

    with open(args.out or args.vcf, "w") as f:
        for line in header:
            f.write(line + "\n")
        for record in out_records:
            f.write(record + "\n")

    print(f"{args.vcf}: {len(out_records)} records shifted"
          + (f", {skipped} skipped (no usable reference)" if skipped else ""))


if __name__ == "__main__":
    sys.exit(main())
