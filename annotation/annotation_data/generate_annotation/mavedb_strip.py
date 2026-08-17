#!/usr/bin/env python3
"""
Generate a stripped, bgzipped + tabixed MaveDB file for use as a VEP plugin
database.

The Ensembl MaveDB download carries 246 columns, and MaveDB.pm's parse_data
turns every record in a variant's span into a Perl hash keyed by all of them.
At a deep mutational scanning locus that is hundreds of thousands of records
for a single long allele, which is enough to OOM the VEP process.

Keeping the columns the plugin emits leaves VEP output identical while making
that hash five keys wide instead of 241.

Usage:
    mavedb_strip.py MaveDB_variants_2026-04-30.tsv.gz

    # keep extra data columns (the 5 key columns are always kept)
    mavedb_strip.py --columns urn,score,nt,pro,doi,hgvs,refseq \\
        MaveDB_variants_2026-04-30.tsv.gz
"""

import argparse
import gzip
import os
import shlex
import subprocess
import sys

# MaveDB.pm parse_data unpacks these positionally and splices them off the header before
# building colnames, so they must stay first and in this order
KEY_COLUMNS = ["#chr", "start", "end", "ref", "alt"]

# The columns MaveDB.pm emits, i.e. its default "cols" (urn:score:nt:pro before VEP 115). Keeping
# exactly these leaves the VEP output unchanged and means no cols= is needed on the plugin.
# VariantGrid itself reads only urn and score - see annotation/vep_columns.py.
DEFAULT_COLUMNS = ["urn", "score", "nt", "pro", "doi"]

# run() also reads hgvs and refseq, each behind a plugin flag VariantGrid turns off
# (single_aminoacid_changes=0, transcript_match=0). Add them via --columns for a file that stays
# correct with those flags on.
FLAG_COLUMNS = ["hgvs", "refseq"]

# Matches the index on the Ensembl download: 1-based start/end, "#" comment char
TABIX_ARGS = ["-s", "1", "-b", "2", "-e", "3"]


def read_header(path: str) -> list[str]:
    with gzip.open(path, "rt") as f:
        line = f.readline().rstrip("\n")
    if not line:
        raise SystemExit(f"Empty header in {path}")
    return line.split("\t")


def column_indices(header: list[str], columns: list[str]) -> list[int]:
    """ 1-based indices suitable for `cut -f`, key columns first then data columns in file order. """
    wanted = KEY_COLUMNS + columns
    if missing := [c for c in wanted if c not in header]:
        raise SystemExit(f"Column(s) not found in MaveDB header: {', '.join(missing)}")
    if header[:len(KEY_COLUMNS)] != KEY_COLUMNS:
        raise SystemExit(f"Expected MaveDB file to start with {KEY_COLUMNS}, got {header[:len(KEY_COLUMNS)]}")
    return sorted(header.index(c) + 1 for c in wanted)


def default_output(input_filename: str) -> str:
    basename = os.path.basename(input_filename)
    stem = basename[:-len(".tsv.gz")] if basename.endswith(".tsv.gz") else basename
    return os.path.join(os.path.dirname(input_filename), f"{stem}.stripped.tsv.gz")


def run(command: str):
    print(f"+ {command}", file=sys.stderr)
    subprocess.run(command, shell=True, check=True, executable="/bin/bash")


def count_data_rows(path: str) -> int:
    with gzip.open(path, "rt") as f:
        return sum(1 for line in f if not line.startswith("#"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="Ensembl MaveDB download, e.g. MaveDB_variants_2026-04-30.tsv.gz")
    ap.add_argument("--columns", default=",".join(DEFAULT_COLUMNS),
                    help=f"Comma-separated data columns to keep (default: {','.join(DEFAULT_COLUMNS)}, "
                         f"the set the plugin emits). The {len(KEY_COLUMNS)} key columns are always "
                         f"kept. Add {','.join(FLAG_COLUMNS)} to support single_aminoacid_changes=1 / "
                         f"transcript_match=1.")
    ap.add_argument("--output", default=None,
                    help="Output file (default: alongside the input, with '.stripped' before '.tsv.gz')")
    ap.add_argument("--verify", action="store_true",
                    help="Re-read both files afterwards and check the record counts match")
    args = ap.parse_args()

    if not os.path.isfile(args.input):
        raise SystemExit(f"Input file not found: {args.input}")

    columns = [c for c in args.columns.split(",") if c]
    header = read_header(args.input)
    indices = column_indices(header, columns)
    out_gz = args.output or default_output(args.input)

    print(f"[MaveDB] {args.input} -> {out_gz}")
    print(f"[MaveDB] keeping {len(indices)} of {len(header)} columns: "
          f"{','.join(header[i - 1] for i in indices)}")

    cut_arg = ",".join(str(i) for i in indices)
    run(f"zcat {shlex.quote(args.input)} | cut -f {cut_arg} | bgzip -c > {shlex.quote(out_gz)}")
    run(f"tabix -f {' '.join(TABIX_ARGS)} {shlex.quote(out_gz)}")

    if args.verify:
        print("[MaveDB] counting records (reads both files in full)")
        before = count_data_rows(args.input)
        after = count_data_rows(out_gz)
        if before != after:
            raise SystemExit(f"Record count changed: {before} in, {after} out")
        print(f"[MaveDB] {after} records, unchanged")

    run(f"ls -lh {shlex.quote(out_gz)} {shlex.quote(out_gz)}.tbi")
    print("[MaveDB] OK")
    if not set(DEFAULT_COLUMNS).issubset(columns):
        # The plugin validates cols= against the file header and defaults to columns no longer there
        print(f"[MaveDB] this set needs an explicit cols={':'.join(columns)} on the plugin, and emits "
              f"fewer fields than the unstripped file")


if __name__ == "__main__":
    main()
