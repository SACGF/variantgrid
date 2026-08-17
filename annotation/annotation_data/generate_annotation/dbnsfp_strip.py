#!/usr/bin/env python3
"""
Generate a stripped, bgzipped + tabixed dbNSFP file for use as a VEP plugin
database.

Replaces dbnsfp_grch37_strip.sh and dbnsfp_grch38_strip.sh, and adds
T2T-CHM13v2.0 support (issue #1131).

The script reads the header of a dbNSFP variant file to look up column indices
by name (so it survives dbNSFP version changes that reorder/add columns), then
emits a shell pipeline that uses ``cut`` / ``sort`` / ``bgzip`` / ``tabix`` for
the bulk work, with one row-by-row stage handled by re-invoking this script as
``--filter``. That stage drops rows unmapped in the target build and collapses
the transcript-parallel fields in ``COLLAPSE_COLUMNS`` to their distinct values
— see that constant for why.

Modes:
* ``--sorted``: pass a single, already-sorted concatenated file (e.g. the
  pre-sorted ``dbNSFP5.3.1a_grch37.gz`` / ``dbNSFP5.3.1a_grch38.gz``
  downloads). Just cuts columns + bgzip + tabix. Fast.
* default (used for T2T): pass per-chrom dbNSFP files. Concats, drops rows
  with no mapping in this build, sorts, bgzips, tabixes.

Usage:
    # GRCh37 / GRCh38 — pre-sorted single download
    dbnsfp_strip.py --build GRCh38 --version 5.3.1a --sorted \\
        dbNSFP5.3.1a_grch38.gz

    # T2T — sort per-chrom inputs (shell expands the glob)
    dbnsfp_strip.py --build T2T-CHM13v2.0 --version 5.3.1a \\
        dbNSFP5.3.1a_variant.chr*.gz

Prints the script to stdout — redirect to a file and run with bash.
"""

import argparse
import gzip
import os
import shlex
import sys
from dataclasses import dataclass

GRCH37 = "GRCh37"
GRCH38 = "GRCh38"
T2T = "T2T-CHM13v2.0"

GENOME_BUILDS = (GRCH37, GRCH38, T2T)


@dataclass(frozen=True)
class BuildSpec:
    short: str            # short tag used in output filenames
    chr_col: str          # dbNSFP column name for chromosome
    pos_col: str          # dbNSFP column name for 1-based position


BUILDS: dict[str, BuildSpec] = {
    GRCH37: BuildSpec("grch37", "hg19_chr", "hg19_pos(1-based)"),
    GRCH38: BuildSpec("grch38", "#chr",     "pos(1-based)"),
    T2T:    BuildSpec("t2t",    "hs1_chr",  "hs1_pos(1-based)"),
}


# Key columns we always need (lookup keys for VEP)
KEY_COLUMNS = ["ref", "alt", "aaref", "aaalt", "Ensembl_transcriptid"]

# Existing rankscore / Aloft / conservation / domain columns we already pull
DBNSFP_COLUMNS = [
    "Aloft_Confidence",
    "Aloft_pred",
    "Aloft_prob_Dominant",
    "Aloft_prob_Recessive",
    "Aloft_prob_Tolerant",
    "AlphaMissense_pred",
    "AlphaMissense_rankscore",
    "AlphaMissense_score",
    "BayesDel_noAF_rankscore",
    "BayesDel_noAF_score",
    "CADD_phred",
    "CADD_raw",
    "CADD_raw_rankscore",
    "ClinPred_pred",
    "ClinPred_rankscore",
    "ClinPred_score",
    "GERP++_RS",
    "Interpro_domain",
    "MPC_score",
    "MetaLR_rankscore",
    "MetaRNN_pred",
    "MetaRNN_score",
    "MutPred2_score",
    "MutPred2_top5_mechanisms",
    "PrimateAI_pred",
    "PrimateAI_score",
    "REVEL_rankscore",
    "REVEL_score",
    "VARITY_ER_score",
    "VARITY_R_score",
    "VEST4_rankscore",
    "VEST4_score",
]


# dbNSFP stores these as ";"-separated arrays parallel to Ensembl_transcriptid. We run VEP with
# --refseq, so the dbNSFP plugin can never match a NM_ stable_id against Ensembl_transcriptid and
# transcript_match is unusable (it would delete every transcript-specific field) - leaving the plugin
# to emit the whole array into every consequence. Collapsing to distinct values here costs us nothing:
# BulkVEPVCFAnnotationInserter already reduces these to a set via remove_empty_multiples/join_uniq.
#
# The value is True where entries additionally repeat internally on "|" - Interpro_domain lists a
# domain once per Interpro predicting database, so GNAS goes from 8464 chars to 61.
COLLAPSE_COLUMNS = {
    "Interpro_domain": True,
    "MutPred2_top5_mechanisms": False,  # "|" separates the ranked top 5, so keep entries intact
}

MISSING = "."


def collapse(value: str, split_pipe: bool) -> str:
    """ Distinct values of a ";"-separated array, first-seen order, dropping the "." dbNSFP uses
        for missing. split_pipe also splits each entry on "|" and rejoins on it. """
    if value in (MISSING, ""):
        return MISSING
    if ";" not in value and not (split_pipe and "|" in value):
        return value

    entries = value.split(";")
    if split_pipe:
        entries = [token for entry in entries for token in entry.split("|")]
    distinct = dict.fromkeys(t for t in (e.strip() for e in entries) if t and t != MISSING)
    if not distinct:
        return MISSING
    return ("|" if split_pipe else ";").join(distinct)


def run_filter(pipe_cols: list[int], entry_cols: list[int], drop_missing_col: int | None):
    """ Row-by-row stage of the generated pipeline: reads stdin, writes stdout. Column numbers are
        1-based positions in the already-cut stream. """
    collapse_cols = [(c - 1, True) for c in pipe_cols] + [(c - 1, False) for c in entry_cols]
    drop_idx = drop_missing_col - 1 if drop_missing_col else None

    write = sys.stdout.write
    for line in sys.stdin:
        if line.startswith("#"):
            write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if drop_idx is not None and fields[drop_idx] == MISSING:
            continue
        changed = False
        for i, split_pipe in collapse_cols:
            collapsed = collapse(fields[i], split_pipe)
            if collapsed != fields[i]:
                fields[i] = collapsed
                changed = True
        write("\t".join(fields) + "\n" if changed else line)


def read_header(path: str) -> list[str]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        line = f.readline().rstrip("\n")
    if not line:
        raise ValueError(f"Empty header in {path}")
    return line.split("\t")


def column_index(header: list[str], name: str) -> int:
    """ 1-based index suitable for `cut -f` / `sort -k` / `tabix -s`. """
    try:
        return header.index(name) + 1
    except ValueError:
        raise SystemExit(f"Column {name!r} not found in dbNSFP header")


def build_pipeline(build: str, version: str, files: list[str], tmp_dir: str,
                   sorted_input: bool) -> str:
    spec = BUILDS[build]

    files = sorted(files)
    if sorted_input and len(files) != 1:
        raise SystemExit(
            f"--sorted expects a single pre-sorted file but got {len(files)}: {files}"
        )
    missing = [f for f in files if not os.path.isfile(f)]
    if missing:
        raise SystemExit(f"Input file(s) not found: {missing}")

    header = read_header(files[0])

    wanted = KEY_COLUMNS + DBNSFP_COLUMNS
    # Always include the build's chr+pos so the stripped file is tabix-able
    wanted = [spec.chr_col, spec.pos_col] + [c for c in wanted if c not in (spec.chr_col, spec.pos_col)]

    cut_indices = sorted({column_index(header, c) for c in wanted})
    cut_arg = ",".join(str(i) for i in cut_indices)

    # Re-derive post-cut positions of chr / pos for sort + tabix
    post_cut_header = [header[i - 1] for i in cut_indices]
    seq_col = post_cut_header.index(spec.chr_col) + 1
    pos_col = post_cut_header.index(spec.pos_col) + 1

    # Post-cut positions of the columns the filter stage collapses
    pipe_cols = []
    entry_cols = []
    for name, split_pipe in COLLAPSE_COLUMNS.items():
        col = post_cut_header.index(name) + 1
        (pipe_cols if split_pipe else entry_cols).append(str(col))

    out_file = f"dbNSFP{version}_{spec.short}.stripped"
    out_gz = f"{out_file}.gz"
    header_file = f"dbnsfp_header_{spec.short}.txt"  # per-build, parallel-safe

    # Leading "#" makes tabix treat it as a comment; the dbNSFP plugin strips it back off
    header_line = "#" + "\\t".join(post_cut_header)

    lines: list[str] = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        "# Loud failure: print exactly which command failed and where",
        "trap 'rc=$?; echo \"ERROR rc=$rc at line $LINENO running: $BASH_COMMAND\" >&2; "
        "exit $rc' ERR",
        "",
        f"echo '[{build} dbNSFP {version}] -> {out_gz}'",
        f"echo '[{build}] selected {len(cut_indices)} of {len(header)} columns"
        + ("  (input pre-sorted, skipping sort)" if sorted_input else "")
        + "'",
        "",
        f"echo '[{build}] writing stripped header -> {header_file}'",
        f"printf '%b\\n' {shlex.quote(header_line)} > {header_file}",
        "",
    ]

    file_list = " ".join(shlex.quote(f) for f in files)

    pipeline = [f"zgrep -h -v '^#chr' {file_list}"]
    cut_cmd = f"cut -f {cut_arg}"

    # Runs before the sort so there is ~30% less data to sort
    this_script = shlex.quote(os.path.abspath(__file__))
    filter_cmd = (f"{shlex.quote(sys.executable)} {this_script} --filter"
                  f" --pipe-cols {','.join(pipe_cols)} --entry-cols {','.join(entry_cols)}")

    if sorted_input:
        stage_desc = f"echo '[{build}] cutting columns + collapsing + bgzipping -> {out_gz}'"
        pipeline += [cut_cmd, filter_cmd]
    else:
        # Rows with no mapping in this build carry "." in the chr column - drop them in the same pass
        filter_cmd += f" --drop-missing-col {seq_col}"
        sort_cmd = f"sort -T \"$TMP_DIR\" -k{seq_col},{seq_col} -k{pos_col},{pos_col}n"
        stage_desc = (f"echo '[{build}] concat + collapse + sort + bgzip "
                      f"({len(files)} input files) -> {out_gz}'")
        pipeline += [cut_cmd, filter_cmd, sort_cmd]
        lines += [
            f"TMP_DIR={shlex.quote(tmp_dir)}",
            'mkdir -p "$TMP_DIR"',
            "",
        ]

    pipeline += [f"cat {header_file} -",
                 f"bgzip -c > {out_gz}"]

    lines += [
        stage_desc,
        " \\\n    | ".join(pipeline),
        "",
        f"echo '[{build}] building tabix index'",
        f"tabix -f -s {seq_col} -b {pos_col} -e {pos_col} {out_gz}",
        "",
        f"echo '[{build}] OK'; ls -lh {out_gz} {out_gz}.tbi",
    ]
    return "\n".join(lines) + "\n"


def filter_main(argv: list[str]):
    ap = argparse.ArgumentParser(
        prog="dbnsfp_strip.py --filter",
        description="Row-by-row stage of the generated pipeline (stdin -> stdout). Column numbers "
                    "are 1-based positions in the already-cut stream."
    )
    ap.add_argument("--pipe-cols", default="", type=_column_list,
                    help="Columns to collapse, also splitting/rejoining entries on '|'")
    ap.add_argument("--entry-cols", default="", type=_column_list,
                    help="Columns to collapse, keeping entries intact")
    ap.add_argument("--drop-missing-col", type=int, default=None,
                    help="Drop rows where this column is '.' (unmapped in this build)")
    args = ap.parse_args(argv)
    run_filter(args.pipe_cols, args.entry_cols, args.drop_missing_col)


def _column_list(value: str) -> list[int]:
    return [int(c) for c in value.split(",") if c]


def main():
    if sys.argv[1:2] == ["--filter"]:
        filter_main(sys.argv[2:])
        return

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", required=True, choices=GENOME_BUILDS,
                    help="Genome build to extract")
    ap.add_argument("--version", required=True,
                    help="dbNSFP version string used in output filename, e.g. 5.3.1a")
    ap.add_argument("inputs", nargs="+",
                    help="dbNSFP input file(s). With --sorted, a single "
                         "pre-sorted file. Otherwise the per-chrom files "
                         "(e.g. dbNSFP5.3.1a_variant.chr*.gz — let the shell "
                         "expand the glob).")
    ap.add_argument("--sorted", dest="sorted_input", action="store_true",
                    help="Input is a single already-sorted file (e.g. the "
                         "pre-sorted GRCh37/GRCh38 dbNSFP downloads). "
                         "Skips the sort step.")
    ap.add_argument("--tmp-dir", default=None,
                    help="Sort scratch dir, only used without --sorted "
                         "(default: /tmp/dbnsfp_<build>)")
    args = ap.parse_args()

    tmp_dir = args.tmp_dir or f"/tmp/dbnsfp_{BUILDS[args.build].short}"
    script = build_pipeline(args.build, args.version, args.inputs,
                            tmp_dir, args.sorted_input)
    sys.stdout.write(script)


if __name__ == "__main__":
    main()
