#!/usr/bin/env python3
"""
Generate a stripped, bgzipped + tabixed dbNSFP file for use as a VEP plugin
database.

Replaces dbnsfp_grch37_strip.sh and dbnsfp_grch38_strip.sh, and adds
T2T-CHM13v2.0 support (issue #1131).

The script reads the header of a dbNSFP variant file to look up column indices
by name (so it survives dbNSFP version changes that reorder/add columns), then
emits a shell pipeline that uses ``cut`` / ``awk`` / ``sort`` / ``bgzip`` /
``tabix`` to do the actual work (much faster than streaming through Python).

The ``awk`` stage collapses the transcript-parallel fields in
``COLLAPSE_COLUMNS`` to their distinct values — see that constant for why.

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

DEDUP_AWK = """BEGIN {
    FS = OFS = "\\t"
    n = split(PIPE_COLS, cols, ",")
    for (i = 1; i <= n; i++) pipe_col[cols[i] + 0] = 1
    n = split(ENTRY_COLS, cols, ",")
    for (i = 1; i <= n; i++) entry_col[cols[i] + 0] = 1
}

$1 ~ /^#/ { print; next }

{
    for (c in pipe_col) $c = collapse($c, 1)
    for (c in entry_col) $c = collapse($c, 0)
    print
}

function trim(s) {
    gsub(/^[ \\t]+/, "", s)
    gsub(/[ \\t]+$/, "", s)
    return s
}

# Distinct values of a ";"-separated array, first-seen order, dropping the "." dbNSFP uses for
# missing. split_pipe also splits each entry on "|" and rejoins on it.
function collapse(v, split_pipe,    n, m, i, j, a, b, tok, out, cnt, sep) {
    if (v == "." || v == "") return "."
    if (index(v, ";") == 0 && !(split_pipe && index(v, "|"))) return v

    delete seen
    out = ""
    cnt = 0
    sep = split_pipe ? "|" : ";"
    n = split(v, a, ";")
    for (i = 1; i <= n; i++) {
        if (split_pipe) {
            m = split(a[i], b, "|")
            for (j = 1; j <= m; j++) {
                tok = trim(b[j])
                if (tok != "" && tok != "." && !(tok in seen)) {
                    seen[tok] = 1
                    out = cnt++ ? out sep tok : tok
                }
            }
        } else {
            tok = trim(a[i])
            if (tok != "" && tok != "." && !(tok in seen)) {
                seen[tok] = 1
                out = cnt++ ? out sep tok : tok
            }
        }
    }
    return cnt ? out : "."
}
"""


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

    # Post-cut positions of the columns the dedup awk rewrites
    pipe_cols = []
    entry_cols = []
    for name, split_pipe in COLLAPSE_COLUMNS.items():
        col = post_cut_header.index(name) + 1
        (pipe_cols if split_pipe else entry_cols).append(str(col))

    out_file = f"dbNSFP{version}_{spec.short}.stripped"
    out_gz = f"{out_file}.gz"
    header_file = f"dbnsfp_header_{spec.short}.txt"  # per-build, parallel-safe
    dedup_awk_file = f"dbnsfp_dedup_{spec.short}.awk"

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
        # Read just the first line via process substitution so when we close the
        # pipe early, zcat's SIGPIPE doesn't trip `set -o pipefail` (it would if
        # we used `zcat | head -n1 | ...`).
        f"read -r _hdr < <(zcat {shlex.quote(files[0])})",
        "printf '%s\\n' \"$_hdr\" \\",
        f"    | cut -f {cut_arg} \\",
        f"    | awk 'BEGIN{{OFS=\"\\t\"}}{{ $1=\"#\"$1; print }}' > {header_file}",
        "",
        f"echo '[{build}] writing dedup awk -> {dedup_awk_file} "
        f"(collapsing {', '.join(COLLAPSE_COLUMNS)})'",
        f"cat > {dedup_awk_file} <<'DEDUP_AWK_EOF'",
        DEDUP_AWK.rstrip("\n"),
        "DEDUP_AWK_EOF",
        "",
    ]

    file_list = " ".join(shlex.quote(f) for f in files)

    pipeline = [f"zgrep -h -v '^#chr' {file_list}"]
    cut_cmd = f"cut -f {cut_arg}"

    # Runs before the sort so there is ~30% less data to sort
    dedup_cmd = (f"awk -v PIPE_COLS={','.join(pipe_cols)} -v ENTRY_COLS={','.join(entry_cols)}"
                 f" -f {dedup_awk_file}")

    if sorted_input:
        stage_desc = f"echo '[{build}] cutting columns + collapsing + bgzipping -> {out_gz}'"
        pipeline += [cut_cmd, dedup_cmd]
    else:
        src_chr_idx = column_index(header, spec.chr_col)
        awk_cmd = f"awk -F'\\t' '${src_chr_idx} != \".\"'"
        sort_cmd = f"sort -T \"$TMP_DIR\" -k{seq_col},{seq_col} -k{pos_col},{pos_col}n"
        stage_desc = (f"echo '[{build}] concat + collapse + sort + bgzip "
                      f"({len(files)} input files) -> {out_gz}'")
        pipeline += [awk_cmd, cut_cmd, dedup_cmd, sort_cmd]
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


def main():
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
