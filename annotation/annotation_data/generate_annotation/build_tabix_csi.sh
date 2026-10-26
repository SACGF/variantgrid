#!/bin/bash
# Build a fine-binned CSI index (.csi) next to the .tbi of every tabix-indexed annotation file (#2499).
#
# htslib opens the .csi in preference to the .tbi, and its small bins (4 kb by default, tabix -m 12) cut the
# per-lookup decode that dominates VEP's --custom/plugin time (the .tbi's linear index is fixed at 16 kb). The .tbi
# stays: VEP plugins check for it, and the column settings for TSV sources are read from it so the .csi indexes the
# same columns.
#
# Usage: build_tabix_csi.sh [--dry-run] [--jobs N] [--min-shift N] [dir ...]
#   dirs default to ${VEP_ANNOTATION_DIR}/annotation_data/{GRCh37,GRCh38,T2T-CHM13v2.0}, VEP_ANNOTATION_DIR being the
#   settings' ANNOTATION_VEP_BASE_DIR as for the sibling scripts here. Every <file>.tbi there without a <file>.csi is
#   indexed, JOBS at a time (default 2), at low CPU and I/O priority. Each build reads the whole file (about 70 MB/s
#   per job), so budget hours for the 300 GB of GRCh38 sources. Progress goes to build_tabix_csi.log
#   in each dir. Safe to re-run: existing .csi files are skipped and a failed build leaves no .csi behind.
#   --min-shift sets the bin size, 2^N bp. The default 12 (4 kb) keeps about 85% of the gain of 10 (1 kb) for a
#   quarter of the index, which htslib holds in memory in every VEP process (+0.34 GB vs +1.4 GB for 20 files).
set -u

JOBS=2; DRY_RUN=0; MIN_SHIFT=12; DIRS=()
while [ $# -gt 0 ]; do
  case $1 in
    --dry-run) DRY_RUN=1;;
    --jobs) JOBS=$2; shift;;
    --min-shift) MIN_SHIFT=$2; shift;;
    -h|--help) sed -n '2,16p' "$0"; exit 0;;
    *) DIRS+=("$1");;
  esac
  shift
done
if [ ${#DIRS[@]} -eq 0 ]; then
  if [ -z "${VEP_ANNOTATION_DIR:-}" ]; then
    echo "Set VEP_ANNOTATION_DIR (the settings' ANNOTATION_VEP_BASE_DIR) or pass the annotation dirs" >&2
    exit 1
  fi
  for build in GRCh37 GRCh38 T2T-CHM13v2.0; do
    [ -d "${VEP_ANNOTATION_DIR}/annotation_data/$build" ] && DIRS+=("${VEP_ANNOTATION_DIR}/annotation_data/$build")
  done
fi

# tabix arguments equivalent to the .tbi's own column settings (format 2 = VCF preset; the UCSC flag = 0-based)
tabix_args_from_tbi() {
  python3 - "$1" <<'EOF'
import gzip, struct, sys
with gzip.open(sys.argv[1], "rb") as f:
    assert f.read(4) == b"TBI\x01", "not a .tbi"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip = struct.unpack("<7i", f.read(28))
if fmt & 0xffff == 2:
    print("-p vcf")
else:
    args = ["-s", str(col_seq), "-b", str(col_beg), "-e", str(col_end or col_beg), "-S", str(skip), "-c", chr(meta)]
    if fmt & 0x10000:
        args.append("-0")
    print(" ".join(args))
EOF
}

build_one() {
  local tbi=$1 f=${1%.tbi} log
  log=$(dirname "$f")/build_tabix_csi.log
  local args; args=$(tabix_args_from_tbi "$tbi") || { echo "SKIP $f: cannot read $tbi" | tee -a "$log"; return 1; }
  if [ "$DRY_RUN" = 1 ]; then echo "tabix -C -m $MIN_SHIFT $args $f"; return 0; fi
  local start=$(date +%s)
  echo "[$(date '+%F %T')] start $f (-m $MIN_SHIFT $args)" >> "$log"
  if nice -n 10 ionice -c2 -n7 tabix -C -m "$MIN_SHIFT" $args "$f" 2>>"$log"; then
    echo "[$(date '+%F %T')] built $f.csi in $(( $(date +%s) - start ))s, $(stat -c %s "$f.csi") bytes, $(tabix -l "$f" | wc -l) contigs" | tee -a "$log"
  else
    rm -f "$f.csi"
    echo "[$(date '+%F %T')] FAILED $f (see $log)" | tee -a "$log"
    return 1
  fi
}
export -f tabix_args_from_tbi build_one
export DRY_RUN MIN_SHIFT

command -v tabix >/dev/null || { echo "tabix not found" >&2; exit 1; }
todo=()
for d in "${DIRS[@]}"; do
  for tbi in "$d"/*.tbi; do
    [ -e "$tbi" ] || continue
    f=${tbi%.tbi}
    [ -e "$f" ] || { echo "SKIP $tbi: no data file"; continue; }
    [ -e "$f.csi" ] && continue
    todo+=("$tbi")
  done
done
[ ${#todo[@]} -eq 0 ] && { echo "nothing to do: every tabix file already has a .csi"; exit 0; }
echo "${#todo[@]} file(s) to index, $JOBS at a time"
printf '%s\n' "${todo[@]}" | xargs -P "$JOBS" -I{} bash -c 'build_one "$1"' _ {}
