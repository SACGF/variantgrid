#!/bin/bash
set -euo pipefail

#
# Strip gnomAD v4.1 joint VCFs down to the INFO fields we need,
# rename _joint -> legacy names, and submit as SLURM jobs (one per chrom).
#
# Usage:
#   bash strip_gnomad41.sh /path/to/input/dir /path/to/output/dir PARTITION
#   bash strip_gnomad41.sh /path/to/input/dir /path/to/output/dir PARTITION --dry-run
#

INPUT_DIR="${1:?Usage: $0 INPUT_DIR OUTPUT_DIR PARTITION [--dry-run]}"
OUTPUT_DIR="${2:?Usage: $0 INPUT_DIR OUTPUT_DIR PARTITION [--dry-run]}"
PARTITION="${3:?Usage: $0 INPUT_DIR OUTPUT_DIR PARTITION [--dry-run]}"
DRY_RUN="${4:-}"

# HPC tool paths (directories containing bcftools, tabix, bgzip etc.)
TOOL_PATHS=(
    /hpcfs/groups/phoenix-hpc-sacgf/tools/bcftools/current
    /hpcfs/groups/phoenix-hpc-sacgf/tools/tabix-0.2.6
)

CHROMOSOMES="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"

# INFO fields to keep (v4.1 _joint names)
KEEP_FIELDS=(
    # Core counts
    AC_joint AN_joint nhomalt_joint
    # Filtering allele frequencies
    faf95_joint faf99_joint fafmax_faf95_max_joint fafmax_faf99_max_joint
    # grpmax
    AF_grpmax_joint AC_grpmax_joint AN_grpmax_joint grpmax_joint nhomalt_grpmax_joint
    # Per-population AC/AN
    AC_joint_afr AN_joint_afr
    AC_joint_amr AN_joint_amr
    AC_joint_asj AN_joint_asj
    AC_joint_eas AN_joint_eas
    AC_joint_fin AN_joint_fin
    AC_joint_nfe AN_joint_nfe
    AC_joint_sas AN_joint_sas
    AC_joint_remaining AN_joint_remaining
    AC_joint_mid AN_joint_mid
    # ChrX sex-specific
    AC_joint_XY AN_joint_XY AF_joint_XY
)

# Build the bcftools keep string: ^INFO/AC_joint,INFO/AN_joint,...
KEEP_STR=$(printf "INFO/%s," "${KEEP_FIELDS[@]}")
KEEP_STR="^${KEEP_STR%,}"  # strip trailing comma

mkdir -p "$OUTPUT_DIR"

# Write rename file: _joint -> legacy names for downstream compatibility
RENAME_FILE="$OUTPUT_DIR/rename_annots.txt"
cat > "$RENAME_FILE" <<'RENAME'
INFO/AC_joint	AC
INFO/AN_joint	AN
INFO/nhomalt_joint	nhomalt
INFO/faf95_joint	faf95
INFO/faf99_joint	faf99
INFO/fafmax_faf95_max_joint	fafmax_faf95_max
INFO/fafmax_faf99_max_joint	fafmax_faf99_max
INFO/AF_grpmax_joint	AF_grpmax
INFO/AC_grpmax_joint	AC_grpmax
INFO/AN_grpmax_joint	AN_grpmax
INFO/grpmax_joint	grpmax
INFO/nhomalt_grpmax_joint	nhomalt_grpmax
INFO/AC_joint_afr	AC_afr
INFO/AN_joint_afr	AN_afr
INFO/AC_joint_amr	AC_amr
INFO/AN_joint_amr	AN_amr
INFO/AC_joint_asj	AC_asj
INFO/AN_joint_asj	AN_asj
INFO/AC_joint_eas	AC_eas
INFO/AN_joint_eas	AN_eas
INFO/AC_joint_fin	AC_fin
INFO/AN_joint_fin	AN_fin
INFO/AC_joint_nfe	AC_nfe
INFO/AN_joint_nfe	AN_nfe
INFO/AC_joint_sas	AC_sas
INFO/AN_joint_sas	AN_sas
INFO/AC_joint_remaining	AC_remaining
INFO/AN_joint_remaining	AN_remaining
INFO/AC_joint_mid	AC_mid
INFO/AN_joint_mid	AN_mid
INFO/AC_joint_XY	AC_XY
INFO/AN_joint_XY	AN_XY
INFO/AF_joint_XY	AF_XY
RENAME

echo "Wrote $RENAME_FILE"

JOBIDS=()
for CHR in $CHROMOSOMES; do
    INPUT_VCF="$INPUT_DIR/gnomad.joint.v4.1.sites.chr${CHR}.vcf.bgz"
    OUTPUT_VCF="$OUTPUT_DIR/gnomad4.1_chr${CHR}.stripped.vcf.gz"

    SLURM_SCRIPT="$OUTPUT_DIR/strip_chr${CHR}.slurm"

    # Build PATH export for SLURM script
    EXTRA_PATH=$(IFS=:; echo "${TOOL_PATHS[*]}")

    cat > "$SLURM_SCRIPT" <<SLURM
#!/bin/bash
#SBATCH --job-name=gnomad41_chr${CHR}
#SBATCH --output=${OUTPUT_DIR}/gnomad41_chr${CHR}_%j.out
#SBATCH --error=${OUTPUT_DIR}/gnomad41_chr${CHR}_%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=${PARTITION}
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail
export PATH="${EXTRA_PATH}:\${PATH}"
echo "Processing chr${CHR} - \$(date)"

bcftools annotate \\
    --remove '${KEEP_STR}' \\
    --rename-annots ${RENAME_FILE} \\
    ${INPUT_VCF} \\
    | sed 's/,Number=A,/,Number=1,/' \\
    | bcftools view -O z -o ${OUTPUT_VCF}

tabix -p vcf ${OUTPUT_VCF}

echo "Done chr${CHR} - \$(date)"
SLURM

    if [ "$DRY_RUN" = "--dry-run" ]; then
        echo "  [dry-run] Would submit $SLURM_SCRIPT"
    else
        JID=$(sbatch --parsable "$SLURM_SCRIPT")
        echo "  Submitted chr${CHR} -> job $JID"
        JOBIDS+=($JID)
    fi
done

echo ""
if [ "$DRY_RUN" = "--dry-run" ]; then
    echo "Dry run complete. 24 scripts written to $OUTPUT_DIR"
    echo "Review a script, then re-run without --dry-run to submit."
else
    echo "Submitted ${#JOBIDS[@]} jobs"
    echo "Monitor with: squeue -u $USER"
fi
