#!/bin/bash
#
# Shared config for stripping gnomAD v4.1 joint VCFs:
#   - KEEP_FIELDS: v4.1 INFO field names we retain
#   - write_rename_file: writes the bcftools --rename-annots TSV
#     mapping v4.1 _joint names to legacy (v4.0-style) names
#
# Source this file; do not execute directly.
#

# v4.1 source fields to keep (everything else is stripped).
# Mapping to legacy names is applied by write_rename_file below.
KEEP_FIELDS=(
    # Core counts (global)
    AC_joint AN_joint AF_joint nhomalt_joint
    # Filtering allele frequencies
    faf95_joint faf99_joint fafmax_faf95_max_joint fafmax_faf99_max_joint
    # grpmax
    AF_grpmax_joint AC_grpmax_joint AN_grpmax_joint grpmax_joint
    # Per-population AC/AN/AF
    AC_joint_afr AN_joint_afr AF_joint_afr
    AC_joint_amr AN_joint_amr AF_joint_amr
    AC_joint_asj AN_joint_asj AF_joint_asj
    AC_joint_eas AN_joint_eas AF_joint_eas
    AC_joint_fin AN_joint_fin AF_joint_fin
    AC_joint_mid AN_joint_mid AF_joint_mid
    AC_joint_nfe AN_joint_nfe AF_joint_nfe
    AC_joint_remaining AN_joint_remaining AF_joint_remaining
    AC_joint_sas AN_joint_sas AF_joint_sas
    # ChrX sex-specific
    AC_joint_XY AN_joint_XY AF_joint_XY
    # ChrX/chrY non-PAR flag (only present on sex chromosomes; Flag type so
    # no _joint suffix — used to derive gnomad_hemi_count downstream)
    non_par
)

# Build the bcftools --remove expression: "^INFO/foo,INFO/bar,..."
# (caret = complement = keep only these, strip everything else)
build_keep_expr() {
    local s
    s=$(printf "INFO/%s," "${KEEP_FIELDS[@]}")
    echo "^${s%,}"
}

# Write the bcftools --rename-annots TSV to $1.
# Maps v4.1 _joint names -> legacy v4.0-style names.
write_rename_file() {
    local out="$1"
    cat > "$out" <<'RENAME'
INFO/AC_joint	AC
INFO/AN_joint	AN
INFO/AF_joint	AF
INFO/nhomalt_joint	nhomalt
INFO/faf95_joint	faf95
INFO/faf99_joint	faf99
INFO/fafmax_faf95_max_joint	fafmax_faf95_max
INFO/fafmax_faf99_max_joint	fafmax_faf99_max
INFO/AF_grpmax_joint	AF_grpmax
INFO/AC_grpmax_joint	AC_grpmax
INFO/AN_grpmax_joint	AN_grpmax
INFO/grpmax_joint	grpmax
INFO/AC_joint_afr	AC_afr
INFO/AN_joint_afr	AN_afr
INFO/AF_joint_afr	AF_afr
INFO/AC_joint_amr	AC_amr
INFO/AN_joint_amr	AN_amr
INFO/AF_joint_amr	AF_amr
INFO/AC_joint_asj	AC_asj
INFO/AN_joint_asj	AN_asj
INFO/AF_joint_asj	AF_asj
INFO/AC_joint_eas	AC_eas
INFO/AN_joint_eas	AN_eas
INFO/AF_joint_eas	AF_eas
INFO/AC_joint_fin	AC_fin
INFO/AN_joint_fin	AN_fin
INFO/AF_joint_fin	AF_fin
INFO/AC_joint_mid	AC_mid
INFO/AN_joint_mid	AN_mid
INFO/AF_joint_mid	AF_mid
INFO/AC_joint_nfe	AC_nfe
INFO/AN_joint_nfe	AN_nfe
INFO/AF_joint_nfe	AF_nfe
INFO/AC_joint_remaining	AC_remaining
INFO/AN_joint_remaining	AN_remaining
INFO/AF_joint_remaining	AF_remaining
INFO/AC_joint_sas	AC_sas
INFO/AN_joint_sas	AN_sas
INFO/AF_joint_sas	AF_sas
INFO/AC_joint_XY	AC_XY
INFO/AN_joint_XY	AN_XY
INFO/AF_joint_XY	AF_XY
RENAME
}
