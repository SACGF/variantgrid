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
    # Witness fields used to derive a synthetic non_par flag (see
    # stamp_non_par_from_witness below). v4.1 dropped the top-level non_par
    # INFO flag that v4.0 had, but the per-sex FAF fields are populated only
    # on non-PAR records (their header descriptions read "...in non-PAR
    # regions of sex chromosomes only"), so their presence is the marker.
    # These are stripped again by stamp_non_par_from_witness after the flag
    # has been set, so they never reach the final VCF.
    faf99_joint_XX faf99_joint_XY
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

# Write the non_par INFO header line to $1, suitable for use with
# `bcftools annotate --header-lines`. Adding the line via bcftools (instead
# of via awk) lets the header land for free during an existing streaming
# stage, which matters on autosomes where no body rewrite is otherwise
# needed.
write_non_par_header_file() {
    local out="$1"
    cat > "$out" <<'HDR'
##INFO=<ID=non_par,Number=0,Type=Flag,Description="Variant (on sex chromosome) falls outside a pseudoautosomal region">
HDR
}

# Stream filter: read VCF on stdin, write VCF on stdout. Used only on
# chrX/chrY (autosomes don't need it).
#   - on data rows, sets non_par in INFO when any witness FAF field
#     (faf99_joint_XX / faf99_joint_XY) is present
#   - strips the witness fields from INFO afterwards (they exist only to
#     drive this stage)
# The non_par INFO header line is added upstream via
# `bcftools annotate --header-lines`, so this filter does not touch the
# header — it just streams meta rows through unchanged.
stamp_non_par_from_witness() {
    awk '
    BEGIN { FS = OFS = "\t" }
    /^#/ { print; next }
    {
        n = split($8, kv, ";")
        is_non_par = 0
        out = ""
        for (i = 1; i <= n; i++) {
            f = kv[i]
            tag = f
            sub(/=.*/, "", tag)
            if (tag == "faf99_joint_XX" || tag == "faf99_joint_XY") {
                is_non_par = 1
                continue
            }
            out = (out == "" ? f : out ";" f)
        }
        if (is_non_par) out = (out == "" ? "non_par" : out ";non_par")
        $8 = out
        print
    }'
}
