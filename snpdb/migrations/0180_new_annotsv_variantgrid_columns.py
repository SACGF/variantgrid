# VariantGridColumn registration for AnnotSV fields - #1040 / #1533 follow-up.
# Backfills the 14 columns shipped without VariantGridColumn entries in #1533,
# and adds the 10 new AnnotSV fields landed alongside columns_version 4.

from django.db import migrations
from django.db.models import F, Max

from library.django_utils import bulk_insert_class_data


_ANNOTSV_LINK = "<a href='https://lbgi.fr/AnnotSV/' target='_blank'>AnnotSV</a>"

_ANNOTSV_COLUMNS = [
    # (grid_column_name, label, description)

    # ---- pre-existing #1533 fields (backfill VariantGridColumn) ----
    ("annotsv_acmg_class",
     "AnnotSV ACMG class",
     f"{_ANNOTSV_LINK} ACMG-style ranking class for SVs (1=benign, 2=likely benign, 3=VUS, 4=likely pathogenic, 5=pathogenic)."),
    ("annotsv_acmg_score",
     "AnnotSV ACMG score",
     f"{_ANNOTSV_LINK} ACMG-style ranking score: sum of points across triggered rules."),
    ("annotsv_re_gene",
     "Regulatory element gene",
     "Regulatory element gene(s) overlapping the SV, per AnnotSV's regulatory-element bundle (RE_gene)."),
    ("annotsv_repeat_type_left",
     "Repeat type (left)",
     "Repeat type at the SV left breakpoint, per UCSC RepeatMasker (Repeat_type_left)."),
    ("annotsv_repeat_type_right",
     "Repeat type (right)",
     "Repeat type at the SV right breakpoint, per UCSC RepeatMasker (Repeat_type_right)."),
    ("annotsv_segdup_left",
     "SegDup (left)",
     "UCSC segmental duplication feature overlapping the SV left breakpoint (SegDup_left)."),
    ("annotsv_segdup_right",
     "SegDup (right)",
     "UCSC segmental duplication feature overlapping the SV right breakpoint (SegDup_right)."),
    ("annotsv_encode_blacklist_left",
     "ENCODE blacklist (left)",
     "ENCODE blacklist region overlapping the SV left breakpoint."),
    ("annotsv_encode_blacklist_right",
     "ENCODE blacklist (right)",
     "ENCODE blacklist region overlapping the SV right breakpoint."),
    ("annotsv_encode_blacklist_characteristics_left",
     "ENCODE blacklist characteristics (left)",
     "ENCODE blacklist characteristic at the SV left breakpoint (e.g. Low Mappability, High Signal Region)."),
    ("annotsv_encode_blacklist_characteristics_right",
     "ENCODE blacklist characteristics (right)",
     "ENCODE blacklist characteristic at the SV right breakpoint."),
    ("annotsv_b_gain_af_max",
     "AnnotSV benign gain AF max",
     f"Max population AF across {_ANNOTSV_LINK}'s benign-SV gain (duplication/insertion-gain) sources (B_gain_AFmax)."),
    ("annotsv_b_loss_af_max",
     "AnnotSV benign loss AF max",
     f"Max population AF across {_ANNOTSV_LINK}'s benign-SV loss (deletion) sources (B_loss_AFmax)."),
    ("annotsv_b_ins_af_max",
     "AnnotSV benign ins AF max",
     f"Max population AF across {_ANNOTSV_LINK}'s benign-SV insertion sources (B_ins_AFmax)."),
    ("annotsv_b_inv_af_max",
     "AnnotSV benign inv AF max",
     f"Max population AF across {_ANNOTSV_LINK}'s benign-SV inversion sources (B_inv_AFmax)."),

    # ---- new in this change ----
    ("annotsv_acmg_criteria",
     "AnnotSV ACMG criteria",
     f"Triggered {_ANNOTSV_LINK} ACMG-SV rules (e.g. 1A, 2H) that summed into the ranking score."),
    ("annotsv_frameshift",
     "AnnotSV frameshift",
     "Whether the SV introduces a frameshift in any overlapped transcript."),
    ("annotsv_exons_spanned",
     "AnnotSV exons spanned",
     "Number of exons fully spanned by the SV."),
    ("annotsv_dist_nearest_ss",
     "Distance to nearest splice site",
     "Distance in bp from the SV breakpoint to the nearest splice site (Dist_nearest_SS)."),
    ("annotsv_nearest_ss_type",
     "Nearest splice site type",
     "Type of the nearest splice site (5' donor or 3' acceptor)."),
    ("annotsv_omim_inheritance",
     "OMIM inheritance",
     "OMIM inheritance pattern(s) for genes overlapped by the SV (AD, AR, XL, etc.)."),
    ("annotsv_omim_morbid",
     "OMIM morbid",
     "Whether any gene overlapped by the SV is listed in the OMIM morbid map."),
    ("annotsv_omim_phenotype",
     "OMIM phenotype",
     "OMIM phenotype text for genes overlapped by the SV."),
    ("annotsv_omim_id",
     "OMIM ID",
     "OMIM MIM number(s) for genes overlapped by the SV."),
    ("annotsv_pathogenic_overlaps",
     "Pathogenic SV overlaps",
     f"Per-event-type (gain/loss/ins/inv) summary of pathogenic SV overlaps from {_ANNOTSV_LINK}'s pathogenic-SV bundle: source (ClinVar / dbVar / ClinGen / OMIM-morbid), reported phenotype, HPO terms, and reference SV coordinates."),
]


def _add_annotsv_columns(apps, _schema_editor):
    rows = [
        {
            "grid_column_name": name,
            "variant_column": f"variantannotation__{name}",
            "annotation_level": "V",
            "width": None,
            "label": label,
            "description": description,
            "model_field": True,
            "queryset_field": True,
        }
        for name, label, description in _ANNOTSV_COLUMNS
    ]
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", rows)])

    # Append to "All columns" after the SV-overlap block. This is the first
    # time any annotsv_ column gets registered, so anchor off the existing
    # gnomad_sv_overlap entries.
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")
    all_columns = CustomColumnsCollection.objects.get(name="All columns")
    sv_qs = all_columns.customcolumn_set.filter(
        column__grid_column_name__startswith="gnomad_sv_overlap",
        column__annotation_level="V",
    )
    sort_order_max = sv_qs.aggregate(Max("sort_order"))["sort_order__max"] or 0
    all_columns.customcolumn_set.filter(sort_order__gt=sort_order_max).update(
        sort_order=F("sort_order") + len(_ANNOTSV_COLUMNS)
    )
    for i, (name, _label, _desc) in enumerate(_ANNOTSV_COLUMNS, start=1):
        CustomColumn.objects.create(
            custom_columns_collection=all_columns,
            sort_order=sort_order_max + i,
            column_id=name,
        )


def _remove_annotsv_columns(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(
        grid_column_name__in=[name for name, _label, _desc in _ANNOTSV_COLUMNS]
    ).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0179_cohort_genotype_stats_backfill"),
    ]

    operations = [
        migrations.RunPython(_add_annotsv_columns, reverse_code=_remove_annotsv_columns),
    ]
