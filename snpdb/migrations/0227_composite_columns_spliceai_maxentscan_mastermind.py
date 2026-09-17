from django.db import migrations

from library.django_utils import bulk_insert_class_data

_KEEP_EVERYTHING = "All columns"

# The visible column each composite cell sits on, the partners it now draws (which leave the system
# collections and stay in the catalogue and 'All columns'), and the column to sit after when a
# collection shows none of the group
_COMPOSITES = [
    ("spliceai_max_ds", ["spliceai_pred_ds_ag", "spliceai_pred_ds_al", "spliceai_pred_ds_dg",
                         "spliceai_pred_ds_dl", "spliceai_pred_dp_ag", "spliceai_pred_dp_al",
                         "spliceai_pred_dp_dg", "spliceai_pred_dp_dl"], None),
    ("maxentscan_percent_diff_ref", ["maxentscan_ref", "maxentscan_alt", "maxentscan_diff"], None),
    ("mastermind_count_1_cdna", ["mastermind_count_2_cdna_prot", "mastermind_count_3_aa_change",
                                 "mastermind_mmid3"], None),
    ("predictions_num_pathogenic", ["predictions_num_benign"], "splice_region"),
    ("total_db_het", ["total_db_hom", "total_db_ref", "total_db_unk"], None),
]

_SPLICEAI_MAX_DS_COLUMN = {
    'grid_column_name': 'spliceai_max_ds',
    'variant_column': 'variantannotation__spliceai_max_ds',
    'annotation_level': 'V',
    'width': 80,
    'label': 'SpliceAI',
    'description': "Highest of the four <a href='https://pubmed.ncbi.nlm.nih.gov/30661751/'>SpliceAI</a> "
                   "delta scores, with a dot at the 0.2 / 0.5 / 0.8 cutoffs. Hover for each prediction "
                   "and its position offset",
    'model_field': True,
    'queryset_field': True,
}

_COLUMN_UPDATES = {
    'maxentscan_percent_diff_ref': (
        90,
        "<a href='http://hollywood.mit.edu/burgelab/maxent/Xmaxentscan_scoreseq.html'>MaxEntScan</a> "
        "change in the 5' splice site score as a percentage of the reference. Negative means a weakened "
        "site. Hover for the reference and alternate scores"),
    'mastermind_count_1_cdna': (
        80,
        "<a href='https://mastermind.genomenon.com/'>Mastermind</a> article counts at increasing "
        "granularity - cDNA, cDNA + protein, amino acid change. Click to open the article list. "
        "Sorts by the cDNA count"),
    'predictions_num_pathogenic': (
        80,
        "Pathogenicity prediction tools that made a call, damaging first - one segment each. "
        "Sorts by the number of damaging calls"),
    'total_db_het': (
        70,
        "Times this variant has been seen in this database - hom &#146; het. Hover for the reference "
        "and unknown counts. Sorts by het count"),
}


def _add_columns(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", [_SPLICEAI_MAX_DS_COLUMN])])
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    for column_id, (width, description) in _COLUMN_UPDATES.items():
        VariantGridColumn.objects.filter(pk=column_id).update(width=width, description=description)


def _update_system_collections(apps, _schema_editor):
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    for ccc in CustomColumnsCollection.objects.filter(user__isnull=True):
        cc_qs = CustomColumn.objects.filter(custom_columns_collection=ccc)
        ordered = list(cc_qs.order_by("sort_order").values_list("column_id", flat=True))
        if ccc.name == _KEEP_EVERYTHING:
            if "spliceai_max_ds" not in ordered:
                ordered.insert(ordered.index("spliceai_pred_ds_ag"), "spliceai_max_ds")
        else:
            for visible, partners, after in _COMPOSITES:
                group = [c for c in ordered if c == visible or c in partners]
                if group:
                    position = ordered.index(group[0])
                elif after in ordered:
                    position = ordered.index(after) + 1
                else:
                    continue
                ordered = [c for c in ordered if c not in partners and c != visible]
                ordered.insert(position, visible)

        for column_id in set(ordered).difference(cc_qs.values_list("column_id", flat=True)):
            CustomColumn.objects.get_or_create(custom_columns_collection=ccc, column_id=column_id,
                                               defaults={"sort_order": 0})
        cc_qs.exclude(column_id__in=ordered).delete()
        for sort_order, column_id in enumerate(ordered):
            cc_qs.filter(column_id=column_id).update(sort_order=sort_order)

        # Historical models skip CustomColumn.save()/delete(), which is what normally bumps the version
        # the node grid definition cache is keyed on
        ccc.version_id += 1
        ccc.save()


def _remove_columns(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="spliceai_max_ds").delete()  # cascades to CustomColumn


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0226_gnomad_af_composite_column'),
    ]

    operations = [
        migrations.RunPython(_add_columns, reverse_code=_remove_columns),
        migrations.RunPython(_update_system_collections, reverse_code=migrations.RunPython.noop),
    ]
