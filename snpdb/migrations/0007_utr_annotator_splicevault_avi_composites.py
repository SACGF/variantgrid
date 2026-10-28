from django.db import migrations

from library.django_utils.seed_data import insert_seed_data

_UTR_ANNOTATOR = "<a href='https://github.com/Ensembl/UTRannotator'>UTRAnnotator</a>"
_SPLICEVAULT = "<a href='https://pubmed.ncbi.nlm.nih.gov/36747048'>SpliceVault</a>"
_AVI = "<a href='https://deepmind.google.com/science/alphagenome'>AlphaGenome Variant Impact</a>"

_COMPOSITES = [
    {'grid_column_name': 'utr_annotator',
     'variant_column': 'utr_annotator',
     'annotation_level': 'T',
     'label': "5' UTR (UTRAnnotator)",
     'description': (f"{_UTR_ANNOTATOR} effect of a 5' UTR variant on upstream ORFs - uAUG gained/lost, uORF stop "
                     "codon gained/lost or uORF frameshift. Hover for each effect's detail (Kozak strength, distance "
                     "to the CDS) and the counts of existing uORFs and of out-of-frame and in-frame oORFs"),
     'model_field': False,
     'queryset_field': False},
    {'grid_column_name': 'splicevault',
     'variant_column': 'splicevault',
     'annotation_level': 'T',
     'label': 'SpliceVault',
     'description': (f"{_SPLICEVAULT} share of the top RNA-seq mis-splicing events that cause a frameshift, at the "
                     "splice site SpliceAI predicts lost, eg Frameshift:2/4 - 3 or more of 4 in-frame events means "
                     "the site is likely splice-rescued. Hover for the site type, SpliceAI delta, GTEx sample "
                     "support, site position and the top events. MANE transcripts only"),
     'model_field': False,
     'queryset_field': False},
    {'grid_column_name': 'avi',
     'variant_column': 'avi',
     'annotation_level': 'V',
     'label': 'AVI',
     'description': (f"{_AVI} PHRED-scaled score - 20 is the top 1% most deleterious substitutions genome-wide - "
                     "with the raw score on hover. SNVs only"),
     'model_field': False,
     'queryset_field': False},
]

# Composite -> (member, in_sort_menu), headline first
_MEMBERS = {
    'utr_annotator': [
        ('utr_annotator_consequence', True),
        ('utr_annotator_annotation', False),
        ('utr_annotator_existing_uorfs', True),
        ('utr_annotator_existing_out_of_frame_oorfs', True),
        ('utr_annotator_existing_inframe_oorfs', True),
    ],
    'splicevault': [
        ('splicevault_out_of_frame_events', True),
        ('splicevault_site_type', True),
        ('splicevault_spliceai_delta', True),
        ('splicevault_site_sample_count', True),
        ('splicevault_site_max_depth', True),
        ('splicevault_site_pos', False),
        ('splicevault_top_events', False),
    ],
    'avi': [
        ('avi_phred', True),
        ('avi_raw', True),
    ],
}

_COMPOSITE_NAMES = [c['grid_column_name'] for c in _COMPOSITES]

# The order snpdb 0003 appended the members to "All columns" in
_MEMBERS_IN_ALL_COLUMNS_ORDER = [
    'utr_annotator_consequence', 'utr_annotator_annotation', 'utr_annotator_existing_uorfs',
    'utr_annotator_existing_out_of_frame_oorfs', 'utr_annotator_existing_inframe_oorfs',
    'splicevault_top_events', 'splicevault_site_sample_count', 'splicevault_site_max_depth',
    'splicevault_out_of_frame_events', 'splicevault_site_pos', 'splicevault_site_type',
    'splicevault_spliceai_delta', 'avi_phred', 'avi_raw',
]


def _new_composites(apps, schema_editor):
    composite_members = [
        {'composite_id': composite, 'column_id': column, 'sort_order': i, 'in_sort_menu': in_sort_menu}
        for composite, members in _MEMBERS.items()
        for i, (column, in_sort_menu) in enumerate(members)
    ]
    insert_seed_data(apps, schema_editor, [
        ("snpdb.VariantGridColumn", _COMPOSITES),
        ("snpdb.CompositeColumnMember", composite_members),
    ])


def _reverse_new_composites(apps, _schema_editor):
    apps.get_model("snpdb", "VariantGridColumn").objects.filter(grid_column_name__in=_COMPOSITE_NAMES).delete()


def _replace_in_all_columns(apps, old_column_ids: list[str], new_column_ids: list[str]):
    """ Swaps old_column_ids out of the public "All columns" for new_column_ids, placed where the first
        old one was (at the end if none were there) """
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    all_columns = CustomColumnsCollection.objects.filter(name='All columns', user__isnull=True).first()
    if all_columns is None:
        return
    rows = list(all_columns.customcolumn_set.order_by("sort_order", "pk"))
    present = {r.column_id for r in rows}
    old_rows = [r for r in rows if r.column_id in old_column_ids]
    insert_at = rows.index(old_rows[0]) if old_rows else len(rows)

    new_rows = [CustomColumn(custom_columns_collection=all_columns, column_id=column_id, sort_order=0)
                for column_id in new_column_ids if column_id not in present]
    kept_before = rows[:insert_at]
    kept_after = [r for r in rows[insert_at:] if r.column_id not in old_column_ids]
    ordered = kept_before + new_rows + kept_after
    for sort_order, row in enumerate(ordered):
        row.sort_order = sort_order

    CustomColumn.objects.filter(pk__in=[r.pk for r in old_rows]).delete()
    CustomColumn.objects.bulk_update(kept_before + kept_after, ["sort_order"])
    CustomColumn.objects.bulk_create(new_rows)
    # CustomColumn.save()/delete() do this for the grid caches, but the bulk calls skip it
    all_columns.version_id += 1
    all_columns.save(update_fields=["version_id"])


def _composites_into_all_columns(apps, _schema_editor):
    _replace_in_all_columns(apps, _MEMBERS_IN_ALL_COLUMNS_ORDER, _COMPOSITE_NAMES)


def _reverse_composites_into_all_columns(apps, _schema_editor):
    _replace_in_all_columns(apps, _COMPOSITE_NAMES, _MEMBERS_IN_ALL_COLUMNS_ORDER)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0006_one_off_dedupe_cohort_genotypes'),
    ]

    operations = [
        migrations.RunPython(_new_composites, reverse_code=_reverse_new_composites),
        migrations.RunPython(_composites_into_all_columns, reverse_code=_reverse_composites_into_all_columns),
    ]
