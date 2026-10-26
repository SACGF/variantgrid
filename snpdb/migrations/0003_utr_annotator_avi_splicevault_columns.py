from django.db import migrations
from django.db.models import Max

from library.django_utils.seed_data import insert_seed_data

_TRANSCRIPT = 'T'
_VARIANT = 'V'
_UTR_ANNOTATOR = 'https://github.com/Ensembl/UTRannotator'
_SPLICEVAULT = 'https://pubmed.ncbi.nlm.nih.gov/36747048'
_AVI = 'https://deepmind.google.com/science/alphagenome'

_NEW_COLUMNS = [
    {'grid_column_name': 'utr_annotator_consequence',
     'variant_column': 'variantannotation__utr_annotator_consequence',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR uORF consequence",
     'description': ("UTRAnnotator effect of a 5' UTR variant on upstream ORFs: uAUG gained/lost, uORF stop "
                     f"codon gained/lost, or uORF frameshift. {_UTR_ANNOTATOR}")},
    {'grid_column_name': 'utr_annotator_annotation',
     'variant_column': 'variantannotation__utr_annotator_annotation',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR uORF details",
     'description': ("UTRAnnotator detail of each uORF effect (&-joined), eg Kozak context and strength, "
                     "uORF type, distance to the CDS and Evidence (whether the uORF data file has it). "
                     f"{_UTR_ANNOTATOR}")},
    {'grid_column_name': 'utr_annotator_existing_uorfs',
     'variant_column': 'variantannotation__utr_annotator_existing_uorfs',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR existing uORFs",
     'description': f"UTRAnnotator count of uORFs with a stop codon within the 5' UTR. {_UTR_ANNOTATOR}"},
    {'grid_column_name': 'utr_annotator_existing_out_of_frame_oorfs',
     'variant_column': 'variantannotation__utr_annotator_existing_out_of_frame_oorfs',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR existing out-of-frame oORFs",
     'description': ("UTRAnnotator count of out-of-frame ORFs in the 5' UTR that overlap the CDS. "
                     f"{_UTR_ANNOTATOR}")},
    {'grid_column_name': 'utr_annotator_existing_inframe_oorfs',
     'variant_column': 'variantannotation__utr_annotator_existing_inframe_oorfs',
     'annotation_level': _TRANSCRIPT,
     'label': "5' UTR existing in-frame oORFs",
     'description': f"UTRAnnotator count of in-frame ORFs in the 5' UTR that overlap the CDS. {_UTR_ANNOTATOR}"},
    {'grid_column_name': 'splicevault_top_events',
     'variant_column': 'variantannotation__splicevault_top_events',
     'annotation_level': _TRANSCRIPT,
     'label': 'SpliceVault top events',
     'description': ("SpliceVault's most common mis-splicing events in RNA-seq (GTEx) at the splice site SpliceAI "
                     "predicts lost, &-joined, each rank:type:transcript_impact:percent_of_samples:frame. Type is "
                     "exon skipping (ES), cryptic donor (CD) or cryptic acceptor (CA). MANE transcripts only. "
                     f"{_SPLICEVAULT}")},
    {'grid_column_name': 'splicevault_site_sample_count',
     'variant_column': 'variantannotation__splicevault_site_sample_count',
     'annotation_level': _TRANSCRIPT,
     'label': 'SpliceVault site sample count',
     'description': f"SpliceVault count of samples with annotated splicing at this splice site. {_SPLICEVAULT}"},
    {'grid_column_name': 'splicevault_site_max_depth',
     'variant_column': 'variantannotation__splicevault_site_max_depth',
     'annotation_level': _TRANSCRIPT,
     'label': 'SpliceVault site max depth',
     'description': ("SpliceVault maximum reads in any one GTEx sample representing annotated splicing at this "
                     f"splice site. {_SPLICEVAULT}")},
    {'grid_column_name': 'splicevault_out_of_frame_events',
     'variant_column': 'variantannotation__splicevault_out_of_frame_events',
     'annotation_level': _TRANSCRIPT,
     'label': 'SpliceVault out-of-frame events',
     'description': ("SpliceVault fraction of the top events that cause a frameshift, eg Frameshift:2/4. Sites with "
                     f"3 or more of 4 in-frame events are likely splice-rescued, not loss of function. {_SPLICEVAULT}")},
    {'grid_column_name': 'splicevault_site_pos',
     'variant_column': 'variantannotation__splicevault_site_pos',
     'annotation_level': _TRANSCRIPT,
     'label': 'SpliceVault site position',
     'description': f"Genomic position of the splice site SpliceAI predicts lost. {_SPLICEVAULT}"},
    {'grid_column_name': 'splicevault_site_type',
     'variant_column': 'variantannotation__splicevault_site_type',
     'annotation_level': _TRANSCRIPT,
     'label': 'SpliceVault site type',
     'description': f"Type of the splice site SpliceAI predicts lost (Donor_loss / Acceptor_loss). {_SPLICEVAULT}"},
    {'grid_column_name': 'splicevault_spliceai_delta',
     'variant_column': 'variantannotation__splicevault_spliceai_delta',
     'annotation_level': _TRANSCRIPT,
     'label': 'SpliceVault SpliceAI delta',
     'description': f"SpliceAI delta score for the lost splice site, as SpliceVault provides it. {_SPLICEVAULT}"},
    {'grid_column_name': 'avi_phred',
     'variant_column': 'variantannotation__avi_phred',
     'annotation_level': _VARIANT,
     'label': 'AVI (phred)',
     'description': ("AlphaGenome Variant Impact (Google DeepMind) PHRED-scaled quantile, 0 to ~80: 20 is the top "
                     f"1% most deleterious substitutions genome-wide. Higher is more deleterious. SNVs only. {_AVI}")},
    {'grid_column_name': 'avi_raw',
     'variant_column': 'variantannotation__avi_raw',
     'annotation_level': _VARIANT,
     'label': 'AVI (raw)',
     'description': ("AlphaGenome Variant Impact (Google DeepMind) raw score, about -1.3 to +4.0. Higher is more "
                     f"deleterious. SNVs only. {_AVI}")},
]


def _new_columns(apps, schema_editor):
    insert_seed_data(apps, schema_editor, [("snpdb.VariantGridColumn", _NEW_COLUMNS)])


def _reverse_new_columns(apps, _schema_editor):
    grid_names = [c['grid_column_name'] for c in _NEW_COLUMNS]
    apps.get_model("snpdb", "VariantGridColumn").objects.filter(grid_column_name__in=grid_names).delete()


def _append_to_all_columns(apps, _schema_editor):
    CustomColumnsCollection = apps.get_model("snpdb", "CustomColumnsCollection")
    CustomColumn = apps.get_model("snpdb", "CustomColumn")

    all_columns = CustomColumnsCollection.objects.filter(name='All columns', user__isnull=True).first()
    if all_columns is None:
        return
    sort_order_max = all_columns.customcolumn_set.aggregate(Max("sort_order"))["sort_order__max"] or 0
    CustomColumn.objects.bulk_create([
        CustomColumn(custom_columns_collection=all_columns, sort_order=sort_order_max + 1 + i,
                     column_id=c['grid_column_name'])
        for i, c in enumerate(_NEW_COLUMNS)
    ])
    # CustomColumn.save() does this for the grid caches, but bulk_create skips it
    all_columns.version_id += 1
    all_columns.save(update_fields=["version_id"])


def _reverse_append_to_all_columns(apps, _schema_editor):
    CustomColumn = apps.get_model("snpdb", "CustomColumn")
    grid_names = [c['grid_column_name'] for c in _NEW_COLUMNS]
    CustomColumn.objects.filter(custom_columns_collection__name='All columns',
                                column_id__in=grid_names).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0004_utr_annotator_avi_splicevault'),
        ('snpdb', '0002_seed_data'),
    ]

    operations = [
        migrations.RunPython(_new_columns, reverse_code=_reverse_new_columns),
        migrations.RunPython(_append_to_all_columns, reverse_code=_reverse_append_to_all_columns),
    ]
