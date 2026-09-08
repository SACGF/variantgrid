"""
Two more composite cells - the non-gnomAD population frequencies and the UniProt gene annotation -
the conservation cell's description now that it draws a dot per score, and widths for the composite
cells measured against what they actually draw.

Reverse is lossy in the same way as snpdb/0238: the composite rows go and a collection keeps only
whatever members it still had.
"""
from django.db import migrations

from library.django_utils import bulk_insert_class_data
from library.django_utils.composite_columns import collapse_into_composite

_VARIANT = 'V'
_UNIPROT = 'U'

_COMPOSITES = [
    {'pk': 'pop_freq_other', 'level': _VARIANT, 'width': 110, 'label': 'Pop Freq Other',
     'description': 'The highest allele frequency among the population databases other than gnomAD - '
                    '1000 Genomes, UK10K and TOPMed - with which one it came from. Hover for each. '
                    'Sorts by the 1000 Genomes frequency',
     'members': ['af_1kg', 'af_uk10k', 'topmed_af']},
    {'pk': 'uniprot', 'level': _UNIPROT, 'width': 300, 'label': 'UniProt',
     'description': '<a href="http://www.uniprot.org">UniProtKB</a> entry for the gene - the accession, '
                    'linked to UniProt, and the function summary. Hover for the pathway, tissue '
                    'specificity and Reactome pathways',
     'members': ['uniprot_id', 'function_from_uniprotkb', 'pathway_from_uniprotkb',
                 'tissue_specificity_from_uniprotkb', 'uniprot_reactome'],
     'detail_only': ['function_from_uniprotkb', 'pathway_from_uniprotkb',
                     'tissue_specificity_from_uniprotkb', 'uniprot_reactome']},
]

_CONSERVATION_DESCRIPTION = (
    'One dot per conservation score - PhyloP, PhastCons and GERP++ - filled where the score reaches '
    'its conserved threshold, hollow where this position has no score. Hover for each score against '
    'its range. Sorts by PhyloP 100 way vertebrate'
)

# Measured over 8000 random annotated variants at the 95th percentile of what the cell leads with.
# gnomAD fits once the popmax population is drawn as its three letter code
_WIDTH_UPDATES = {
    'gnomad': 170,
    'conservation': 90,
    'aloft': 90,
}

_NEW_COMPOSITE_COLUMNS = [
    {'grid_column_name': c['pk'], 'variant_column': c['pk'], 'annotation_level': c['level'],
     'width': c['width'], 'label': c['label'], 'description': c['description'],
     'model_field': False, 'queryset_field': False}
    for c in _COMPOSITES
]


def _add_composite_columns(apps, _schema_editor):
    bulk_insert_class_data(apps, "snpdb", [("VariantGridColumn", _NEW_COMPOSITE_COLUMNS)])
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    VariantGridColumn.objects.filter(pk="conservation").update(description=_CONSERVATION_DESCRIPTION)
    for column_id, width in _WIDTH_UPDATES.items():
        VariantGridColumn.objects.filter(pk=column_id).update(width=width)


def _add_members(apps, _schema_editor):
    CompositeColumnMember = apps.get_model("snpdb", "CompositeColumnMember")
    members = []
    for composite in _COMPOSITES:
        detail_only = set(composite.get('detail_only', []))
        for sort_order, column_id in enumerate(composite['members']):
            members.append(CompositeColumnMember(composite_id=composite['pk'], column_id=column_id,
                                                 sort_order=sort_order,
                                                 in_sort_menu=column_id not in detail_only))
    CompositeColumnMember.objects.bulk_create(members)


def _collapse_collections(apps, _schema_editor):
    for composite in _COMPOSITES:
        collapse_into_composite(apps, composite['pk'])


def _remove_composite_columns(apps, _schema_editor):
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    # Cascades to the members and to the CustomColumns showing the composite
    VariantGridColumn.objects.filter(pk__in=[c['pk'] for c in _COMPOSITES]).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0240_open_targets_l2g_scores_backfill_note'),
    ]

    operations = [
        migrations.RunPython(_add_composite_columns, reverse_code=_remove_composite_columns),
        migrations.RunPython(_add_members, reverse_code=migrations.RunPython.noop),
        migrations.RunPython(_collapse_collections, reverse_code=migrations.RunPython.noop),
    ]
