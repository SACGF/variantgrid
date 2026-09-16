from django.db import migrations

# sacgf/variantgrid#1714 - gene level content travels to every variant in the gene, and germline
# concepts stay out of somatic records

_GENE_KEYS = [
    'condition_incidence',
    'gene_constraint',
    'gene_disease_validity',
    'gene_penetrance',
    'h_summary',
    'mechanism_of_disease',
    'mode_of_inheritance',
    'pubmed_gene_search_count',
    'variant_penetrance',
]
"""
The curated gene content - identical for every variant in the gene, and what gene level reuse is for.
The annotation-derived gene metrics of the same category went to copy_consensus=False in #1713 and stay
NONE, they arrive from autopopulate. condition stays at ALLELE - for somatic it is the tumour type, so
reuse is an explicit human decision rather than a silent pre-fill. literature and search_terms stay at
ALLELE until literature is split into gene and variant level content (variantgrid_private#1102).
"""

_GENE_KEYS_PREVIOUSLY_UNCOPIED = [
    'disease_onset',
]
"""
Curated gene content like the above, but uncopyable until now - gene level is where it belongs.
"""

_GERMLINE_KEYS = [
    'a_other',
    'a_summary',
    'condition_incidence',
    'd_other',
    'd_summary',
    'denovo_points',
    'gene_penetrance',
    'match_maker_exchange',
    'mode_of_inheritance',
    'proband_count',
    's_other',
    's_summary',
    'segregation',
    'segregation_affectedcarriers',
    'segregation_affectednoncarriers',
    'segregation_bayes',
    'segregation_lod',
    'segregation_meioses',
    'segregation_unaffectedcarriers',
    'variant_penetrance',
]
"""
Segregation, de novo and allelic data are germline concepts with no namespace to filter on, so nothing
else stops them being copied into a somatic record. The ones that are also gene level carry both -
copying then needs the gene to match and the record to be germline.
"""


def _set_copy_scope_and_allele_origin(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key__in=_GENE_KEYS + _GENE_KEYS_PREVIOUSLY_UNCOPIED).update(copy_scope="G")
    EvidenceKey.objects.filter(key__in=_GERMLINE_KEYS).update(copy_allele_origin="G")


def _unset_copy_scope_and_allele_origin(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key__in=_GENE_KEYS).update(copy_scope="A")
    EvidenceKey.objects.filter(key__in=_GENE_KEYS_PREVIOUSLY_UNCOPIED).update(copy_scope="N")
    EvidenceKey.objects.filter(key__in=_GERMLINE_KEYS).update(copy_allele_origin="A")


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0176_evidence_key_copy_scope'),
    ]

    operations = [
        migrations.RunPython(_set_copy_scope_and_allele_origin, reverse_code=_unset_copy_scope_and_allele_origin),
    ]
