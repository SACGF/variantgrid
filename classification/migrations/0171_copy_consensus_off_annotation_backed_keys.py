from django.db import migrations

# sacgf/variantgrid#1713 - annotation-derived keys that were left copyable when their siblings weren't

_TRANSCRIPT_LEVEL_KEYS = [
    'alphamissense',
    'bayesdel',
    'clinpred',
    'early_stop_codon',
    'interpro_protein_domain',
    'mavedb',
    'metarnn',
    'metarnn_pred',
    'mpc',
    'mutpred2',
    'nmd_escaping_variant',
    'pfam_protein_domain',
    'primateai',
    'varity_er',
    'varity_r',
    'vest',
]
"""
Copy consensus is scoped to the allele, so the source record can be curated on a different transcript -
these values then describe a transcript the target record isn't using, and autopopulate leaves them in
place when it has nothing to overwrite them with. revel, cadd, sift, polyphen2, fathmm, gerp,
mutation_taster, mutation_assessor and molecular_consequence are all transcript level and already False.
"""

_VARIANT_LEVEL_KEYS = [
    'alphamissense_rankscore',
    'aloft',
    'gnomad_allele_count',
    'gnomad_allele_number',
    'gnomad_hemi_count',
    'gnomad_hom_count',
    'gnomad_popmax_allele_count',
    'gnomad_popmax_allele_number',
]
"""
Counts carried over from whichever gnomAD version annotated the source record, alongside gnomad_af,
gnomad_popmax_af, topmed_af, uk10k_af and 1000_genomes_af which are already False. aloft is built in
code rather than from a variantgrid_column, as are its already False siblings phastcons, phylop and
spliceai.
"""

_GENE_LEVEL_KEYS = [
    'essential_gene_crispr',
    'essential_gene_crispr2',
    'essential_gene_gene_trap',
    'gene_damage_index_score',
    'gene_indispensability_score',
    'ghis',
    'gnomad_oe_lof',
    'gnomad_pli',
    'gnomad_pnull',
    'gnomad_prec',
    'hipred_score',
    'phi',
    'prec',
]
"""
Gene constraint metrics that autopopulate fills for every classification from the same tables. Where it
has a value the copy is redundant, and where it doesn't the gap is the honest answer rather than
something to fill from another record. loftool is gene level and already False.
"""

_KEYS = _TRANSCRIPT_LEVEL_KEYS + _VARIANT_LEVEL_KEYS + _GENE_LEVEL_KEYS


def _copy_consensus_off(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key__in=_KEYS).update(copy_consensus=False)


def _copy_consensus_on(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key__in=_KEYS).update(copy_consensus=True)


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0170_copy_consensus_off_patient_and_report_keys'),
    ]

    operations = [
        migrations.RunPython(_copy_consensus_off, reverse_code=_copy_consensus_on),
    ]
