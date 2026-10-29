from django.db import migrations

# SACGF/variantgrid#2563 - evidence keys autopopulated from the AVI, AlphaGenome splicing, SpliceVault and
# UTRAnnotator annotation (#2426 / #2433). Sources and score meanings: #2561

_AVI_SEE = 'https://www.medrxiv.org/content/10.64898/2026.09.16.26363192v1.full'
_ALPHAGENOME_SEE = 'https://www.biorxiv.org/content/10.1101/2025.06.25.661532v1.full'
_SPLICEVAULT_SEE = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9925382/'
_UTR_ANNOTATOR_SEE = 'https://doi.org/10.1093/bioinformatics/btaa783'

_NEW_KEYS = [
    {
        'key': 'avi',
        'label': 'AVI',
        'sub_label': 'Meta-predictor',
        'description': 'AlphaGenome Variant Impact (AVI) PHRED score. <p> Score range: 0 to ~80 <p> A genome-wide rank '
                       'over all SNVs: 10 is the top 10% of predictions, 20 the top 1%. Higher is more deleterious. '
                       'Combines AlphaGenome regulatory predictions with AlphaMissense, VEP truncation features and '
                       'conservation. The authors give no pathogenicity cutoff (they used the top 1% for rare disease '
                       'triage), and it has no ClinGen PP3/BP4 calibration. GRCh38 SNVs only.',
        'examples': [4.2, 23.1],
        'see': _AVI_SEE,
        'value_type': 'L',
        'order': 5,
        'variantgrid_column_id': 'avi_phred',
    },
    {
        'key': 'alphagenome_splicing',
        'label': 'AlphaGenome splicing',
        'sub_label': 'Splicing',
        'description': 'AlphaGenome merged splicing score. <p> Score range: 0 to ~3 <p> The sum of its splice site, '
                       'splice site usage and splice junction (divided by 5) predictions, each 0 to 1. Higher is a '
                       'bigger predicted effect on splicing. No pathogenicity cutoff has been proposed or calibrated. '
                       'GRCh38 SNVs only. Licensed for non-commercial use only, and not for diagnostic procedures.',
        'examples': [0.02, 1.4],
        'see': _ALPHAGENOME_SEE,
        'value_type': 'L',
        'order': 21,
        'variantgrid_column_id': 'alphagenome_splicing',
    },
    {
        'key': 'splicevault_top_events',
        'label': 'SpliceVault top events',
        'sub_label': 'Splicing',
        'description': "SpliceVault's most common mis-splicing events seen in RNA-seq at the splice site SpliceAI "
                       'predicts lost, each rank:type;position;percent of samples;frame. Type is exon skipping (ES, '
                       'position is the exons skipped), cryptic donor (CD) or cryptic acceptor (CA). <p> It predicts the kind of mis-splicing, not '
                       'whether splicing is disrupted. The authors suggest PVS1 at Very Strong only when intron '
                       'retention and all top events are consistent with a null outcome. MANE transcripts only.',
        'examples': ['Top1:CD;-35;31%;Frameshift, Top2:ES;1;5%;Inframe'],
        'see': _SPLICEVAULT_SEE,
        'order': 23,
        'variantgrid_column_id': 'splicevault_top_events',
    },
    {
        'key': 'splicevault_out_of_frame_events',
        'label': 'SpliceVault out-of-frame events',
        'sub_label': 'Splicing',
        'description': 'SpliceVault fraction of the top mis-splicing events that cause a frameshift. All out of '
                       'frame is consistent with a null outcome; 3 or more in-frame events suggest splice rescue '
                       '(PM4 rather than PVS1). MANE transcripts only.',
        'examples': ['Frameshift:4/4', 'Frameshift:1/4'],
        'see': _SPLICEVAULT_SEE,
        'order': 23,
        'variantgrid_column_id': 'splicevault_out_of_frame_events',
    },
    {
        'key': 'splicevault_site_sample_count',
        'label': 'SpliceVault site sample count',
        'sub_label': 'Splicing',
        'description': 'Number of RNA-seq samples with annotated splicing at the splice site, the denominator of the '
                       'top events\' percent of samples. The authors found no sample count or depth threshold '
                       'separated true from false events. MANE transcripts only.',
        'examples': [16054],
        'see': _SPLICEVAULT_SEE,
        'value_type': 'I',
        'order': 23,
        'variantgrid_column_id': 'splicevault_site_sample_count',
    },
    {
        'key': 'utr_annotator_consequence',
        'label': "5' UTR uORF consequence",
        'sub_label': 'UTRAnnotator',
        'description': "UTRAnnotator effect of a 5' UTR variant on upstream ORFs: uAUG gained or lost, uORF stop "
                       'codon gained or lost, or uORF frameshift. <p> uAUG gains creating an out-of-frame overlapping '
                       'ORF are under selection similar to missense variants (Whiffin 2020).',
        'examples': ['5_prime_UTR_premature_start_codon_gain_variant'],
        'see': _UTR_ANNOTATOR_SEE,
        'order': 25,
        'variantgrid_column_id': 'utr_annotator_consequence',
    },
    {
        'key': 'utr_annotator_annotation',
        'label': "5' UTR uORF details",
        'sub_label': 'UTRAnnotator',
        'description': "The type of ORF each 5' UTR effect creates or disrupts and its Kozak strength. <p> Martin-Geary "
                       '2025 prioritised uAUG gains and uORF stop losses creating an out-of-frame overlapping ORF '
                       '(oORF) with strong or moderate Kozak, and uAUG losses with strong Kozak. An in-frame oORF '
                       '(N-terminal extension) is not necessarily loss of function.',
        'examples': ['OutOfFrame_oORF (Strong Kozak)'],
        'see': _UTR_ANNOTATOR_SEE,
        'order': 25,
        'variantgrid_column_id': 'utr_annotator_annotation',
    },
]


def _add_evidence_keys(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    COMPUTATIONAL_AND_PREDICTIVE_DATA = 'CP'
    COPY_SCOPE_NONE = 'N'

    EvidenceKey.objects.bulk_create([
        EvidenceKey(evidence_category=COMPUTATIONAL_AND_PREDICTIVE_DATA, copy_scope=COPY_SCOPE_NONE,
                    options=[], **key_data)
        for key_data in _NEW_KEYS
    ], ignore_conflicts=True)


def _remove_evidence_keys(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key__in=[k['key'] for k in _NEW_KEYS]).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0006_recalc_initiation_codon_p_hgvs'),
        ('snpdb', '0004_alphagenome_splicing_column'),
    ]

    operations = [
        migrations.RunPython(_add_evidence_keys, reverse_code=_remove_evidence_keys),
    ]
