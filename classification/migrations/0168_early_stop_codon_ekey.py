from django.db import migrations

# sacgf/variantgrid#579 - predicted stop codon position for frameshift variants


def _add_early_stop_codon_ekey(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")

    INTEGER = 'I'
    CP = 'CP'  # COMPUTATIONAL_AND_PREDICTIVE_DATA

    early_stop_codon = EvidenceKey(
        key='early_stop_codon',
        label='Early stop codon',
        sub_label='Loss of function',
        description='Codons from the first changed residue to the frameshift\'s new premature '
                    'termination codon, read from the HGVS p. "fsTer&lt;n&gt;". 1 means the changed '
                    'codon is itself the stop. Drives the last-exon / 50nt NMD rules for PVS1.',
        examples=[1, 49],
        options=[],
        see='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6185798/',
        evidence_category=CP,
        value_type=INTEGER,
        order=10,
        mandatory=False,
        max_share_level='logged_in_users',
        copy_consensus=True,
        variantgrid_column_id='ptc_distance_codons',
    )
    # ignore_conflicts: a restored DB may already carry this key from a later snapshot
    EvidenceKey.objects.bulk_create([early_stop_codon], ignore_conflicts=True)


def _reverse_early_stop_codon_ekey(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(pk='early_stop_codon').delete()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0167_one_off_rematch_alignment_gap_allele_infos'),
        ('snpdb', '0200_ptc_columns'),
    ]

    operations = [
        migrations.RunPython(_add_early_stop_codon_ekey, reverse_code=_reverse_early_stop_codon_ekey),
    ]
