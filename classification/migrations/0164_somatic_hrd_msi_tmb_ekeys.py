from django.db import migrations

# sacgf/variantgrid_private#3840


def _add_somatic_hrd_msi_tmb_ekeys(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")

    FLOAT = 'L'
    SELECT = 'S'
    TEXT_AREA = 'T'

    TEST = 'HT'          # EvidenceCategory.HEADER_TEST
    INTERPRETATION = 'HI'  # EvidenceCategory.INTERPRETATION

    new_ekeys = [
        EvidenceKey(
            key='somatic:tmb_value',
            label='TMB',
            description='Tumour Mutational Burden in mutations per megabase (mut/Mb)',
            value_type=FLOAT,
            evidence_category=TEST,
            max_share_level='logged_in_users',
            mandatory=False,
            options=[],
        ),
        EvidenceKey(
            key='somatic:tmb_status',
            label='TMB Status',
            description='Tumour Mutational Burden classification category',
            value_type=SELECT,
            evidence_category=TEST,
            max_share_level='logged_in_users',
            mandatory=False,
            options=[
                {'key': 'low', 'index': 1, 'label': 'Low'},
                {'key': 'high', 'index': 2, 'label': 'High'},
                {'key': 'unable_to_be_determined', 'index': 3, 'label': 'Unable to be determined'},
            ],
        ),
        EvidenceKey(
            key='somatic:msi_value',
            label='MSI %',
            description='Microsatellite Instability percentage of unstable sites',
            value_type=FLOAT,
            evidence_category=TEST,
            max_share_level='logged_in_users',
            mandatory=False,
            options=[],
        ),
        EvidenceKey(
            key='somatic:msi_status',
            label='MSI Status',
            description='Microsatellite Instability status classification',
            value_type=SELECT,
            evidence_category=TEST,
            max_share_level='logged_in_users',
            mandatory=False,
            options=[
                {'key': 'msi_high', 'index': 1, 'label': 'MSI-H (High)'},
                {'key': 'msi_low', 'index': 2, 'label': 'MSI-L (Low)'},
                {'key': 'ms_stable', 'index': 3, 'label': 'MS-S (Stable)'},
            ],
        ),
        EvidenceKey(
            key='somatic:hrd_status',
            label='HRD Status',
            description='Homologous Recombination Deficiency status',
            value_type=SELECT,
            evidence_category=TEST,
            max_share_level='logged_in_users',
            mandatory=False,
            options=[
                {'key': 'positive', 'index': 1, 'label': 'Positive'},
                {'key': 'negative', 'index': 2, 'label': 'Negative'},
                {'key': 'unable_to_be_determined', 'index': 3, 'label': 'Unable to be determined'},
            ],
        ),
        # summary_interpretation very likely to contain patient-identifiable information.
        # We will insert TSO500 data into here in VariantGrid but not send to Shariant
        # May change after review by molecular oncology
        EvidenceKey(
            key='somatic:summary_interpretation',
            label='Summary Interpretation',
            description='Free-text summary of all classifications in a report',
            value_type=TEXT_AREA,
            evidence_category=INTERPRETATION,
            max_share_level='lab',
            mandatory=False,
            options=[],
        ),
    ]
    EvidenceKey.objects.bulk_create(new_ekeys)


def _remove_somatic_hrd_msi_tmb_ekeys(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key__in=[
        'somatic:tmb_value',
        'somatic:tmb_status',
        'somatic:msi_value',
        'somatic:msi_status',
        'somatic:hrd_status',
        'somatic:summary_interpretation',
    ]).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0163_alter_classification_share_level_and_more'),
    ]

    operations = [
        migrations.RunPython(_add_somatic_hrd_msi_tmb_ekeys, reverse_code=_remove_somatic_hrd_msi_tmb_ekeys)
    ]
