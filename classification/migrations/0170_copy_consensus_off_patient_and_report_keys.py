from django.db import migrations

# sacgf/variantgrid#1713 - copy consensus was copying facts about the patient's own tumour and
# free text written about a different record

_TEST_LEVEL_KEYS = [
    'somatic:tmb_value',
    'somatic:tmb_status',
    'somatic:msi_value',
    'somatic:msi_status',
    'somatic:hrd_status',
    'testing_context',
]
"""
Measurements of the specimen being reported rather than properties of the allele - copying them imports
another patient's numbers. The five somatic: ones reach a classification from SpecimenMeasure (#1559).
"""

_FREE_TEXT_KEYS = [
    'somatic:summary_interpretation',
    'review_comment',
]
"""
somatic:summary_interpretation summarises every classification in a report and is max_share_level='lab'
because it is likely to carry patient-identifiable information; review_comment is a reviewer's comment
about a different record.
"""

_KEYS = _TEST_LEVEL_KEYS + _FREE_TEXT_KEYS


def _copy_consensus_off(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key__in=_KEYS).update(copy_consensus=False)


def _copy_consensus_on(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(key__in=_KEYS).update(copy_consensus=True)


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0169_one_off_revalidate_g_hgvs_only_allele_infos'),
    ]

    operations = [
        migrations.RunPython(_copy_consensus_off, reverse_code=_copy_consensus_on),
    ]
