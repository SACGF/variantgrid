from django.db import migrations

# SACGF/variantgrid#444 - the caller's copy ratio, beside copy_number.
# CN is an absolute count, SM/FC are ratios against the normal, so they are different quantities and
# get different keys - @see library/genomics/vcf_enums.py VCFConstant.COPY_NUMBER_FIELD_IS_RATIO


def _create_fold_change_ekey(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")

    FLOAT = 'L'
    V = 'V'  # VARIANT
    NONE = 'N'  # CopyScope.NONE - it is the caller's measurement, so it never copies

    fold_change = EvidenceKey(
        key='fold_change',
        label='Fold change',
        description='Copy ratio of the call against the normal (matches VCF SM or FC). '
                    'A whole-gene amplification prints this alongside the copy number.',
        examples=[0.53, 1.9, 4.31],
        options=[],
        evidence_category=V,
        value_type=FLOAT,
        order=14,
        mandatory=False,
        max_share_level='public',
        copy_scope=NONE,
        variantgrid_column_id=None,
    )
    EvidenceKey.objects.bulk_create([fold_change], ignore_conflicts=True)


def _delete_fold_change_ekey(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(pk='fold_change').delete()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0179_default_case_report_template'),
    ]

    operations = [
        migrations.RunPython(_create_fold_change_ekey, reverse_code=_delete_fold_change_ekey),
    ]
