from django.db import migrations

# SACGF/variantgrid#444 - what the report calls a splicing event ("AR-V7", "MET exon 14 skipping").
# A SpliceGirl call imports as a <DEL> variant with coordinates and a c.HGVS, so nothing in the record
# says what the event is named - the scientist names it, as they typed it into the legacy report


def _create_splice_label_ekey(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")

    FREE_ENTRY = 'F'
    V = 'V'  # VARIANT
    ALLELE = 'A'  # CopyScope.ALLELE - the same splice event is called the same thing every time
    SOMATIC = 'S'  # CopyAlleleOrigin.SOMATIC

    splice_label = EvidenceKey(
        key='splice_label',
        label='Splice label',
        description='What the report calls this splicing event. Printed in the Splicing Variants '
                    'table and sent as the JSON description, in place of a protein change.',
        examples=["AR-V7", "MET exon 14 skipping", "EGFRvIII"],
        options=[],
        evidence_category=V,
        value_type=FREE_ENTRY,
        order=15,
        mandatory=False,
        max_share_level='public',
        copy_scope=ALLELE,
        copy_allele_origin=SOMATIC,
        variantgrid_column_id=None,
    )
    EvidenceKey.objects.bulk_create([splice_label], ignore_conflicts=True)


def _delete_splice_label_ekey(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(pk='splice_label').delete()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0181_remove_json_template'),
    ]

    operations = [
        migrations.RunPython(_create_splice_label_ekey, reverse_code=_delete_splice_label_ekey),
    ]
