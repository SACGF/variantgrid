from django.db import migrations, models

# sacgf/variantgrid#1714 - copy_consensus (bool) becomes copy_scope, so a value can say how far it
# travels, plus copy_allele_origin so germline concepts stay out of somatic records


def _copy_consensus_to_copy_scope(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(copy_consensus=True).update(copy_scope="A")
    EvidenceKey.objects.filter(copy_consensus=False).update(copy_scope="N")


def _copy_scope_to_copy_consensus(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    EvidenceKey.objects.filter(copy_scope="N").update(copy_consensus=False)
    EvidenceKey.objects.exclude(copy_scope="N").update(copy_consensus=True)


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0175_ekey_gene_fusion_options'),
    ]

    operations = [
        migrations.AddField(
            model_name='evidencekey',
            name='copy_scope',
            field=models.CharField(choices=[('N', 'None'), ('A', 'Allele'), ('G', 'Gene')], default='A',
                                   max_length=1),
        ),
        migrations.AddField(
            model_name='evidencekey',
            name='copy_allele_origin',
            field=models.CharField(choices=[('A', 'Any'), ('G', 'Germline'), ('S', 'Somatic')], default='A',
                                   max_length=1),
        ),
        migrations.RunPython(_copy_consensus_to_copy_scope, reverse_code=_copy_scope_to_copy_consensus),
        migrations.RemoveField(
            model_name='evidencekey',
            name='copy_consensus',
        ),
    ]
