from django.db import migrations, models

# #2596 - the PTC escapes NMD checkbox becomes a choice: a ticked box is PTCNMDFilter.ESCAPING
PTC_NMD_ESCAPING = 'E'


def _ptc_nmd_escaping_to_ptc_nmd(apps, _schema_editor):
    DamageNode = apps.get_model("analysis", "DamageNode")
    DamageNode.objects.filter(ptc_nmd_escaping=True).update(ptc_nmd=PTC_NMD_ESCAPING)


def _ptc_nmd_to_ptc_nmd_escaping(apps, _schema_editor):
    DamageNode = apps.get_model("analysis", "DamageNode")
    DamageNode.objects.filter(ptc_nmd__isnull=False).update(ptc_nmd_escaping=True)


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0008_analysisnode_node_type_not_null'),
    ]

    operations = [
        migrations.AddField(
            model_name='damagenode',
            name='ptc_nmd',
            field=models.CharField(blank=True, choices=[('E', 'Escapes NMD'), ('A', 'Escapes NMD or NMD predicted')], max_length=1, null=True),
        ),
        migrations.RunPython(_ptc_nmd_escaping_to_ptc_nmd, reverse_code=_ptc_nmd_to_ptc_nmd_escaping),
        migrations.RemoveField(
            model_name='damagenode',
            name='ptc_nmd_escaping',
        ),
        migrations.RenameField(
            model_name='damagenode',
            old_name='ptc_nmd_escaping_required',
            new_name='ptc_nmd_required',
        ),
    ]
