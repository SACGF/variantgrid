from django.db import migrations

# SACGF/variantgrid#2583 - the key's help text showed SpliceVault's data file format (fields joined by ';'), but
# VEP's VCF output turns ';' into ':', which is what the column holds and the key is autopopulated with.
_KEY = 'splicevault_top_events'
_RANK_TYPE_POSITION = "each rank:type;position;percent of samples;frame."
_FIXED_RANK_TYPE_POSITION = "each rank:type:position:percent of samples:frame."
_EXAMPLES = ['Top1:CD;-35;31%;Frameshift, Top2:ES;1;5%;Inframe']
_FIXED_EXAMPLES = ['Top1:CD:-35:31%:Frameshift, Top2:ES:1:5%:inFrame']


def _set_help_text(apps, old: str, new: str, examples: list[str]):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    if ekey := EvidenceKey.objects.filter(key=_KEY).first():
        ekey.description = ekey.description.replace(old, new)
        ekey.examples = examples
        ekey.save(update_fields=["description", "examples"])


def _fix_separator(apps, _schema_editor):
    _set_help_text(apps, _RANK_TYPE_POSITION, _FIXED_RANK_TYPE_POSITION, _FIXED_EXAMPLES)


def _restore_separator(apps, _schema_editor):
    _set_help_text(apps, _FIXED_RANK_TYPE_POSITION, _RANK_TYPE_POSITION, _EXAMPLES)


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0007_avi_splicevault_utr_annotator_ekeys'),
    ]

    operations = [
        migrations.RunPython(_fix_separator, reverse_code=_restore_separator),
    ]
