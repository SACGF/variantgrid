from django.db import migrations

SPLICING_VARIANT_SO = "SO:0001568"


def _add_splicing_variant_options(apps, _schema_editor):
    """ The gene-level annotation pipeline writes splicing_variant as a splice call's consequence and
        variant class, and autopopulate copies both straight into a classification - so the options
        list needs to carry it. @see annotation.gene_level_annotation """

    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    for key, label in [("molecular_consequence", "Splicing variant"), ("variant_class", "Splicing variant")]:
        ekey = EvidenceKey.objects.filter(pk=key).first()
        if ekey is None:
            continue
        options = ekey.options or []
        if any(o.get("key") == "splicing_variant" for o in options):
            continue
        index = max((o.get("index") or 0) for o in options) + 1 if options else 1
        options.append({"so": SPLICING_VARIANT_SO, "key": "splicing_variant", "index": index, "label": label})
        ekey.options = options
        ekey.save()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0182_splice_label_ekey'),
    ]

    operations = [
        migrations.RunPython(_add_splicing_variant_options, migrations.RunPython.noop),
    ]
