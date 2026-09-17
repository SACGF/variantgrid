from django.db import migrations

GENE_FUSION_SO = "SO:0001565"


def _add_gene_fusion_options(apps, _schema_editor):
    """ The gene-level annotation pipeline writes gene_fusion as a fusion's consequence and variant
        class, and autopopulate copies both straight into a classification - so the options list needs
        to carry it. @see annotation.gene_level_annotation """

    EvidenceKey = apps.get_model("classification", "EvidenceKey")
    for key, label in [("molecular_consequence", "Gene fusion"), ("variant_class", "Gene fusion")]:
        ekey = EvidenceKey.objects.filter(pk=key).first()
        if ekey is None:
            continue
        options = ekey.options or []
        if any(o.get("key") == "gene_fusion" for o in options):
            continue
        index = max((o.get("index") or 0) for o in options) + 1 if options else 1
        options.append({"so": GENE_FUSION_SO, "key": "gene_fusion", "index": index, "label": label})
        ekey.options = options
        ekey.save()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0174_reclassificationevent_event_type_and_more'),
    ]

    operations = [
        migrations.RunPython(_add_gene_fusion_options, migrations.RunPython.noop),
    ]
