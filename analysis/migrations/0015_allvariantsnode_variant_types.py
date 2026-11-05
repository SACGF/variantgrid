import django.contrib.postgres.fields
from django.db import migrations, models

import snpdb.variant_filters
from snpdb.variant_filters import VariantType, get_all_variant_types, get_symbolic_variant_types

# Each checkbox the node had and the type tokens it stood for - "structural variants" was every symbolic alt
TYPES_FOR_FIELD = {
    "reference": lambda: [VariantType.REFERENCE],
    "snps": lambda: [VariantType.SNV],
    "indels": lambda: [VariantType.INDEL],
    "complex_subsitution": lambda: [VariantType.COMPLEX],
    "structural_variants": get_symbolic_variant_types,
    "fusions": lambda: [VariantType.FUSION],
    "copy_number": lambda: [VariantType.COPY_NUMBER],
    "splicing": lambda: [VariantType.SPLICE],
}


def booleans_to_variant_types(ticked: dict[str, bool], all_types: list[str]) -> list[str]:
    """ The node's old checkboxes to the All Variants grid's types (#2680), in that list's order """
    selected = {vt for field, types_for in TYPES_FOR_FIELD.items() if ticked[field] for vt in types_for()}
    return [vt for vt in all_types if vt in selected]


def forwards(apps, _schema_editor):
    AllVariantsNode = apps.get_model("analysis", "AllVariantsNode")
    all_types = get_all_variant_types(include_reference=True)
    for node in AllVariantsNode.objects.all():
        ticked = {field: getattr(node, field) for field in TYPES_FOR_FIELD}
        node.variant_types = booleans_to_variant_types(ticked, all_types)
        node.save(update_fields=["variant_types"])


class Migration(migrations.Migration):

    dependencies = [
        ("analysis", "0014_damagenode_variant_types"),
    ]

    operations = [
        migrations.AddField(
            model_name="allvariantsnode",
            name="variant_types",
            field=django.contrib.postgres.fields.ArrayField(
                base_field=models.CharField(max_length=20),
                blank=True,
                default=snpdb.variant_filters.get_all_variant_types,
                size=None,
            ),
        ),
        migrations.RunPython(forwards, migrations.RunPython.noop),
        migrations.RemoveField(model_name="allvariantsnode", name="complex_subsitution"),
        migrations.RemoveField(model_name="allvariantsnode", name="copy_number"),
        migrations.RemoveField(model_name="allvariantsnode", name="fusions"),
        migrations.RemoveField(model_name="allvariantsnode", name="indels"),
        migrations.RemoveField(model_name="allvariantsnode", name="reference"),
        migrations.RemoveField(model_name="allvariantsnode", name="snps"),
        migrations.RemoveField(model_name="allvariantsnode", name="splicing"),
        migrations.RemoveField(model_name="allvariantsnode", name="structural_variants"),
    ]
