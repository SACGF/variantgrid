from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_dispatchable_variant_annotation_versions(apps):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(status__in=["NEW", "ACTIVE"]).exists()


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0016_annotsv_sort_regulatory_elements"),
        ("manual", "0002_seed_data"),
    ]

    operations = [
        ManualOperation.operation_manage(["one_off_reannotate_converted_symbolic_variants"],
                                         test=_test_has_dispatchable_variant_annotation_versions),
    ]
