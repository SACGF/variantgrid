"""
Tags created while a variant had duplicate VariantAlleles never got an allele - _liftover_variant_tag died on
MultipleObjectsReturned before it could assign one, so the tag was invisible in the other build (#1361).

snpdb.0259 deduped the VariantAlleles, so now there's exactly one allele to give each of them.
"""
from django.db import migrations


def _backfill_variant_tag_allele(apps, schema_editor):
    VariantTag = apps.get_model("analysis", "VariantTag")
    VariantAllele = apps.get_model("snpdb", "VariantAllele")

    for tag in VariantTag.objects.filter(allele__isnull=True).iterator(chunk_size=2000):
        # A tag whose variant has no VariantAllele at all is left for variant_tag_created_task to pick up
        if va := VariantAllele.objects.filter(variant_id=tag.variant_id,
                                              genome_build_id=tag.genome_build_id).first():
            tag.allele_id = va.allele_id
            tag.save(update_fields=["allele_id"])


class Migration(migrations.Migration):
    dependencies = [
        ("analysis", "0141_varianttag_varianttag_one_per_sample_in_analysis"),
        ("snpdb", "0259_variantallele_unique_variant_build"),
    ]

    operations = [
        migrations.RunPython(_backfill_variant_tag_allele, migrations.RunPython.noop),
    ]
