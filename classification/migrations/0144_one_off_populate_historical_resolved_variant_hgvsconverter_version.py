# Generated by Django 4.2.9 on 2024-06-26 04:18
from importlib import metadata

from django.conf import settings
from django.db import migrations

from genes.hgvs import HGVSConverterType


def _one_off_populate_historical_resolved_variant_hgvsconverter_version(apps, _schema_editor):
    HGVSConverterVersion = apps.get_model('classification', 'HGVSConverterVersion')
    ResolvedVariantInfo = apps.get_model('classification', 'ResolvedVariantInfo')

    rvi_qs = ResolvedVariantInfo.objects.all()
    if rvi_qs.exists():
        hgvs_converter_type = HGVSConverterType[settings.HGVS_DEFAULT_METHOD.upper()]
        if hgvs_converter_type == HGVSConverterType.BIOCOMMONS_HGVS:
            version = metadata.version('hgvs')
        else:
            version = metadata.version('pyhgvs')

        hcv = HGVSConverterVersion.objects.get_or_create(hgvs_converter_type=hgvs_converter_type.name,
                                                         version=f"Legacy set from settings: '{version}' or earlier",
                                                         method='unknown tool or clingen fallback',
                                                         code_git_hash='not-a-real-git-hash')[0]
        rvi_qs.update(c_hgvs_converter_version=hcv)


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0143_hgvsconverterversion_and_more'),
    ]

    operations = [
        migrations.RunPython(_one_off_populate_historical_resolved_variant_hgvsconverter_version)
    ]
