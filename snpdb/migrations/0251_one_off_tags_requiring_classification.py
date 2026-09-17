from django.conf import settings
from django.db import migrations

from classification.enums import AlleleOriginBucket

SOMATIC_REPORTABLE = "SomaticReportable"


def _set_tags_requiring_classification(apps, _schema_editor):
    """ The tags that were already being used as 'this needs classifying' - RequiresClassification via
        settings, and SomaticReportable which somatic labs keep on the variant after classifying """
    Tag = apps.get_model("snpdb", "Tag")
    Tag.objects.filter(pk__in=[settings.TAG_REQUIRES_CLASSIFICATION, SOMATIC_REPORTABLE]) \
        .update(requires_classification=True)
    Tag.objects.filter(pk=SOMATIC_REPORTABLE, allele_origin_bucket=AlleleOriginBucket.UNKNOWN) \
        .update(allele_origin_bucket=AlleleOriginBucket.SOMATIC)


def _clear_tags_requiring_classification(apps, _schema_editor):
    Tag = apps.get_model("snpdb", "Tag")
    Tag.objects.update(requires_classification=False)


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0250_tag_requires_classification"),
    ]

    operations = [
        migrations.RunPython(_set_tags_requiring_classification, _clear_tags_requiring_classification),
    ]
