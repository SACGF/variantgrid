# Generated by Django 4.2.10 on 2024-04-10 04:17
import operator
from functools import reduce

from django.db import migrations
from django.db.models import Q

from classification.models import ImportedAlleleInfoStatus
from manual.operations.manual_operations import ManualOperation


def _test_has_iai_with_provided_ref(apps):
    ImportedAlleleInfo = apps.get_model("classification", "ImportedAlleleInfo")
    regex_span_plus_provided = r"\d+_\d+(del|dup)[GATC]+"
    filters = [
        Q(imported_c_hgvs__regex=regex_span_plus_provided),
        Q(imported_g_hgvs__regex=regex_span_plus_provided)
    ]
    q = ~Q(status=ImportedAlleleInfoStatus.FAILED) & reduce(operator.or_, filters)
    return ImportedAlleleInfo.objects.filter(q).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0134_new_ekeys_tags_and_mave'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_legacy_imported_allele_info_hgvs_error_provided_reference_length"]), test=_test_has_iai_with_provided_ref)
    ]
