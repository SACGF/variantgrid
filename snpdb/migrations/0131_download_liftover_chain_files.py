# Generated by Django 4.2.10 on 2024-05-08 01:08
import os

from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


LIFTOVER_37_TO_38 = settings.ANNOTATION["GRCh37"]["liftover"]["GRCh38"]
LIFTOVER_38_TO_37 = settings.ANNOTATION["GRCh38"]["liftover"]["GRCh37"]


def _test_has_missing_37_to_38_chain_file(apps):
    return not os.path.exists(LIFTOVER_37_TO_38)


def _test_has_missing_38_to_37_chain_file(apps):
    return not os.path.exists(LIFTOVER_38_TO_37)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0130_one_off_legacy_populate_allele_liftover'),
    ]

    operations = [
        ManualOperation.operation_other(args=[
            f"Download 37 -> 38 chain file go to {os.path.dirname(LIFTOVER_37_TO_38)} and run: 'wget http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz'"],
            test=_test_has_missing_37_to_38_chain_file),
        ManualOperation.operation_other(args=[
            f"Download 38 -> 37 chain file go to {os.path.dirname(LIFTOVER_38_TO_37)} and run: 'wget http://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh38_to_GRCh37.chain.gz'"],
            test=_test_has_missing_37_to_38_chain_file),
    ]
