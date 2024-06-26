# Generated by Django 4.2.10 on 2024-03-25 00:09
import os

from django.conf import settings
from django.db import migrations

from library.utils import file_sha256sum
from manual.operations.manual_operations import ManualOperation


def _structuralvariantoverlap_hash_check(apps):
    """ Return True if file hash is different """

    # We renamed the setting in 0097 switching from plugin -> custom VCF
    SV_SETTINGS = ["structuralvariantoverlap", "gnomad_sv"]

    for sv_setting_name in SV_SETTINGS:
        try:
            sv_basename = settings.ANNOTATION['GRCh37']['vep_config'][sv_setting_name]
            sv_filename = os.path.join(settings.ANNOTATION_VEP_BASE_DIR, sv_basename)
            # If this file is not found, it should be caught in "deployment_check" and the new one will be downloaded
            if os.path.exists(sv_filename):
                print(f"Checking hash of '{sv_filename}'")
                sv_hash = file_sha256sum(sv_filename)
                if sv_hash != "01e72bcf6fa9efb0346f4baf437fdd82bceaee7b6aad3921da78da4e581df6c5":
                    return True
        except KeyError:
            pass

        try:
            sv_basename = settings.ANNOTATION['GRCh38']['vep_config'][sv_setting_name]
            sv_filename = os.path.join(settings.ANNOTATION_VEP_BASE_DIR, sv_basename)
            # If this file is not found, it should be caught in "deployment_check" and the new one will be downloaded
            if os.path.exists(sv_filename):
                print(f"Checking hash of '{sv_filename}'")
                sv_hash = file_sha256sum(sv_filename)
                if sv_hash != "a5ffd9a0c8369e2fdaeec1ae23989021d851f1e90dec9cabc0ac26a29f665775":
                    return True
        except KeyError:
            pass

    return False


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0094_gnomad_sv_vep_fields'),
    ]

    operations = [
        ManualOperation.operation_other(args=[
            "Update copies of 'gnomad_v2.1_sv.sites.grch37.converted.vcf.gz' and 'gnomad_v2.1_sv.sites.grch37.converted.vcf.gz' (and .tbis)"],
                                        test=_structuralvariantoverlap_hash_check),

    ]
