import glob
import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import List
from django.conf import settings

import celery

from library.django_utils.django_file_utils import get_import_processing_dir
from library.utils import execute_cmd
from snpdb.models import VCF, Cohort, Trio, SomalierVCFExtract, SomalierConfig, SomalierAncestryRun, ProcessingStatus, \
    SomalierCohortRelate, SomalierRelate
from snpdb.variants_to_vcf import vcf_export_to_file


def _write_somalier_vcf(cfg: SomalierConfig, processing_dir, vcf: VCF):
    sites_name = os.path.basename(cfg.get_sites(vcf.genome_build))
    sites_vcf_kwargs = {"name": sites_name, "genome_build": vcf.genome_build}
    try:
        sites_vcf = VCF.objects.get(**sites_vcf_kwargs)
    except VCF.DoesNotExist as dne:
        print(f"Expected single VCF loaded via: {sites_vcf_kwargs}")
        raise dne

    sites_qs = sites_vcf.get_variant_qs()
    exported_vcf_filename = os.path.join(processing_dir, f"vcf_{vcf.pk}.vcf.bgz")
    vcf_export_to_file(vcf, exported_vcf_filename, original_qs=sites_qs)

    tabix_command = ["tabix", exported_vcf_filename]
    return_code, stdout, stderr = execute_cmd(tabix_command)
    print(f"return_code: {return_code}, stdout: {stdout}, stderr: {stderr}")
    return exported_vcf_filename

@celery.task
def somalier_vcf(vcf_id: int):
    """ Extract, Ancestry, and Relate (Cohort) """
    vcf = VCF.objects.get(pk=vcf_id)
    vcf_extract = SomalierVCFExtract.objects.create(vcf=vcf)
    processing_dir = get_import_processing_dir(vcf_extract.pk, "somalier_vcf_extract")
    cfg = SomalierConfig()
    vcf_filename = _write_somalier_vcf(cfg, processing_dir, vcf)

    # Extract
    somalier_bin = cfg.get_annotation("command")
    extract_cmd = [somalier_bin, "extract",
                   "--out-dir", vcf_extract.get_somalier_dir(),
                   "--sites", cfg.get_sites(vcf.genome_build),
                   "--fasta", vcf.genome_build.reference_fasta,
                   vcf_filename]
    vcf_extract.execute(extract_cmd)

    # Ancestry
    ancestry_run = SomalierAncestryRun.objects.create(vcf_extract=vcf_extract)
    compare_samples = glob.glob1(cfg.get_annotation("ancestry_somalier_dir"), "*.somalier")
    ancestry_cmd = [somalier_bin, "ancestry",
                    "--out-dir", ancestry_run.get_report_dir(),
                    "--labels", cfg.get_annotation("ancestry_labels")] \
                   + compare_samples + ["++"] + ancestry_run.get_sample_somalier_filenames()
    ancestry_run.execute(ancestry_cmd)

    # TODO: Load SomalierAncestry for each

    # Relate
    relate = SomalierCohortRelate.objects.create(cohort=vcf.cohort, cohort_version=vcf.cohort.version)
    _somalier_relate(relate)

    if settings.VCF_IMPORT_DELETE_TEMP_FILES_ON_SUCCESS:
        shutil.rmtree(processing_dir)


def _somalier_relate(somalier_relate: SomalierRelate):
    cfg = SomalierConfig()
    somalier_bin = cfg.get_annotation("command")
    processing_dir = get_import_processing_dir(somalier_relate.pk, "somalier_relate")

    ped_filename = os.path.join(processing_dir, "temp.ped")
    somalier_relate.write_ped_file(ped_filename)

    with Path(cfg.related_dir(somalier_relate.uuid)):
        command = [somalier_bin, "relate", "--ped", ped_filename] \
                  + somalier_relate.get_sample_somalier_filenames()
        somalier_relate.execute(command)

    shutil.rmtree(processing_dir)


@celery.task
def somalier_cohort_relate(cohort_id: int):
    cohort = Cohort.objects.get(pk=cohort_id)


@celery.task
def somalier_trio_relate(trio_id: int):
    trio = Trio.objects.get(pk=trio_id)


# TODO: Do for Pedigree