import os
from pathlib import Path
from typing import List
from django.conf import settings

import celery

from library.utils import execute_cmd
from snpdb.models import VCF, Cohort, Trio, GenomeBuild, SomalierVCFExtract


class SomalierConfig:
    def __init__(self):
        self.settings = settings.SOMALIER

    def _annotation_dir(self, dirname):
        return os.path.join(self.settings["annotation_base_dir"], dirname)

    def get_annotation(self, key):
        return self._annotation_dir(self.settings["annotation"][key])

    def report_dir(self, uuid):
        return os.path.join(self.settings["report_base_dir"], uuid)

    def get_sites(self, genome_build: GenomeBuild):
        sites = self.settings["annotation"]["sites"][genome_build.name]
        return self._annotation_dir(sites)

    def __getitem__(self, key):
        return self.settings[key]


@celery.task
def somalier_vcf(vcf_id: int):
    """ Extract, Ancestry, and Relate (Cohort) """
    vcf = VCF.objects.get(pk=vcf_id)
    sve = SomalierVCFExtract.objects.create(vcf=vcf)

    uf = vcf.uploadedvcf.uploaded_file
    vcf_filename = uf.get_filename()
    if not os.path.exists(vcf_filename):
        raise FileNotFoundError(vcf_filename)

    # The VCF needs to be tabix indexed
    # Check if can tabix. If it fails then create a temp file
    # tabix that (remember to delete temp dir)
    tabix_command = ["tabix", vcf_filename]
    execute_cmd(tabix_command)

    # Extract
    cfg = SomalierConfig()
    somalier_bin = cfg.get_annotation("command")
    command = [somalier_bin, "extract",
               "--out-dir", sve.get_dir(),
               "--sites", cfg.get_sites(vcf.genome_build),
               "--fasta", vcf.genome_build.reference_fasta,
               vcf_filename]
    execute_cmd(command)

    # Ancestry

    # Relate


def _somalier_relate(uuid: str, ped_filename: str, samples: List[str]):
    cfg = SomalierConfig()
    somalier_bin = cfg.get_annotation("command")
    with Path(cfg.report_dir(uuid)):
        command = [somalier_bin, "relate", "--ped", ped_filename] + samples
        execute_cmd(command)


@celery.task
def somalier_cohort_relate(cohort_id: int):
    cohort = Cohort.objects.get(pk=cohort_id)


@celery.task
def somalier_trio_relate(trio_id: int):
    trio = Trio.objects.get(pk=trio_id)


# TODO: Do for Pedigree