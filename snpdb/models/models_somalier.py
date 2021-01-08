import logging
import os
import uuid
from subprocess import CalledProcessError
from typing import List, Iterable

from django.conf import settings
from django.db import models
from django.db.models import CASCADE
from django.utils.text import slugify
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager

from library.utils import execute_cmd
from pedigree.ped.export_ped import write_unrelated_ped, write_trio_ped
from snpdb.models import Sample, VCF, Cohort, Trio, SuperPopulationCode, ImportStatus
from snpdb.models.models_enums import ProcessingStatus
from patients.models_enums import Sex


class AbstractSomalierModel(TimeStampedModel):
    status = models.CharField(max_length=1, choices=ProcessingStatus.choices, default=ProcessingStatus.CREATED)
    error_exception = models.TextField(null=True, blank=True)

    class Meta:
        abstract = True

    def get_samples(self) -> Iterable[Sample]:
        raise NotImplementedError()

    def get_sample_somalier_filenames(self) -> List[str]:
        return [AbstractSomalierModel.sample_filename(s) for s in self.get_samples()]

    def execute(self, command: List[str], **kwargs):
        """ Executes code and handles saving errors """

        cmd = " ".join(command)
        logging.info('About to call %s', cmd)

        self.status = ProcessingStatus.PROCESSING
        self.save()

        return_code, stdout, stderr = execute_cmd(command, **kwargs)
        if return_code != 0:
            self.error_exception = f"return_code: {return_code}. stdout: {stdout}, stderr: {stderr}"
            self.status = ProcessingStatus.ERROR
        else:
            self.status = ProcessingStatus.SUCCESS
        self.save()

        if return_code != 0:
            raise CalledProcessError(returncode=return_code, cmd=cmd, output=self.error_exception)

    @staticmethod
    def sample_name(sample: Sample):
        # Add PK as suffix so they're all unique
        return f"{slugify(sample.vcf_sample_name)}_{sample.pk}"

    @staticmethod
    def sample_filename(sample: Sample):
        vcf_dir = sample.vcf.somaliervcfextract.get_somalier_dir()
        return os.path.join(vcf_dir, AbstractSomalierModel.sample_name(sample) + ".somalier")

    @staticmethod
    def sample_name_to_id(sample_name: str):
        """ Sample ID is stored at the end """
        return sample_name.rsplit("_", 1)[-1]

    @staticmethod
    def media_url(file_path):
        # Need to use a slash, so that later joins don't have absolute path
        media_root_with_slash = os.path.join(settings.MEDIA_ROOT, "")
        if not file_path.startswith(media_root_with_slash):
            raise ValueError(f"'{file_path}' must start with MEDIA_ROOT: {media_root_with_slash}")

        return os.path.join(settings.MEDIA_URL, file_path[len(media_root_with_slash):])


class SomalierVCFExtract(AbstractSomalierModel):
    vcf = models.OneToOneField(VCF, on_delete=CASCADE)

    def get_somalier_dir(self):
        cfg = SomalierConfig()
        return os.path.join(cfg["vcf_base_dir"], str(self.vcf.pk))

    def get_samples(self) -> Iterable[Sample]:
        return self.vcf.sample_set.filter(no_dna_control=False).order_by("pk")


class SomalierSampleExtract(models.Model):
    vcf_extract = models.ForeignKey(SomalierVCFExtract, on_delete=CASCADE)
    sample = models.OneToOneField(Sample, on_delete=CASCADE)
    ref_count = models.IntegerField(default=0)
    het_count = models.IntegerField(default=0)
    hom_count = models.IntegerField(default=0)
    unk_count = models.IntegerField(default=0)


class SomalierAncestryRun(AbstractSomalierModel):
    """ We do a run against a whole VCF """
    vcf_extract = models.OneToOneField(SomalierVCFExtract, on_delete=CASCADE)
    uuid = models.UUIDField(default=uuid.uuid4, editable=False)  # code to hide directories in media_root

    def get_report_dir(self):
        cfg = SomalierConfig()
        return cfg.ancestry_dir(self.uuid)

    def get_samples(self) -> Iterable[Sample]:
        return self.vcf_extract.get_samples()

    @property
    def url(self):
        report_dir = self.get_report_dir()
        return self.media_url(os.path.join(report_dir, "somalier-ancestry.somalier-ancestry.html"))


class SomalierAncestry(TimeStampedModel):
    ancestry_run = models.ForeignKey(SomalierAncestryRun, on_delete=CASCADE)
    sample_extract = models.OneToOneField(SomalierSampleExtract, on_delete=CASCADE)
    predicted_ancestry = models.CharField(max_length=1, choices=SuperPopulationCode.choices)
    EAS_prob = models.FloatField()
    AFR_prob = models.FloatField()
    AMR_prob = models.FloatField()
    SAS_prob = models.FloatField()
    EUR_prob = models.FloatField()


class SomalierRelate(AbstractSomalierModel):
    objects = InheritanceManager()
    uuid = models.UUIDField(default=uuid.uuid4, editable=False)  # code to hide directories in media_root

    class Meta:
        abstract = True

    def write_ped_file(self, filename):
        """ Sample IDs have to match samples provided in get_samples() """
        write_unrelated_ped(filename, [AbstractSomalierModel.sample_name(s) for s in self.get_samples()])

    @property
    def url(self):
        cfg = SomalierConfig()
        return self.media_url(os.path.join(cfg.related_dir(self.uuid), "somalier.html"))


class SomalierCohortRelate(SomalierRelate):
    cohort = models.OneToOneField(Cohort, on_delete=CASCADE)
    cohort_version = models.IntegerField()

    def get_samples(self) -> Iterable[Sample]:
        return self.cohort.get_samples().filter(no_dna_control=False)


class SomalierTrioRelate(SomalierRelate):
    trio = models.OneToOneField(Trio, on_delete=CASCADE)

    def get_samples(self) -> Iterable[Sample]:
        return self.trio.get_samples()

    def write_ped_file(self, filename):
        proband = AbstractSomalierModel.sample_name(self.trio.proband.sample)
        father = AbstractSomalierModel.sample_name(self.trio.father.sample)
        mother = AbstractSomalierModel.sample_name(self.trio.mother.sample)
        proband_sex = Sex.UNKNOWN
        if patient := self.trio.proband.sample.patient:
            proband_sex = patient.sex
        write_trio_ped(filename, proband, proband_sex,
                       father, self.trio.father_affected, mother, self.trio.mother_affected)


class SomalierAllSamplesRelate(SomalierRelate):
    def get_sample_somalier_filenames(self) -> List[str]:
        cfg = SomalierConfig()
        return [f"{cfg['vcf_base_dir']}/**/*.somalier"]  # Wild card

    def get_samples(self) -> Iterable[Sample]:
        return Sample.objects.filter(import_status=ImportStatus.SUCCESS)


class SomalierRelatePairs(models.Model):
    relate = models.ForeignKey(SomalierAllSamplesRelate, on_delete=CASCADE)
    # Sample A always has a lower PK than B
    sample_a = models.ForeignKey(Sample, on_delete=CASCADE, related_name="somalierrelatepairs_a")
    sample_b = models.ForeignKey(Sample, on_delete=CASCADE, related_name="somalierrelatepairs_b")
    relatedness = models.FloatField()
    ibs0 = models.IntegerField()
    ibs2 = models.IntegerField()
    hom_concordance = models.FloatField()
    hets_a = models.IntegerField()
    hets_b = models.IntegerField()
    hets_ab = models.IntegerField()
    shared_hets = models.IntegerField()
    hom_alts_a = models.IntegerField()
    hom_alts_b = models.IntegerField()
    shared_hom_alts = models.IntegerField()
    n = models.IntegerField()
    x_ibs0 = models.IntegerField()
    x_ibs2 = models.IntegerField()

    class Meta:
        unique_together = ('sample_a', 'sample_b')


class SomalierConfig:
    def __init__(self):
        self.settings = settings.SOMALIER

    def _annotation_dir(self, dirname):
        return os.path.join(self.settings["annotation_base_dir"], dirname)

    def get_annotation(self, key):
        return self._annotation_dir(self.settings["annotation"][key])

    def report_dir(self, *args):
        return os.path.join(self.settings["report_base_dir"], *map(str, args))

    def ancestry_dir(self, subdir):
        return self.report_dir("ancestry", subdir)

    def related_dir(self, subdir):
        return self.report_dir("related", subdir)

    def get_sites(self, genome_build: 'GenomeBuild'):
        sites = self.settings["annotation"]["sites"][genome_build.name]
        return self._annotation_dir(sites)

    def get_sites_vcf(self, genome_build: 'GenomeBuild'):
        sites_name = os.path.basename(self.get_sites(genome_build))
        sites_vcf_kwargs = {"name": sites_name, "genome_build": genome_build}
        try:
            return VCF.objects.get(**sites_vcf_kwargs)
        except VCF.DoesNotExist as dne:
            print(f"Expected single VCF loaded via: {sites_vcf_kwargs}")
            raise dne

    def __getitem__(self, key):
        return self.settings[key]