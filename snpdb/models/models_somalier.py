import logging
import os
import uuid
from typing import List, Iterable

from django.conf import settings
from django.db import models
from django.db.models import CASCADE
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager

from library.utils import execute_cmd
from pedigree.ped.export_ped import write_unrelated_ped, write_trio_ped
from snpdb.models import Sample, VCF, Cohort, Trio, SuperPopulationCode
from snpdb.models.models_enums import ProcessingStatus
from patients.models_enums import Sex


class AbstractSomalierModel(TimeStampedModel):
    status = models.CharField(max_length=1, choices=ProcessingStatus.CHOICES, default=ProcessingStatus.CREATED)
    error_exception = models.TextField(null=True, blank=True)

    class Meta:
        abstract = True

    def get_samples(self) -> Iterable[Sample]:
        raise NotImplementedError()

    def get_sample_somalier_filenames(self) -> List[str]:
        return [AbstractSomalierModel.sample_filename(s) for s in self.get_samples()]

    def execute(self, command: List[str], **kwargs):
        """ Executes code and handles saving errors """

        logging.info('About to call %s', " ".join(command))

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
            raise RuntimeError(self.error_exception)

    @staticmethod
    def sample_filename(sample: Sample):
        vcf_dir = sample.vcf.somaliervcfextract.get_somalier_dir()
        return os.path.join(vcf_dir, sample.vcf_sample_name + ".somalier")


class SomalierVCFExtract(AbstractSomalierModel):
    vcf = models.OneToOneField(VCF, on_delete=CASCADE)

    def get_somalier_dir(self):
        cfg = SomalierConfig()
        return os.path.join(cfg["vcf_base_dir"], str(self.vcf.pk))

    def get_samples(self) -> Iterable[Sample]:
        return self.vcf.sample_set.all().order_by("pk")


class SomalierAncestryRun(AbstractSomalierModel):
    """ We do a run against a whole VCF """
    vcf_extract = models.OneToOneField(SomalierVCFExtract, on_delete=CASCADE)
    uuid = models.UUIDField(default=uuid.uuid4, editable=False)  # code to hide directories in media_root

    def get_report_dir(self):
        cfg = SomalierConfig()
        return cfg.ancestry_dir(self.uuid)

    def get_samples(self) -> Iterable[Sample]:
        return self.vcf_extract.get_samples()


class SomalierAncestry(TimeStampedModel):
    ancestry_run = models.ForeignKey(SomalierAncestryRun, on_delete=CASCADE)
    sample = models.OneToOneField(Sample, on_delete=CASCADE)
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
        raise NotImplementedError()


class SomalierCohortRelate(SomalierRelate):
    cohort = models.OneToOneField(Cohort, on_delete=CASCADE)
    cohort_version = models.IntegerField()

    def get_samples(self) -> Iterable[Sample]:
        return self.cohort.get_samples()

    def write_ped_file(self, filename):
        write_unrelated_ped(filename, [s.vcf_sample_name for s in self.get_samples()])


class SomalierTrioRelate(SomalierRelate):
    trio = models.OneToOneField(Trio, on_delete=CASCADE)

    def get_samples(self) -> Iterable[Sample]:
        return self.trio.get_samples()

    def write_ped_file(self, filename):
        proband_sample = self.trio.proband.sample
        proband = proband_sample.vcf_sample_name
        father = self.trio.father.sample.vcf_sample_name
        mother = self.trio.mother.sample.vcf_sample_name
        proband_sex = Sex.UNKNOWN
        if patient := proband_sample.patient:
            proband_sex = patient.sex
        write_trio_ped(filename, proband, proband_sex,
                       father, self.trio.father_affected, mother, self.trio.mother_affected)


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

    def __getitem__(self, key):
        return self.settings[key]