import os

from django.conf import settings
from django.db import models
from django.db.models import CASCADE
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager

from snpdb.models import Sample, VCF, Cohort, Trio, SuperPopulationCode


class SomalierVCFExtract(TimeStampedModel):
    vcf = models.OneToOneField(VCF, on_delete=CASCADE)

    def get_dir(self):
        somalier_cfg = settings.SOMALIER
        return os.path.join(somalier_cfg["vcf_base_dir"], self.vcf.pk)

    @staticmethod
    def sample_filename(sample: Sample):
        vcf_dir = sample.vcf.somaliervcfextract.get_dir()
        return os.path.join(vcf_dir, sample.vcf_sample_name + ".somalier")


class SomalierAncestryRun(TimeStampedModel):
    """ We do a run against a whole VCF """
    vcf_extract = models.ForeignKey(SomalierVCFExtract, on_delete=CASCADE)


class SomalierAncestry(TimeStampedModel):
    ancestry_run = models.ForeignKey(SomalierAncestryRun, on_delete=CASCADE)
    sample = models.OneToOneField(Sample, on_delete=CASCADE)
    predicted_ancestry = models.CharField(max_length=1, choices=SuperPopulationCode.choices)
    EAS_prob = models.FloatField()
    AFR_prob = models.FloatField()
    AMR_prob = models.FloatField()
    SAS_prob = models.FloatField()
    EUR_prob = models.FloatField()


class SomalierRelate(TimeStampedModel):
    objects = InheritanceManager()
    uuid = models.UUIDField()  # code to hide directories in media_root

    class Meta:
        abstract = True


class SomalierCohortRelate(SomalierRelate):
    cohort = models.OneToOneField(Cohort, on_delete=CASCADE)
    cohort_version = models.IntegerField()


class SomalierTrioRelate(SomalierRelate):
    trio = models.OneToOneField(Trio, on_delete=CASCADE)
