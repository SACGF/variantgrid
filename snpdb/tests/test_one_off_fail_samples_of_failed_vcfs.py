from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase
from django.utils import timezone

from snpdb.models import VCF, GenomeBuild, ImportStatus, Sample


class FailSamplesOfFailedVCFsTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        user = User.objects.get_or_create(username="testuser_failed_vcfs")[0]
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.samples = {}
        for vcf_status, sample_status in [(ImportStatus.ERROR, ImportStatus.SUCCESS),
                                          (ImportStatus.ERROR, ImportStatus.MARKED_FOR_DELETION),
                                          (ImportStatus.SUCCESS, ImportStatus.SUCCESS)]:
            vcf = VCF.objects.create(name=f"{vcf_status}_{sample_status}", genome_build=grch37, user=user,
                                     genotype_samples=1, import_status=vcf_status, date=timezone.now())
            sample = Sample.objects.create(name=vcf.name, vcf=vcf, import_status=sample_status)
            cls.samples[(vcf_status, sample_status)] = sample.pk

    def _status(self, key):
        return Sample.objects.get(pk=self.samples[key]).import_status

    def test_only_open_samples_of_failed_vcfs_are_failed(self):
        call_command("one_off_fail_samples_of_failed_vcfs")

        self.assertEqual(self._status((ImportStatus.ERROR, ImportStatus.SUCCESS)), ImportStatus.ERROR)
        self.assertEqual(self._status((ImportStatus.ERROR, ImportStatus.MARKED_FOR_DELETION)),
                         ImportStatus.MARKED_FOR_DELETION)
        self.assertEqual(self._status((ImportStatus.SUCCESS, ImportStatus.SUCCESS)), ImportStatus.SUCCESS)
