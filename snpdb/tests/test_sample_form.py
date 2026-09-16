from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from snpdb.forms import SampleForm
from snpdb.models import VCF, GenomeBuild, ImportStatus, Sample


class TestSampleFormUploadedBy(TestCase):
    """ Sample has no owner of its own - it's derived from the VCF that uploaded it """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.uploader = User.objects.get_or_create(username="vcf_uploader")[0]
        cls.viewer = User.objects.get_or_create(username="sample_viewer")[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        vcf = VCF.objects.create(name="test_sample_form_vcf", genotype_samples=1, genome_build=genome_build,
                                 import_status=ImportStatus.SUCCESS, user=cls.uploader, date=timezone.now())
        cls.sample = Sample.objects.create(name="proband", vcf=vcf, import_status=ImportStatus.SUCCESS)

    def test_uploaded_by_initial_comes_from_vcf_user(self):
        form = SampleForm(user=self.viewer, instance=self.sample)
        self.assertEqual(form.fields["uploaded_by"].initial, self.uploader)
