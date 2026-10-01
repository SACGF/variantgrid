from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from guardian.shortcuts import assign_perm

from library.django_utils.unittest_utils import prevent_request_warnings
from patients.models import Patient
from snpdb.fake_data import create_fake_cohort
from snpdb.models.models_enums import SampleFileType
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_vcf import VCF, SampleFilePath
from snpdb.sample_file_path import (
    create_sample_file_paths,
    resolve_sample_file_paths,
    validate_sample_file_path_pattern,
)


class TestVCFSampleFilePaths(TestCase):

    @classmethod
    def setUpTestData(cls):
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.user = User.objects.get_or_create(username="vcf_sample_file_paths_owner")[0]
        cls.vcf = create_fake_cohort(cls.user, grch37, name="vcf_sample_file_paths").vcf
        cls.samples = {s.name: s for s in cls.vcf.sample_set.all()}
        for name, sample in cls.samples.items():
            sample.vcf_sample_name = f"{name}_vcf"
            sample.save()
        cls.samples["proband"].patient = Patient.objects.create(patient_code="PAT-1")
        cls.samples["proband"].save()
        cls.samples["mother"].patient = Patient.objects.create(first_name="no", last_name="code")
        cls.samples["mother"].save()

    def _paths_and_errors(self, pattern) -> dict:
        resolutions = resolve_sample_file_paths(self.vcf, pattern, SampleFileType.BAM)
        return {r.sample.name: (r.file_path, r.error) for r in resolutions}

    def test_resolve_pattern(self):
        self.assertEqual(self._paths_and_errors("/data/%(vcf_sample_name)s.bam"), {
            "proband": ("/data/proband_vcf.bam", None),
            "mother": ("/data/mother_vcf.bam", None),
            "father": ("/data/father_vcf.bam", None),
        })
        self.assertEqual(self._paths_and_errors("/data/%(patient_code)s.bam"), {
            "proband": ("/data/PAT-1.bam", None),
            "mother": (None, "No patient_code"),
            "father": (None, "No patient_code"),
        })
        for bad_pattern in ["/data/sample.bam", "/data/%(unknown)s.bam", "/data/%(sample)d.bam", "/data/%(sample)s%"]:
            with self.assertRaises(ValueError, msg=bad_pattern):
                validate_sample_file_path_pattern(bad_pattern)
        validate_sample_file_path_pattern("/data/100%%/%(sample_id)05d.bam")

    def test_create_is_additive(self):
        proband = self.samples["proband"]
        mother = self.samples["mother"]
        SampleFilePath.objects.create(sample=proband, file_type=SampleFileType.BAM, label="old",
                                      file_path="/data/proband_vcf.bam")
        SampleFilePath.objects.create(sample=mother, file_type=SampleFileType.BED, file_path="/data/mother_vcf.bam")

        pattern = "/data/%(vcf_sample_name)s.bam"
        for _ in range(2):
            resolutions = resolve_sample_file_paths(self.vcf, pattern, SampleFileType.BAM)
            create_sample_file_paths(resolutions, SampleFileType.BAM, "new")

        proband_files = list(proband.samplefilepath_set.values_list("file_type", "label"))
        self.assertEqual(proband_files, [(SampleFileType.BAM, "old")])
        mother_files = set(mother.samplefilepath_set.values_list("file_type", "label"))
        self.assertEqual(mother_files, {(SampleFileType.BED, None), (SampleFileType.BAM, "new")})
        self.assertEqual(self.samples["father"].samplefilepath_set.count(), 1)

    @prevent_request_warnings
    def test_write_permission_required(self):
        reader = User.objects.create(username="vcf_sample_file_paths_reader")
        assign_perm(VCF.get_read_perm(), reader, self.vcf)
        self.client.force_login(reader)

        proband = self.samples["proband"]
        vcf_url = reverse("vcf_sample_files_tab", kwargs={"vcf_id": self.vcf.pk})
        sample_url = reverse("sample_files_tab", kwargs={"sample_id": proband.pk})
        for url in [vcf_url, sample_url]:
            self.assertEqual(self.client.get(url).status_code, 200, msg=url)

        data = {"pattern": "/data/%(vcf_sample_name)s.bam", "file_type": SampleFileType.BAM, "action": "save"}
        response = self.client.post(vcf_url, data)
        self.assertEqual(response.status_code, 403)

        formset_data = {
            "samplefilepath_set-TOTAL_FORMS": 1, "samplefilepath_set-INITIAL_FORMS": 0,
            "samplefilepath_set-0-file_type": SampleFileType.BAM, "samplefilepath_set-0-file_path": "/data/x.bam",
        }
        response = self.client.post(sample_url, formset_data)
        self.assertEqual(response.status_code, 403)
        self.assertFalse(SampleFilePath.objects.filter(sample__vcf=self.vcf).exists())
