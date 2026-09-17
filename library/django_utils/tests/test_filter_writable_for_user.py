"""
filter_writable_for_user() is the batch form of can_write() - the grids' delete columns resolve a
whole page through it. These check the two agree wherever a model overrides one of them.
"""
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls.base import reverse
from django.utils import timezone
from guardian.shortcuts import assign_perm

from analysis.models import Analysis, AnalysisLock, VariantTag
from annotation.fake_annotation import create_fake_variants
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import ExternalModelManager, ExternalPK, Patient
from snpdb.models import VCF, Cohort, GenomeBuild, ImportStatus, Sample, Tag, Trio, Variant
from snpdb.tests.utils.fake_cohort_data import create_fake_trio


class FilterWritableForUserTest(TestCase):
    """ The owner can write their own objects, another user can't """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.get_or_create(username="writable_owner")[0]
        cls.other = User.objects.get_or_create(username="writable_other")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.trio = create_fake_trio(cls.owner, cls.grch37)
        create_fake_variants(cls.grch37)
        cls.sample = cls.trio.get_samples()[0]
        cls.vcf = cls.sample.vcf
        cls.analysis = Analysis.objects.create(name="writable analysis", user=cls.owner,
                                               genome_build=cls.grch37)
        assign_permission_to_user_and_groups(cls.owner, cls.analysis)

    def assert_agrees_with_can_write(self, klass, obj, user, expected: bool):
        self.assertEqual(obj.can_write(user), expected, f"{klass.__name__}.can_write for {user}")
        writable_pks = set(klass.filter_writable_for_user(user).values_list("pk", flat=True))
        self.assertEqual(obj.pk in writable_pks, expected,
                         f"{klass.__name__}.filter_writable_for_user for {user}")

    def _tag_analysis(self, tag_id: str) -> VariantTag:
        tag = Tag.objects.get_or_create(pk=tag_id)[0]
        return VariantTag.objects.create(variant=Variant.objects.first(), tag=tag,
                                         analysis=self.analysis, user=self.owner,
                                         genome_build=self.grch37)

    def test_vcf(self):
        self.assert_agrees_with_can_write(VCF, self.vcf, self.owner, True)
        self.assert_agrees_with_can_write(VCF, self.vcf, self.other, False)

    def test_sample_through_vcf(self):
        """ Sample.can_write falls back to the VCF's permission """
        self.assert_agrees_with_can_write(Sample, self.sample, self.owner, True)
        self.assert_agrees_with_can_write(Sample, self.sample, self.other, False)

    def test_trio_through_cohort(self):
        """ Trio delegates its permissions to its cohort """
        self.assert_agrees_with_can_write(Trio, self.trio, self.owner, True)
        self.assert_agrees_with_can_write(Trio, self.trio, self.other, False)

    def test_cohort_behind_vcf(self):
        """ A VCF's cohort takes the VCF's permission """
        self.assert_agrees_with_can_write(Cohort, self.vcf.cohort, self.owner, True)
        self.assert_agrees_with_can_write(Cohort, self.vcf.cohort, self.other, False)

    def test_analysis(self):
        self.assert_agrees_with_can_write(Analysis, self.analysis, self.owner, True)
        self.assert_agrees_with_can_write(Analysis, self.analysis, self.other, False)

    def test_locked_analysis(self):
        """ Analysis.can_write refuses a locked analysis even to its owner """
        AnalysisLock.objects.create(analysis=self.analysis, locked=True, user=self.owner,
                                    date=timezone.now())
        self.assert_agrees_with_can_write(Analysis, self.analysis, self.owner, False)

    def test_unlocked_analysis(self):
        """ The last lock set is the current status """
        AnalysisLock.objects.create(analysis=self.analysis, locked=True, user=self.owner,
                                    date=timezone.now())
        AnalysisLock.objects.create(analysis=self.analysis, locked=False, user=self.owner,
                                    date=timezone.now())
        self.assert_agrees_with_can_write(Analysis, self.analysis, self.owner, True)

    def test_patient(self):
        patient = Patient.objects.create(patient_code="writable patient")
        assign_permission_to_user_and_groups(self.owner, patient)
        self.assert_agrees_with_can_write(Patient, patient, self.owner, True)
        self.assert_agrees_with_can_write(Patient, patient, self.other, False)

    def test_externally_managed_patient(self):
        """ A record an external system owns is read only, whatever Guardian says """
        manager = ExternalModelManager.objects.create(name="writable test manager", can_modify=False)
        external_pk = ExternalPK.objects.create(code="EXT1", external_type="patient",
                                                external_manager=manager)
        patient = Patient.objects.create(patient_code="external patient", external_pk=external_pk)
        assign_permission_to_user_and_groups(self.owner, patient)
        self.assert_agrees_with_can_write(Patient, patient, self.owner, False)

    def test_variant_tag_through_analysis(self):
        """ VariantTag.can_write delegates to the analysis when the tag was made in one """
        variant_tag = self._tag_analysis("writable-test-tag")
        self.assert_agrees_with_can_write(VariantTag, variant_tag, self.owner, True)
        self.assert_agrees_with_can_write(VariantTag, variant_tag, self.other, False)

    def test_variant_tag_locked_analysis(self):
        """ ...including the analysis' lock """
        variant_tag = self._tag_analysis("writable-locked-tag")
        AnalysisLock.objects.create(analysis=self.analysis, locked=True, user=self.owner,
                                    date=timezone.now())
        self.assert_agrees_with_can_write(VariantTag, variant_tag, self.owner, False)


class DeleteColumnPermissionTest(TestCase):
    """ The delete column offers a link only where the user can write the row """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.get_or_create(username="delete_column_owner")[0]
        cls.reader = User.objects.get_or_create(username="delete_column_reader")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        create_fake_trio(cls.owner, cls.grch37)
        VCF.objects.all().update(import_status=ImportStatus.SUCCESS)
        cls.vcf = VCF.objects.get()
        assign_perm(VCF.get_read_perm(), cls.reader, cls.vcf)

    def _delete_links(self, user) -> list:
        self.client.force_login(user)
        response = self.client.get(reverse("vcfs_datatable"), {"draw": 1, "start": 0, "length": 100})
        self.assertEqual(response.status_code, 200)
        return [row["delete"] for row in response.json()["data"]]

    def test_owner_gets_delete_link(self):
        self.assertEqual(self._delete_links(self.owner),
                         [reverse('group_permissions_object_delete',
                                  kwargs={'class_name': 'snpdb.models.models_vcf.VCF',
                                          'primary_key': self.vcf.pk})])

    def test_read_only_user_gets_none(self):
        self.assertEqual(self._delete_links(self.reader), [None])
