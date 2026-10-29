""" A user who can only read a VCF, its samples or a tag config collection can't change them through
the views that edit them """
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone
from guardian.shortcuts import assign_perm

from genes.models import SampleGeneList
from library.guardian_utils import DjangoPermission
from patients.models import Patient
from snpdb.enums import ImportStatus
from snpdb.forms import SampleForm
from snpdb.models import (
    VCF,
    GenomeBuild,
    Sample,
    TagConfigCollection,
    VariantZygosityCountCollection,
)


class VCFWritePermissionTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user(username="vcf_write_owner")
        cls.reader = User.objects.create_user(username="vcf_write_reader")
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.vcf = VCF.objects.create(name="vcf_write_test", genotype_samples=1, genome_build=genome_build,
                                     import_status=ImportStatus.SUCCESS, user=cls.owner, date=timezone.now())
        cls.sample = Sample.objects.create(name="proband", vcf=cls.vcf, import_status=ImportStatus.SUCCESS)
        for user in (cls.owner, cls.reader):
            assign_perm(DjangoPermission.perm(cls.vcf, DjangoPermission.READ), user, cls.vcf)
        assign_perm(DjangoPermission.perm(cls.vcf, DjangoPermission.WRITE), cls.owner, cls.vcf)

    def test_view_vcf_post_needs_write(self):
        self.client.force_login(self.reader)
        url = reverse("view_vcf", kwargs={"vcf_id": self.vcf.pk})
        response = self.client.post(url, {"name": "renamed"})
        self.assertEqual(response.status_code, 403)
        self.vcf.refresh_from_db()
        self.assertEqual(self.vcf.name, "vcf_write_test")

    def test_change_zygosity_count_needs_write_and_known_operation(self):
        vzcc = VariantZygosityCountCollection.objects.create(name="vcf_write_test_vzcc")

        def change_zygosity_count(user, operation):
            self.client.force_login(user)
            kwargs = {"vcf_id": self.vcf.pk, "vzcc_id": vzcc.pk, "operation": operation}
            return self.client.post(reverse("vcf_change_zygosity_count", kwargs=kwargs))

        self.assertEqual(change_zygosity_count(self.reader, "add").status_code, 403)
        self.assertEqual(change_zygosity_count(self.owner, "multiply").status_code, 400)

    def test_vcf_being_deleted_is_not_found(self):
        self.vcf.import_status = ImportStatus.MARKED_FOR_DELETION
        self.vcf.save()
        self.client.force_login(self.owner)
        response = self.client.get(reverse("view_vcf", kwargs={"vcf_id": self.vcf.pk}))
        self.assertEqual(response.status_code, 404)

    def test_sample_gene_list_post_needs_write(self):
        self.client.force_login(self.reader)
        url = reverse("sample_gene_lists_tab", kwargs={"sample_id": self.sample.pk})
        response = self.client.post(url, {"name": "Sample Gene List", "custom_gene_list_text": "BRCA1"})
        self.assertEqual(response.status_code, 403)
        self.assertFalse(SampleGeneList.objects.filter(sample=self.sample).exists())

    def test_sample_patient_choices_are_readable_or_current(self):
        visible = Patient.objects.create(patient_code="visible")
        hidden = Patient.objects.create(patient_code="hidden")
        current = Patient.objects.create(patient_code="current")
        assign_perm(DjangoPermission.perm(visible, DjangoPermission.READ), self.owner, visible)
        self.sample.patient = current
        self.sample.save()

        form = SampleForm(user=self.owner, instance=self.sample)
        patients = set(form.fields["patient"].queryset.filter(pk__in=[visible.pk, hidden.pk, current.pk]))
        self.assertEqual(patients, {visible, current})


class SetTagColorTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user(username="tag_color_owner")
        cls.reader = User.objects.create_user(username="tag_color_reader")
        cls.tcc = TagConfigCollection.objects.create(name="tag_color_test", user=cls.owner)
        assign_perm(DjangoPermission.perm(cls.tcc, DjangoPermission.READ), cls.reader, cls.tcc)

    def _set_color(self, user, rgb, tag="ResultTagColor"):
        self.client.force_login(user)
        url = reverse("set_tag_color", kwargs={"tag_config_collection_id": self.tcc.pk})
        return self.client.post(url, {"tag": tag, "rgb": rgb})

    def test_needs_write(self):
        self.assertEqual(self._set_color(self.reader, "#ff0000").status_code, 403)
        self.assertFalse(self.tcc.tagconfig_set.exists())

    def test_rgb_must_be_hex_colour(self):
        self.assertEqual(self._set_color(self.owner, "red;}").status_code, 400)
        self.assertFalse(self.tcc.tagconfig_set.exists())

    def test_unknown_tag_not_found(self):
        self.assertEqual(self._set_color(self.owner, "#ff0000", tag="no_such_tag_2512").status_code, 404)
        self.assertFalse(self.tcc.tagconfig_set.exists())
