""" The VCF locus filter editor panel (analysis/views/views_node.py:vcf_locus_filters) - who may load which VCF """
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from guardian.shortcuts import assign_perm

from analysis.models import Analysis, TrioNode
from analysis.models.enums import TrioInheritance
from annotation.fake_data import get_fake_annotation_version
from library.django_utils.unittest_utils import prevent_request_warnings
from library.guardian_utils import DjangoPermission
from snpdb.fake_data import create_fake_trio
from snpdb.models import GenomeBuild
from snpdb.models.models_vcf import VCF


class TestVCFLocusFilterViews(TestCase):

    @classmethod
    def setUpTestData(cls):
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)
        owner = User.objects.create_user(username="vlf_owner")
        cls.viewer = User.objects.create_user(username="vlf_viewer")

        cls.trio = create_fake_trio(owner, grch37, name="vlf")
        cls.other_trio = create_fake_trio(owner, grch37, name="vlf_other")

        cls.analysis = Analysis(genome_build=grch37)
        cls.analysis.set_defaults_and_save(owner)
        assign_perm(DjangoPermission.perm(cls.analysis, DjangoPermission.READ), cls.viewer, cls.analysis)
        cls.node = TrioNode.objects.create(analysis=cls.analysis, trio=cls.trio, inheritance=TrioInheritance.RECESSIVE)

        # Analysis owner has left - every one of their permissions now fails (#2325)
        owner.is_active = False
        owner.save()

    def setUp(self):
        self.client.force_login(self.viewer)

    def _vcf_url(self, vcf: VCF) -> str:
        return reverse("vcf_locus_filters", kwargs={"analysis_id": self.analysis.pk, "node_id": self.node.pk,
                                                    "vcf_id": vcf.pk})

    def test_nodes_own_vcf_comes_with_the_analysis(self):
        vcf = self.trio.cohort.vcf
        self.assertFalse(self.viewer.has_perm(DjangoPermission.perm(vcf, DjangoPermission.READ), vcf))

        response = self.client.get(self._vcf_url(vcf))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "YOUSHALLNOTPASS")

        cohort_url = reverse("cohort_vcf_locus_filters", kwargs={"analysis_id": self.analysis.pk,
                                                                 "node_id": self.node.pk,
                                                                 "cohort_id": self.trio.cohort.pk})
        self.assertEqual(self.client.get(cohort_url).status_code, 200)

    @prevent_request_warnings
    def test_any_other_vcf_is_checked_as_the_analysis_owner(self):
        other_vcf = self.other_trio.cohort.vcf
        assign_perm(DjangoPermission.perm(other_vcf, DjangoPermission.READ), self.viewer, other_vcf)
        self.assertEqual(self.client.get(self._vcf_url(other_vcf)).status_code, 403)

        User.objects.filter(pk=self.analysis.user_id).update(is_active=True)
        self.assertEqual(self.client.get(self._vcf_url(other_vcf)).status_code, 200)
