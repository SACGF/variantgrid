""" analysis/related_analyses.py via the related_analyses_for_samples tag - what a sample or patient page lists """
from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import Analysis
from analysis.models.nodes.sources import CohortNode, SampleNode
from analysis.models.nodes.sources.duo_node import DuoNode
from analysis.models.nodes.sources.quad_node import QuadNode
from analysis.templatetags.related_analyses_tags import related_analyses_for_samples
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Patient
from snpdb.enums import ImportStatus
from snpdb.fake_data import create_fake_cohort, create_fake_duo, create_fake_quad
from snpdb.models import Cohort, CohortSample, Country, GenomeBuild, Lab, Organization


class RelatedAnalysesForSamplesTest(TestCase):
    """ A sample's page lists analyses built on a duo or quad the sample is a member of (#1889) """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='related_analyses_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.duo = create_fake_duo(cls.user, cls.grch37)
        cls.quad = create_fake_quad(cls.user, cls.grch37)

    def _analysis(self) -> Analysis:
        analysis = Analysis(genome_build=self.grch37)
        analysis.set_defaults_and_save(self.user)
        return analysis

    def test_shows_duo_and_quad_analyses(self):
        duo_node = DuoNode.objects.create(analysis=self._analysis(), duo=self.duo)
        quad_node = QuadNode.objects.create(analysis=self._analysis(), quad=self.quad)
        samples = [self.duo.relative.sample, self.quad.sibling.sample]
        context = related_analyses_for_samples({"user": self.user}, samples, show_sample_info=True)
        self.assertEqual(dict(context["analysis_details"]),
                         {duo_node.analysis: f"Duo: {self.duo}", quad_node.analysis: f"Quad: {self.quad}"})


class RelatedAnalysesForSamplesPermissionsTest(TestCase):
    """ A patient's samples can come from another lab's VCF - the table only names those the user can view """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        organization = Organization.objects.get_or_create(name="Fake Org", group_name="fake_org")[0]
        country = Country.objects.get_or_create(name="Australia")[0]
        cohorts = []
        for lab_name in ["lab_a", "lab_b"]:
            lab = Lab.objects.get_or_create(name=lab_name, city="Adelaide", country=country,
                                            organization=organization, group_name=f"fake_org/{lab_name}")[0]
            user = User.objects.get_or_create(username=f"related_analyses_{lab_name}_user")[0]
            lab.group.user_set.add(user)
            cohorts.append(create_fake_cohort(user, cls.grch37, name=lab_name))
        cls.user, other_lab_user = [cohort.user for cohort in cohorts]
        cls.cohort, cls.other_lab_cohort = cohorts

        cls.patient = Patient.objects.create(first_name="RELATED", last_name="ANALYSES")
        cls.sample = cls.cohort.cohortsample_set.get(sample__name="proband").sample
        cls.other_lab_sample = cls.other_lab_cohort.cohortsample_set.get(sample__name="proband").sample
        for sample in [cls.sample, cls.other_lab_sample]:
            sample.patient = cls.patient
            sample.save()

        # The other lab's own cohort built over this lab's sample
        cls.other_lab_sub_cohort = Cohort.objects.create(name="lab_b_sub_cohort", user=other_lab_user,
                                                         genome_build=cls.grch37, import_status=ImportStatus.SUCCESS)
        CohortSample.objects.create(cohort=cls.other_lab_sub_cohort, sample=cls.sample,
                                    cohort_genotype_packed_field_index=0, sort_order=0)
        assign_permission_to_user_and_groups(other_lab_user, cls.other_lab_sub_cohort)

    def _analysis(self) -> Analysis:
        analysis = Analysis(genome_build=self.grch37)
        analysis.set_defaults_and_save(self.user)
        return analysis

    def test_only_samples_and_cohorts_user_can_view(self):
        sample_node = SampleNode.objects.create(analysis=self._analysis(), sample=self.sample)
        other_lab_analysis = self._analysis()
        SampleNode.objects.create(analysis=other_lab_analysis, sample=self.other_lab_sample)
        CohortNode.objects.create(analysis=other_lab_analysis, cohort=self.other_lab_cohort)
        CohortNode.objects.create(analysis=other_lab_analysis, cohort=self.other_lab_sub_cohort)

        context = related_analyses_for_samples({"user": self.user}, self.patient.get_samples(), show_sample_info=True)
        self.assertEqual(context["samples"], [self.sample])
        self.assertEqual(dict(context["analysis_details"]), {sample_node.analysis: "Sample: proband"})
