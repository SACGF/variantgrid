"""  The cohort's own phenotype text (issue #1426).

     Owns: that saving a Cohort matches its phenotype text to ontology terms through CohortTextPhenotype,
     and that the cohort and VCF pages carry the editor and the seed-from-patients button.
"""
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from annotation.models.models_phenotype_match import CohortTextPhenotype
from annotation.phenotype_matcher import PhenotypeMatcher
from library.guardian_utils import assign_permission_to_user_and_groups
from ontology.tests.test_data_ontology import (
    create_ontology_test_data,
    create_test_ontology_version,
)
from patients.models import Patient
from snpdb.models import Cohort, GenomeBuild, ImportStatus
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort

RAISED_TSH = "Raised TSH"
RAISED_TSH_HPO = "HP:0002925"


class CohortPhenotypeTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()
        create_test_ontology_version()

        cls.user = User.objects.get_or_create(username='cohort_phenotype_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.vcf_cohort = create_fake_cohort(cls.user, cls.grch37)
        cls.vcf = cls.vcf_cohort.vcf

        cls.cohort = Cohort.objects.create(name="custom cohort", user=cls.user, genome_build=cls.grch37,
                                           import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.user, cls.cohort)

    def test_saving_phenotype_text_matches_ontology_terms(self):
        cohort = Cohort.objects.get(pk=self.cohort.pk)
        cohort.phenotype = RAISED_TSH
        cohort.save(phenotype_matcher=PhenotypeMatcher(), phenotype_approval_user=self.user)

        link = CohortTextPhenotype.objects.get(cohort=cohort)
        self.assertEqual(link.approved_by, self.user)
        self.assertEqual(link.phenotype_description.original_text, RAISED_TSH)
        self.assertIn(RAISED_TSH_HPO, cohort.get_ontology_term_ids())

    def test_cohort_page_saves_the_phenotype(self):
        self.client.force_login(self.user)
        url = reverse("view_cohort", kwargs={"cohort_id": self.cohort.pk})
        response = self.client.post(url, {"name": "custom cohort", "phenotype": RAISED_TSH})
        self.assertEqual(response.status_code, 200)

        cohort = Cohort.objects.get(pk=self.cohort.pk)
        self.assertEqual(cohort.phenotype, RAISED_TSH)
        self.assertEqual(CohortTextPhenotype.objects.get(cohort=cohort).approved_by, self.user)

    def test_vcf_page_edits_its_cohorts_phenotype(self):
        self.client.force_login(self.user)
        response = self.client.get(reverse("view_vcf", kwargs={"vcf_id": self.vcf.pk}))
        self.assertEqual(response.status_code, 200)
        # Prefixed so it doesn't collide with the create patient dialog's own phenotype field
        self.assertEqual(response.context["cohort_phenotype_form"].instance, self.vcf_cohort)
        self.assertContains(response, 'id="id_cohort-phenotype"')

    def test_seed_button_needs_a_visible_patient_with_phenotype_text(self):
        self.client.force_login(self.user)
        url = reverse("view_vcf", kwargs={"vcf_id": self.vcf.pk})
        self.assertNotContains(self.client.get(url), 'id="add-sample-patient-phenotypes"')

        patient = Patient(phenotype=RAISED_TSH)
        patient.save(phenotype_matcher=PhenotypeMatcher())
        assign_permission_to_user_and_groups(self.user, patient)
        sample = self.vcf_cohort.get_samples()[0]
        sample.patient = patient
        sample.save()

        self.assertContains(self.client.get(url), 'id="add-sample-patient-phenotypes"')
