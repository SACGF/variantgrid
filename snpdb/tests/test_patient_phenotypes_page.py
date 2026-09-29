"""  The patient phenotype column on the VCF / cohort pages and the row on the sample page (issue #1426).

     Owns: that the page data covers the page's samples' patients in one query and leaves out patients
     the user cannot view.
"""
from django.contrib.auth.models import User
from django.db import connection
from django.test import TestCase
from django.test.utils import CaptureQueriesContext
from django.urls import reverse

from annotation.phenotype_matcher import PhenotypeMatcher
from library.guardian_utils import assign_permission_to_user_and_groups
from ontology.tests.test_data_ontology import (
    create_ontology_test_data,
    create_test_ontology_version,
)
from patients.models import Patient
from snpdb.fake_data import create_fake_cohort
from snpdb.models import GenomeBuild
from snpdb.views.vcf_cohort_page import vcf_cohort_page_context


class PatientPhenotypesPageTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()
        create_test_ontology_version()
        phenotype_matcher = PhenotypeMatcher()

        cls.user = User.objects.get_or_create(username='phenotype_page_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, cls.grch37)
        cls.vcf = cls.cohort.vcf

        cls.patient = Patient(phenotype="Raised TSH")
        cls.patient.save(phenotype_matcher=phenotype_matcher)
        assign_permission_to_user_and_groups(cls.user, cls.patient)

        cls.other_lab_patient = Patient(phenotype="Raised TSH")
        cls.other_lab_patient.save(phenotype_matcher=phenotype_matcher)

        proband, mother = cls.cohort.get_samples()[:2]
        proband.patient = cls.patient
        proband.save()
        mother.patient = cls.other_lab_patient
        mother.save()

    def test_vcf_page_data_covers_visible_patients_in_one_query(self):
        with CaptureQueriesContext(connection) as ctx:
            context = vcf_cohort_page_context(self.user, self.cohort, True, vcf=self.vcf)
        patient_phenotypes = context["patient_phenotypes"]

        self.assertEqual(list(patient_phenotypes), [self.patient.pk],
                         "Only patients the user can view are in the page data")
        hpo_terms = patient_phenotypes[self.patient.pk]["terms"]["HPO"]
        self.assertEqual([t["id"] for t in hpo_terms], ["HP:0002925"])

        phenotype_queries = [q for q in ctx.captured_queries if "annotation_textphenotypematch" in q["sql"]]
        self.assertEqual(len(phenotype_queries), 1, "Phenotype terms come from one query for the page")

    def test_sample_page_data(self):
        sample = self.cohort.get_samples()[0]
        self.client.force_login(self.user)
        response = self.client.get(reverse("view_sample", kwargs={"sample_id": sample.pk}))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(list(response.context["patient_phenotypes"]), [self.patient.pk])
