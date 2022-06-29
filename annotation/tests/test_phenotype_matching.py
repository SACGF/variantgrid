from django.test import TestCase

from annotation.phenotype_matcher import PhenotypeMatcher
from ontology.models import OntologyService
from ontology.tests.test_data_ontology import create_ontology_test_data, create_test_ontology_version
from patients.models import Patient


class TestPhenotypeMatching(TestCase):

    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()
        create_test_ontology_version()
        cls.phenotype_matcher = PhenotypeMatcher()

    def create_patient_match_phenotypes(self, phenotype):
        patient = Patient(phenotype=phenotype)
        patient.save(phenotype_matcher=self.phenotype_matcher)

        patient.process_phenotype_if_changed(phenotype_matcher=self.phenotype_matcher)
        return patient.patient_text_phenotype.phenotype_description.get_results()

    def check_expected_results_for_description(self, expected_results_by_description):
        for phenotype, expected_results in expected_results_by_description.items():
            results = self.create_patient_match_phenotypes(phenotype)
            if expected_results:
                expected_ontology_service, expected_pk = expected_results
                if results:
                    result = results[0]
                    self.assertEqual(result["ontology_service"], expected_ontology_service, "Ontology Service")
                    self.assertEqual(result["pk"], expected_pk, "Match PK")
                else:
                    self.fail(f"No results for '{phenotype}', expected: {expected_results}")
            else:
                self.assertListEqual(results, [], "Empty results")

    def test_acronyms(self):
        HARDCODED_LOOKUPS = {"FTT": (OntologyService.HPO, "HP:0001508")}

        self.check_expected_results_for_description(HARDCODED_LOOKUPS)

    def test_aliases(self):
        CASE_INSENSITIVE_LOOKUPS = {"Raised TSH": (OntologyService.HPO, "HP:0002925"),
                                    "MEN type 1": (OntologyService.OMIM, "OMIM:131100"),
                                    "Hypoplastic right ventricle": (OntologyService.HPO, "HP:0004762"),
                                    "bowel polyps": (OntologyService.HPO, "HP:0200063")}

        self.check_expected_results_for_description(CASE_INSENSITIVE_LOOKUPS)

    def test_mispellings(self):
        TYPOS = {"Maplem Syrup Urine Disease": (OntologyService.OMIM, "OMIM:248600"),
                 "XENTEROCYTE COBALAMIN MALABSORPTION": (OntologyService.OMIM, "OMIM:261100")}  # Alias

        self.check_expected_results_for_description(TYPOS)

    def test_syndrome(self):
        SYNDROME_ABBREV = {"IMERSLUND-GRSBECK SYNDROME 1": (OntologyService.OMIM, "OMIM:261100"),  # Alias
                           "IMERSLUND-GRSBECK SYNDROMES 1": (OntologyService.OMIM, "OMIM:261100"),  # Alias
                           "IMERSLUND-GRSBECK SYND 1": (OntologyService.OMIM, "OMIM:261100")}  # Alias

        self.check_expected_results_for_description(SYNDROME_ABBREV)

    def test_commas(self):
        COMMA_OMIM = {
            "PLATELET DISORDER, FAMILIAL, WITH ASSOCIATED MYELOID MALIGNANCY": (OntologyService.OMIM, "OMIM:601399"),
            "LEUKEMIA, ACUTE MYELOID": (OntologyService.OMIM, "OMIM:601626"),
        }

        self.check_expected_results_for_description(COMMA_OMIM)
