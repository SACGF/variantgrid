from django.test import TestCase

from annotation.models.models_phenotype_match import PhenotypeMatchTypes
from annotation.phenotype_matching import default_lookup_factory, cached_lookup_factory
from ontology.tests.test_data_ontology import create_ontology_test_data
from patients.models import Patient


class TestPhenotypeMatching(TestCase):

    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()

        lookups = default_lookup_factory()
        cls.lookup_factory = cached_lookup_factory(*lookups)

    def create_patient_match_phenotypes(self, phenotype):
        patient = Patient.objects.create(phenotype=phenotype)
        patient.process_phenotype_if_changed(lookup_factory=self.lookup_factory)

        return patient.patient_text_phenotype.phenotype_description.get_results()

    def check_expected_results_for_description(self, expected_results_by_description):
        for phenotype, expected_results in expected_results_by_description.items():
            results = self.create_patient_match_phenotypes(phenotype)
            if expected_results:
                expected_match_type, expected_pk = expected_results
                if results:
                    result = results[0]
                    self.assertEqual(result["match_type"], expected_match_type, "Match type")
                    self.assertEqual(result["pk"], expected_pk, "Match PK")
                else:
                    self.fail(f"No results for '{phenotype}', expected: {expected_results}")
            else:
                self.assertListEqual(results, [], "Empty results")

    def test_acronyms(self):
        HARDCODED_LOOKUPS = {"FTT": (PhenotypeMatchTypes.HPO, "HP:0001508"),
                             "LQTS": (PhenotypeMatchTypes.HPO, "HP:0001657")}

        self.check_expected_results_for_description(HARDCODED_LOOKUPS)

    def test_aliases(self):
        CASE_INSENSITIVE_LOOKUPS = {"Raised TSH": (PhenotypeMatchTypes.HPO, "HP:0002925"),
                                    "MEN type 1": (PhenotypeMatchTypes.OMIM, "OMIM:131100"),
                                    "Hypoplastic right ventricle": (PhenotypeMatchTypes.HPO, "HP:0004762"),
                                    "bowel polyps": (PhenotypeMatchTypes.HPO, "HP:0200063")}

        self.check_expected_results_for_description(CASE_INSENSITIVE_LOOKUPS)

    def test_mispellings(self):
        TYPOS = {"Maplem Syrup Urine Disease": (PhenotypeMatchTypes.OMIM, "OMIM:248600"),
                 "XENTEROCYTE COBALAMIN MALABSORPTION": (PhenotypeMatchTypes.OMIM, "OMIM:261100")}  # Alias

        self.check_expected_results_for_description(TYPOS)

    def test_syndrome(self):
        SYNDROME_ABBREV = {"IMERSLUND-GRSBECK SYNDROME 1": (PhenotypeMatchTypes.OMIM, "OMIM:261100"),  # Alias
                           "IMERSLUND-GRSBECK SYNDROMES 1": (PhenotypeMatchTypes.OMIM, "OMIM:261100"),  # Alias
                           "IMERSLUND-GRSBECK SYND 1": (PhenotypeMatchTypes.OMIM, "OMIM:261100")}  # Alias

        self.check_expected_results_for_description(SYNDROME_ABBREV)

    def test_commas(self):
        COMMA_OMIM = {
            "PLATELET DISORDER, FAMILIAL, WITH ASSOCIATED MYELOID MALIGNANCY": (PhenotypeMatchTypes.OMIM, "OMIM:601399"),
            "LEUKEMIA, ACUTE MYELOID": (PhenotypeMatchTypes.OMIM, "OMIM:601626"),
        }

        self.check_expected_results_for_description(COMMA_OMIM)
