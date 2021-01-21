from django.test import TestCase


from annotation.phenotype_matching import default_lookup_factory, \
    cached_lookup_factory
from annotation.tests.test_data_mim_hpo import create_mim_hpo_test_data
from patients.models import Patient


class TestPhenotypeMatching(TestCase):

    @classmethod
    def setUpTestData(cls):
        create_mim_hpo_test_data()

        lookups = default_lookup_factory()
        cls.lookup_factory = cached_lookup_factory(*lookups)

    def create_patient_match_phenotypes(self, phenotype):
        patient = Patient.objects.create(phenotype=phenotype)
        patient.process_phenotype_if_changed(lookup_factory=self.lookup_factory)

        return patient.patient_text_phenotype.phenotype_description.get_results()

    def check_expected_results_for_description(self, expected_results_by_description):
        for acronym, expected_results in expected_results_by_description.items():
            results = self.create_patient_match_phenotypes(acronym)
            if expected_results:
                expected_match_type, expected_pk = expected_results
                if results:
                    result = results[0]
                    self.assertEqual(result["match_type"], expected_match_type, "Match type")
                    self.assertEqual(result["pk"], expected_pk, "Match PK")
                else:
                    self.fail(f"No results, expected: {expected_results}")
            else:
                self.assertListEqual(results, [], "Empty results")

    def test_acronyms(self):
        HARDCODED_LOOKUPS = {"FTT": (PhenotypeMatchTypes.HPO, 1508),
                             "LQTS": (PhenotypeMatchTypes.HPO, 1657)}

        self.check_expected_results_for_description(HARDCODED_LOOKUPS)

    def test_aliases(self):
        CASE_INSENSITIVE_LOOKUPS = {"Raised TSH": (PhenotypeMatchTypes.HPO, 2925),
                                    "MEN type 1": (PhenotypeMatchTypes.OMIM, 131100),
                                    "Hypoplastic right ventricle": (PhenotypeMatchTypes.HPO, 4762),
                                    "bowel cancer": (PhenotypeMatchTypes.OMIM, 114500),
                                    "bowel polyps": (PhenotypeMatchTypes.HPO, 200063)}

        self.check_expected_results_for_description(CASE_INSENSITIVE_LOOKUPS)

    def test_mispellings(self):
        TYPOS = {"Maplem Syrup Urine Disease": (PhenotypeMatchTypes.OMIM, 248600),
                 "IMERSLUND-GRSBECK SYNDROME": (PhenotypeMatchTypes.OMIM, 53)}  # Alias

        self.check_expected_results_for_description(TYPOS)

    def test_syndrome(self):
        SYNDROME_ABBREV = {"IMERSLUND-GRSBECK SYNDROME": (PhenotypeMatchTypes.OMIM, 53),  # Alias
                           "IMERSLUND-GRSBECK SYNDROMES": (PhenotypeMatchTypes.OMIM, 53),  # Alias
                           "IMERSLUND-GRSBECK SYND": (PhenotypeMatchTypes.OMIM, 53),  # Alias
                           "IMERSLUND-GRSBECK SYND.": (PhenotypeMatchTypes.OMIM, 53)}  # Alias

        self.check_expected_results_for_description(SYNDROME_ABBREV)

    def test_commas(self):
        COMMA_OMIM = {
            "PLATELET DISORDER, FAMILIAL, WITH ASSOCIATED MYELOID MALIGNANCY": (PhenotypeMatchTypes.OMIM, 601399),
            "LEUKEMIA, ACUTE MYELOID": (PhenotypeMatchTypes.OMIM, 79929),
        }

        self.check_expected_results_for_description(COMMA_OMIM)
