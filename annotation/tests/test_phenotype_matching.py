from unittest import mock

from django.test import TestCase

from annotation.models.models_phenotype_match import PhenotypeMatchTypes
from annotation.phenotype_matcher import PhenotypeMatcher
from ontology.tests.test_data_ontology import create_ontology_test_data
from patients.models import Patient


class TestPhenotypeMatching(TestCase):

    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()
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
        HARDCODED_LOOKUPS = {"FTT": (PhenotypeMatchTypes.HPO, "HP:0001508")}

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

    def test_ambiguous_acronym_flagged_and_excluded(self):
        """A phenotype text whose lowercased form is in the denylist should:
        - emit `ambiguous_alias` on the to_dict() payload
        - be excluded from get_ontology_term_ids() so downstream gene-list
          computation doesn't silently use the wrong term."""
        denylist = frozenset({"failure to thrive"})
        with mock.patch(
            "annotation.phenotype_matcher.get_ambiguous_acronym_denylist",
            return_value=denylist,
        ):
            patient = Patient(phenotype="Failure to thrive")
            patient.save(phenotype_matcher=self.phenotype_matcher)
            patient.process_phenotype_if_changed(phenotype_matcher=self.phenotype_matcher)

            pd = patient.patient_text_phenotype.phenotype_description
            results = pd.get_results()
            self.assertTrue(results, "Expected at least one match for 'Failure to thrive'")
            self.assertTrue(
                any(r.get("ambiguous_alias") for r in results),
                f"Expected ambiguous_alias flag on results: {results}",
            )

            term_ids = list(pd.get_ontology_term_ids())
            self.assertEqual(
                term_ids, [],
                "Ambiguous-acronym matches must be excluded from get_ontology_term_ids",
            )

    def test_hardcoded_override_wins_over_denylist(self):
        """If a key has a hardcoded lookup (e.g. FTT), the public denylist
        accessor must filter it out so the match isn't falsely flagged."""
        from annotation.phenotype_matcher import (
            _build_ambiguous_acronym_denylist,
            get_ambiguous_acronym_denylist,
        )
        # Raw includes "ftt"; effective denylist should not.
        with mock.patch(
            "annotation.phenotype_matcher._build_ambiguous_acronym_denylist",
            return_value=frozenset({"ftt", "some_truly_ambiguous_token"}),
        ):
            # bust cache_memoize so our patch is used
            _build_ambiguous_acronym_denylist.invalidate(0)
            effective = get_ambiguous_acronym_denylist()
            self.assertNotIn("ftt", effective, "FTT has a HARDCODED_LOOKUP - must not be flagged")
            self.assertIn("some_truly_ambiguous_token", effective)
