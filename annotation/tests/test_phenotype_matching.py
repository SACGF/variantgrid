from unittest import mock

from django.test import TestCase, override_settings

from annotation.models.models_phenotype_match import PatientTextPhenotype
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

    @override_settings(PATIENT_PHENOTYPE_EXCLUDE_STRING="----needs human review")
    def test_exclude_string_skips_persistence(self):
        phenotype = "Raised TSH\n----needs human review"
        patient = Patient(phenotype=phenotype)
        patient.save(phenotype_matcher=self.phenotype_matcher)
        self.assertFalse(PatientTextPhenotype.objects.filter(patient=patient).exists(),
                         "Exclude marker should prevent persisting matches")

    def test_commas(self):
        COMMA_OMIM = {
            "PLATELET DISORDER, FAMILIAL, WITH ASSOCIATED MYELOID MALIGNANCY": (OntologyService.OMIM, "OMIM:601399"),
            "LEUKEMIA, ACUTE MYELOID": (OntologyService.OMIM, "OMIM:601626"),
        }

        self.check_expected_results_for_description(COMMA_OMIM)

    def test_ambiguous_acronym_flagged_and_excluded(self):
        """A phenotype text whose lowercased form is in the denylist should:
        - emit `ambiguous_alias` and `ambiguous_alias_candidates` on the
          to_dict() payload so the UI can list the conflicting concepts
        - be excluded from get_ontology_term_ids() so downstream gene-list
          computation doesn't silently use the wrong term."""
        denylist = {
            "failure to thrive": (
                ("HP:0001508", "Failure to thrive"),
                ("OMIM:000000", "Some other thing called FTT"),
            ),
        }
        with mock.patch(
            "annotation.models.models_phenotype_match.get_ambiguous_acronym_denylist",
            return_value=denylist,
        ):
            patient = Patient(phenotype="Failure to thrive")
            patient.save(phenotype_matcher=self.phenotype_matcher)
            patient.process_phenotype_if_changed(phenotype_matcher=self.phenotype_matcher)

            pd = patient.patient_text_phenotype.phenotype_description
            results = pd.get_results()
            self.assertTrue(results, "Expected at least one match for 'Failure to thrive'")
            flagged = [r for r in results if r.get("ambiguous_alias")]
            self.assertTrue(flagged, f"Expected ambiguous_alias flag on results: {results}")
            candidates = flagged[0].get("ambiguous_alias_candidates")
            self.assertEqual(
                candidates,
                [
                    {"accession": "HP:0001508", "name": "Failure to thrive"},
                    {"accession": "OMIM:000000", "name": "Some other thing called FTT"},
                ],
                "Expected the conflicting candidates to be exposed for UI display",
            )

            # Cached function - bust the per-instance cache by clearing
            pd.get_ontology_term_ids.invalidate(pd)
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
        raw = {
            "ftt": (("HP:0001508", "Failure to thrive"),),
            "some_truly_ambiguous_token": (("HP:0000001", "All"), ("MONDO:0000001", "disease")),
        }
        with mock.patch(
            "annotation.phenotype_matcher._build_ambiguous_acronym_denylist",
            return_value=raw,
        ):
            # bust cache_memoize so our patch is used
            _build_ambiguous_acronym_denylist.invalidate(0)
            effective = get_ambiguous_acronym_denylist()
            self.assertNotIn("ftt", effective, "FTT has a HARDCODED_LOOKUP - must not be flagged")
            self.assertIn("some_truly_ambiguous_token", effective)
