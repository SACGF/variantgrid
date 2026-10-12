from unittest import mock

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from ontology.models import OntologyImport, OntologyService, OntologyTerm, OntologyVersion
from ontology.tests.test_data_ontology import (
    create_ontology_test_data,
    create_test_ontology_version,
)
from patients.models import Patient
from patients.models.models_phenotype import (
    PatientTextPhenotype,
    PhenotypeMatchVersion,
    TextPhenotype,
    TextPhenotypeMatch,
    patient_phenotype_terms,
)
from patients.phenotype_matcher import (
    PHENOTYPE_MATCHER_VERSION,
    PhenotypeMatcher,
    _build_ambiguous_acronym_denylist,
    get_ambiguous_acronym_denylist,
)
from patients.phenotype_matching import (
    bulk_patient_phenotype_matching,
    create_phenotype_description,
    requeue_sentences,
)
from snpdb.models import Cohort, GenomeBuild


class TestPhenotypeMatching(TestCase):

    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()
        create_test_ontology_version()
        cls._create_common_word_test_terms()
        cls.phenotype_matcher = PhenotypeMatcher()

    @staticmethod
    def _create_common_word_test_terms():
        """ Real terms that a common word matches (variantgrid_com#60, #2125) """
        ontology_import, _ = OntologyImport.objects.get_or_create(import_source="test", filename="test_skip_words",
                                                                  context="test",
                                                                  defaults={"processed_date": timezone.now()})
        terms = [
            ("OMIM:614813", OntologyService.OMIM, "SHORT STATURE, ONYCHODYSPLASIA, FACIAL DYSMORPHISM, AND HYPOTRICHOSIS",
             ["SOFT", "SOFT SYNDROME"]),
            ("HP:0032198", OntologyService.HPO, "Decreased prothrombin time", ["Decreased INR"]),
            ("HP:0031915", OntologyService.HPO, "Stable", []),
            ("HP:0000256", OntologyService.HPO, "Macrocephaly", []),
            ("MONDO:0018971", OntologyService.MONDO, "isolated oxycephaly", ["acrocephaly"]),
            ("HP:0001903", OntologyService.HPO, "Anemia", ["Anaemia"]),
            ("MONDO:0002280", OntologyService.MONDO, "anemia (disease)", ["anemia"]),
            ("HP:0008151", OntologyService.HPO, "Prolonged prothrombin time", ["Prolonged PT"]),
            ("HP:0000010", OntologyService.HPO, "Recurrent urinary tract infections", ["Recurrent UTIs"]),
            ("HP:0002373", OntologyService.HPO, "Febrile seizure", ["Febrile seizures"]),
            ("HP:0012514", OntologyService.HPO, "Lower limb pain", ["Leg pain"]),
        ]
        for term_id, ontology_service, name, aliases in terms:
            OntologyTerm.objects.get_or_create(id=term_id, defaults={
                "ontology_service": ontology_service, "index": int(term_id.split(":")[1]), "name": name,
                "aliases": aliases, "from_import": ontology_import})

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

    def test_skip_words(self):
        """ "soft" is an exact alias of SOFT syndrome, "decreased in" is 1 edit from "Decreased INR" but in a
            2-letter word, which is an abbreviation rather than a typo """
        SKIP_WORDS = {"soft": None,
                      "decreased in": None,
                      "SOFT syndrome": (OntologyService.OMIM, "OMIM:614813"),
                      "Decreased INR": (OntologyService.HPO, "HP:0032198")}

        self.check_expected_results_for_description(SKIP_WORDS)

    def _match_ids(self, text, phenotype_matcher=None) -> set[str]:
        phenotype_matcher = phenotype_matcher or self.phenotype_matcher
        return set(phenotype_matcher.get_matches([(w, None) for w in text.split()]))

    def test_dictionary_word_not_fuzzy_matched(self):
        """ "table" is 1 edit from HPO "Stable" but is spelled correctly, so isn't a typo """
        self.assertEqual(self._match_ids("table"), set())
        self.assertEqual(self._match_ids("macrocephaky"), {"HP:0000256"})

    def test_fuzzy_match_is_one_typo_in_one_word(self):
        """ Without the special-case overrides, 1 edit in an abbreviation or a negation prefix isn't a typo #2130 """
        matcher = PhenotypeMatcher()
        matcher.hardcoded_lookups = {}
        matcher.case_insensitive_lookups = {}
        matcher.disease_families = {}
        expected_ids_by_phrase = {
            "prolonged qt": set(),
            "recurrent urtis": set(),
            "afebrile seizures": set(),
            "decreased in": set(),
            "leg pains": {"HP:0012514"},
            "febrile seizres": {"HP:0002373"},
            "recurrent urinary tract infectons": {"HP:0000010"},
        }
        for phrase, expected_ids in expected_ids_by_phrase.items():
            with self.subTest(phrase=phrase):
                self.assertEqual(self._match_ids(phrase, matcher), expected_ids)

    def test_exact_match_stops_fuzzy_match_in_other_ontologies(self):
        """ macrocephaly is 1 edit from MONDO "acrocephaly", but other spellings of the same term still match """
        self.assertEqual(self._match_ids("macrocephaly"), {"HP:0000256"})
        self.assertEqual(self._match_ids("anaemia"), {"HP:0001903", "MONDO:0002280"})
        self.assertEqual(self._match_ids("anaemias"), {"HP:0001903", "MONDO:0002280"})

    @override_settings(PATIENT_PHENOTYPE_EXCLUDE_STRING="----needs human review")
    def test_exclude_string_skips_persistence(self):
        phenotype = "Raised TSH\n----needs human review"
        patient = Patient(phenotype=phenotype)
        patient.save(phenotype_matcher=self.phenotype_matcher)
        self.assertFalse(PatientTextPhenotype.objects.filter(patient=patient).exists(),
                         "Exclude marker should prevent persisting matches")

    def test_matcher_only_built_when_there_is_new_text_to_match(self):
        """ PhenotypeMatcher loads the whole ontology so takes seconds to build - saving a patient
            shouldn't pay for that unless there are sentences we haven't already matched """
        already_matched_text = "Failure to thrive"
        Patient(phenotype=already_matched_text).save(phenotype_matcher=self.phenotype_matcher)

        with mock.patch("patients.phenotype_matching.PhenotypeMatcher") as mock_matcher:
            Patient(patient_code="no phenotype").save()
            Patient(phenotype=already_matched_text).save()
            mock_matcher.assert_not_called()

    def test_commas(self):
        COMMA_OMIM = {
            "PLATELET DISORDER, FAMILIAL, WITH ASSOCIATED MYELOID MALIGNANCY": (OntologyService.OMIM, "OMIM:601399"),
            "LEUKEMIA, ACUTE MYELOID": (OntologyService.OMIM, "OMIM:601626"),
        }

        self.check_expected_results_for_description(COMMA_OMIM)

    def test_ambiguous_acronym_flagged_and_excluded(self):
        """A phenotype text whose lowercased form is in the denylist should:
        - NOT create TextPhenotypeMatch rows (so downstream Django queries
          through PATIENT_TPM_PATH can't pick up the wrong concept)
        - emit a synthetic `ambiguous_alias` + `ambiguous_alias_candidates`
          entry on get_results() so the UI can list the conflicting concepts
        - be excluded from get_ontology_term_ids()."""
        denylist = {
            "failure to thrive": (
                ("HP:0001508", "Failure to thrive"),
                ("OMIM:000000", "Some other thing called FTT"),
            ),
        }
        with mock.patch(
            "patients.models.models_phenotype.get_ambiguous_acronym_denylist",
            return_value=denylist,
        ):
            # Rebuild the matcher inside the patch so its ambiguous_acronyms
            # picks up the patched denylist (the class-level matcher was built
            # in setUpTestData against the real ontology test data).
            matcher = PhenotypeMatcher()
            patient = Patient(phenotype="Failure to thrive")
            patient.save(phenotype_matcher=matcher)
            patient.process_phenotype_if_changed(phenotype_matcher=matcher)

            pd = patient.patient_text_phenotype.phenotype_description

            saved = TextPhenotypeMatch.objects.filter(
                text_phenotype__textphenotypesentence__phenotype_description=pd,
            )
            self.assertFalse(
                saved.exists(),
                "Ambiguous-acronym text must not produce TextPhenotypeMatch rows",
            )

            results = pd.get_results()
            self.assertTrue(results, "Expected a synthetic warning result for 'Failure to thrive'")
            flagged = [r for r in results if r.get("ambiguous_alias")]
            self.assertTrue(flagged, f"Expected ambiguous_alias flag on results: {results}")
            self.assertEqual(
                flagged[0].get("ambiguous_alias_candidates"),
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

    def _create_matched_patient(self, phenotype: str) -> Patient:
        patient = Patient(phenotype=phenotype)
        patient.save(phenotype_matcher=self.phenotype_matcher)
        patient.process_phenotype_if_changed(phenotype_matcher=self.phenotype_matcher)
        return patient

    def test_patient_phenotype_terms_agrees_with_single_patient_path(self):
        patient = self._create_matched_patient("Raised TSH")
        no_phenotype = Patient(patient_code="no phenotype")
        no_phenotype.save()

        phenotype_terms = patient_phenotype_terms([patient, no_phenotype])
        self.assertNotIn(no_phenotype.pk, phenotype_terms, "Patient without phenotype text is absent")
        terms = phenotype_terms[patient.pk]
        self.assertEqual(terms.text, "Raised TSH")
        term_ids = [term.pk for terms_for_service in terms.terms.values() for term in terms_for_service]
        self.assertEqual(sorted(term_ids), patient.get_ontology_term_ids())
        self.assertEqual([(t["id"], t["match_text"]) for t in terms.to_json()["terms"]["HPO"]],
                         [("HP:0002925", "Raised TSH")])

    def test_patient_phenotype_terms_excludes_ambiguous_acronyms(self):
        """ Rows matched before a term joined the denylist stay in the DB - both paths drop them """
        patient = self._create_matched_patient("Raised TSH")
        denylist = {
            "raised tsh": (("HP:0002925", "Raised TSH"), ("OMIM:000000", "Something else")),
        }
        with mock.patch(
            "patients.models.models_phenotype.get_ambiguous_acronym_denylist",
            return_value=denylist,
        ):
            phenotype_terms = patient_phenotype_terms([patient])
            phenotype_description = patient.patient_text_phenotype.phenotype_description
            phenotype_description.get_ontology_term_ids.invalidate(phenotype_description)
            self.assertEqual(phenotype_description.get_ontology_term_ids(), [])

        self.assertEqual(phenotype_terms[patient.pk].terms, {},
                         "Ambiguous acronym matches must not become terms")

    def test_hardcoded_override_wins_over_denylist(self):
        """If a key has a hardcoded lookup (e.g. FTT), the public denylist
        accessor must filter it out so the match isn't falsely flagged."""
        # Raw includes "ftt"; effective denylist should not.
        raw = {
            "ftt": (("HP:0001508", "Failure to thrive"),),
            "some_truly_ambiguous_token": (("HP:0000001", "All"), ("MONDO:0000001", "disease")),
        }
        with mock.patch(
            "patients.phenotype_matcher._build_ambiguous_acronym_denylist",
            return_value=raw,
        ):
            # bust cache_memoize so our patch is used
            _build_ambiguous_acronym_denylist.invalidate(0)
            effective = get_ambiguous_acronym_denylist()
            self.assertNotIn("ftt", effective, "FTT has a HARDCODED_LOOKUP - must not be flagged")
            self.assertIn("some_truly_ambiguous_token", effective)


class TestPhenotypeMatcherVersion(TestCase):
    """ Sentences record what they were matched with, so stale ones can be rematched in place (#2131) """

    @classmethod
    def setUpTestData(cls):
        create_ontology_test_data()
        cls.ontology_version = create_test_ontology_version()
        cls.phenotype_matcher = PhenotypeMatcher()

    def _stamp(self, text, matcher_version, ontology_version):
        match_version, _ = PhenotypeMatchVersion.objects.get_or_create(matcher_version=matcher_version,
                                                                       ontology_version=ontology_version)
        TextPhenotype.objects.filter(text=text).update(match_version=match_version)

    def test_matched_sentences_are_stamped(self):
        create_phenotype_description("Raised TSH", self.phenotype_matcher)
        create_phenotype_description("...")  # Nothing to match, so marked processed without a matcher
        for text in ["Raised TSH", "..."]:
            text_phenotype = TextPhenotype.objects.get(text=text)
            self.assertTrue(text_phenotype.processed)
            self.assertEqual(text_phenotype.match_version.matcher_version, PHENOTYPE_MATCHER_VERSION)
            self.assertEqual(text_phenotype.match_version.ontology_version, self.ontology_version)
        self.assertEqual(PhenotypeMatchVersion.objects.count(), 1)
        self.assertFalse(TextPhenotype.stale_qs().exists())

    def test_stale_qs(self):
        for text in ["current", "old matcher", "old ontology", "never stamped"]:
            TextPhenotype.objects.create(text=text, processed=True)
        TextPhenotype.objects.create(text="unprocessed")
        imports = {f: getattr(self.ontology_version, f) for f in OntologyVersion.ONTOLOGY_IMPORTS}
        imports["gencc_import"] = OntologyImport.objects.create(import_source="test", filename="older_gencc",
                                                                processed_date=timezone.now())
        other_ontology_version = OntologyVersion.objects.create(**imports)
        self._stamp("current", PHENOTYPE_MATCHER_VERSION, self.ontology_version)
        self._stamp("old matcher", PHENOTYPE_MATCHER_VERSION - 1, self.ontology_version)
        self._stamp("old ontology", PHENOTYPE_MATCHER_VERSION, other_ontology_version)
        stale = set(TextPhenotype.stale_qs().values_list("text", flat=True))
        self.assertEqual(stale, {"old matcher", "old ontology", "never stamped"})

    def test_requeue_keeps_links_and_approvals(self):
        text = "Raised TSH"
        user = User.objects.get_or_create(username="phenotype_approver")[0]
        patient = Patient(phenotype=text)
        patient.save(phenotype_matcher=self.phenotype_matcher, phenotype_approval_user=user)
        cohort = Cohort(name="phenotype cohort", user=user, genome_build=GenomeBuild.get_name_or_alias("GRCh37"),
                        phenotype=text)
        cohort.save(phenotype_matcher=self.phenotype_matcher)
        patient_description = patient.phenotype_description
        cohort_description = cohort.phenotype_description
        expected_term_ids = patient_description.get_ontology_term_ids()
        self.assertTrue(expected_term_ids)

        TextPhenotype.objects.filter(text=text).update(match_version=None)  # Matched before #2131
        self.assertEqual(requeue_sentences(TextPhenotype.stale_qs()), 1)
        self.assertEqual(patient_description.get_ontology_term_ids(), [])  # memo invalidated

        bulk_patient_phenotype_matching(patients=[patient])

        text_phenotype = TextPhenotype.objects.get(text=text)
        self.assertEqual(text_phenotype.match_version, PhenotypeMatchVersion.get_or_create_current())
        patient_text_phenotype = PatientTextPhenotype.objects.get(patient=patient)
        self.assertEqual(patient_text_phenotype.phenotype_description, patient_description)
        self.assertEqual(patient_text_phenotype.approved_by, user)
        self.assertEqual(Cohort.objects.get(pk=cohort.pk).phenotype_description, cohort_description)
        self.assertEqual(patient_description.get_ontology_term_ids(), expected_term_ids)
