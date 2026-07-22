from unittest.mock import patch

from django.test import TestCase, override_settings

from mme.serializers.patient_profile import (
    classification_ontology_slots,
    classification_genomic_feature,
    build_patient_profile,
)
from mme.tests.fakes import (
    FakeClassification, FakeSubmission, FakeCoordinate, FakeVariant, make_term,
)
from ontology.models import OntologyService


@override_settings(MME_ONTOLOGY_SNAKE_EXACT=True, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False,
                   MME_CONTACT={"name": "Test", "href": "mailto:t@t.org"},
                   MME_ENABLED_PRODUCTION_SUBMIT=False)
class PatientProfileTestCase(TestCase):

    def test_hpo_condition_term_routes_to_features(self):
        term = make_term("HP:0001250", OntologyService.HPO, 1250, "Seizure")
        classification = FakeClassification(terms=[term])
        features, disorders = classification_ontology_slots(classification)
        self.assertEqual([f["id"] for f in features], ["HP:0001250"])
        self.assertEqual(features[0]["observed"], "yes")
        self.assertEqual(disorders, [])

    def test_omim_condition_term_routes_to_disorders(self):
        term = make_term("OMIM:143100", OntologyService.OMIM, 143100, "Huntington disease")
        classification = FakeClassification(terms=[term])
        features, disorders = classification_ontology_slots(classification)
        self.assertEqual(disorders, [{"id": "MIM:143100"}])
        self.assertEqual(features, [])

    def test_mondo_condition_term_crosswalked_to_omim(self):
        mondo = make_term("MONDO:0007739", OntologyService.MONDO, 7739, "Huntington disease")
        omim = make_term("OMIM:143100", OntologyService.OMIM, 143100, "Huntington disease")
        classification = FakeClassification(terms=[mondo])
        with patch("mme.serializers.patient_profile.OntologyTermRelation.as_omim", return_value=omim):
            features, disorders = classification_ontology_slots(classification)
        self.assertEqual(disorders, [{"id": "MIM:143100"}])
        self.assertEqual(features, [])

    def test_mixed_ontology_condition_splits_across_slots(self):
        hpo = make_term("HP:0001250", OntologyService.HPO, 1250, "Seizure")
        omim = make_term("OMIM:143100", OntologyService.OMIM, 143100, "Huntington disease")
        classification = FakeClassification(terms=[hpo, omim])
        features, disorders = classification_ontology_slots(classification)
        self.assertEqual([f["id"] for f in features], ["HP:0001250"])
        self.assertEqual(disorders, [{"id": "MIM:143100"}])

    def test_genomic_feature_coordinate_zero_based_conversion(self):
        variant = FakeVariant(FakeCoordinate("17", 43044295, "A", "G"))
        classification = FakeClassification(gene_symbol="BRCA1", variant=variant)
        gf = classification_genomic_feature(classification)
        self.assertEqual(gf[0]["gene"], {"id": "BRCA1"})
        # our 1-based 43044295 -> MME 0-based 43044294
        self.assertEqual(gf[0]["variant"]["start"], 43044294)
        self.assertEqual(gf[0]["variant"]["assembly"], "GRCh38")
        self.assertEqual(gf[0]["variant"]["referenceName"], "17")
        self.assertEqual(gf[0]["variant"]["referenceBases"], "A")
        self.assertEqual(gf[0]["variant"]["alternateBases"], "G")

    def test_genomic_feature_strips_chr_prefix(self):
        variant = FakeVariant(FakeCoordinate("chr17", 100, "A", "G"))
        classification = FakeClassification(gene_symbol="BRCA1", variant=variant)
        gf = classification_genomic_feature(classification)
        self.assertEqual(gf[0]["variant"]["referenceName"], "17")

    def test_genomic_feature_symbolic_variant_gene_only(self):
        variant = FakeVariant(FakeCoordinate("17", 100, "N", "<DEL>", is_symbolic=True))
        classification = FakeClassification(gene_symbol="BRCA1", variant=variant)
        gf = classification_genomic_feature(classification)
        self.assertEqual(gf[0]["gene"], {"id": "BRCA1"})
        self.assertNotIn("variant", gf[0])

    def test_genomic_feature_none_when_no_gene_or_variant(self):
        self.assertIsNone(classification_genomic_feature(FakeClassification()))

    def test_classification_without_patient_still_submits(self):
        # condition term + variant, no linked sample/patient
        term = make_term("OMIM:143100", OntologyService.OMIM, 143100, "Huntington disease")
        variant = FakeVariant(FakeCoordinate("4", 3074877, "C", "T"))
        classification = FakeClassification(terms=[term], gene_symbol="HTT", variant=variant)
        submission = FakeSubmission(classification=classification)
        profile = build_patient_profile(submission)
        self.assertIn("disorders", profile)
        self.assertIn("genomicFeatures", profile)
        self.assertEqual(profile["id"], "vg:1")

    def test_build_profile_raises_without_features_or_genomic(self):
        # a lone disorder is not enough: MME requires features or genomicFeatures
        term = make_term("OMIM:143100", OntologyService.OMIM, 143100, "Huntington disease")
        classification = FakeClassification(terms=[term])   # no gene, no variant, no HPO
        submission = FakeSubmission(classification=classification)
        with self.assertRaises(ValueError):
            build_patient_profile(submission)

    def test_test_flag_present_until_production_submit_enabled(self):
        term = make_term("HP:0001250", OntologyService.HPO, 1250, "Seizure")
        submission = FakeSubmission(classification=FakeClassification(terms=[term]))
        profile = build_patient_profile(submission)
        self.assertTrue(profile.get("test"))

    @override_settings(MME_ENABLED_PRODUCTION_SUBMIT=True)
    def test_test_flag_absent_when_production_submit_enabled(self):
        term = make_term("HP:0001250", OntologyService.HPO, 1250, "Seizure")
        submission = FakeSubmission(classification=FakeClassification(terms=[term]))
        profile = build_patient_profile(submission)
        self.assertNotIn("test", profile)
