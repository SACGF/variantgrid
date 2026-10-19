from datetime import date

from django.contrib.auth.models import User
from django.test import TestCase

from classification.enums import SpecialEKeys, SubmissionSource
from classification.models import (
    Classification,
    ClassificationConsensus,
    ClassificationModification,
    ConsensusComparison,
)
from classification.models.classification import COPY_SCOPES_ALL, COPY_SCOPES_GENE
from classification.models.classification_variant_fields_validation import (
    apply_somatic_tier_from_amp_level,
)
from classification.models.evidence_key import resolve_default_value
from classification.tests.models.test_utils import ClassificationTestUtils
from genes.fake_data import create_fake_gene_version
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.models.models_enums import AlleleOriginBucket


class ClassificationConsensusPatchTestCase(TestCase):
    """ What copy_scope and copy_allele_origin let through when copying from a previous classification """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _consensus_patch(self, allele_origin_bucket: AlleleOriginBucket, **kwargs) -> dict:
        classification = Classification.objects.create(
            lab=self.lab,
            user=self.user,
            lab_record_id=f"consensus_{allele_origin_bucket}",
            allele_origin_bucket=allele_origin_bucket
        )
        modification = ClassificationModification.objects.create(
            classification=classification,
            user=self.user,
            source=SubmissionSource.API,
            is_last_published=True,
            published_evidence={
                "condition": {"value": "MONDO:0007947"},  # ALLELE
                "h_summary": {"value": "gene summary"},  # GENE
                "segregation": {"value": "co-segregates"},  # GERMLINE
                "somatic:tmb_value": {"value": 12},  # NONE - this patient's tumour
            }
        )
        return ClassificationConsensus(modification=modification, **kwargs).consensus_patch

    def test_germline_source_copies_everything_but_the_tumour_measurement(self):
        patch = self._consensus_patch(AlleleOriginBucket.GERMLINE)
        self.assertEqual(patch["condition"], {"value": "MONDO:0007947"})
        self.assertEqual(patch["h_summary"], {"value": "gene summary"})
        self.assertEqual(patch["segregation"], {"value": "co-segregates"})
        self.assertNotIn("somatic:tmb_value", patch)

    def test_somatic_source_leaves_germline_only_keys_behind(self):
        patch = self._consensus_patch(AlleleOriginBucket.SOMATIC)
        self.assertEqual(patch["condition"], {"value": "MONDO:0007947"})
        self.assertEqual(patch["h_summary"], {"value": "gene summary"})
        self.assertNotIn("segregation", patch)
        self.assertNotIn("somatic:tmb_value", patch)

    def test_gene_scope_copy_leaves_the_allele_level_content_behind(self):
        """ The gene box copies what is the same for every variant in the gene, and nothing else """
        patch = self._consensus_patch(AlleleOriginBucket.SOMATIC, copy_scopes=COPY_SCOPES_GENE)
        self.assertEqual(patch["h_summary"], {"value": "gene summary"})
        self.assertNotIn("condition", patch)


class GeneConsensusGroupsTestCase(TestCase):
    """ The deduplicated gene content candidates the create page, the classify dialog and the in-form box share """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        create_fake_gene_version(genome_build, "ENSG00000159216", "RUNX1", AnnotationConsortium.ENSEMBL)

    def _classify(self, lab_record_id: str, allele_origin: str, h_summary: str, curation_date: str,
                  condition: str = "Acute myeloid leukaemia") -> Classification:
        classification = Classification.create(
            user=self.user, lab=self.lab, lab_record_id=lab_record_id, source=SubmissionSource.VARIANT_GRID,
            data={
                SpecialEKeys.GENE_SYMBOL: "RUNX1",
                SpecialEKeys.ALLELE_ORIGIN: allele_origin,
                SpecialEKeys.CURATION_DATE: curation_date,
                SpecialEKeys.CONDITION: condition,
                "h_summary": h_summary,
                "interpretation_summary": f"about {lab_record_id}",  # variant level, so not part of the grouping
            })
        classification.publish_latest(self.user)
        return classification

    def _groups(self, allele_origin_bucket=AlleleOriginBucket.SOMATIC):
        return ClassificationConsensus.gene_consensus_groups(gene_symbol="RUNX1", user=self.user,
                                                             allele_origin_bucket=allele_origin_bucket)

    def test_records_with_the_same_gene_content_are_one_row(self):
        self._classify("older_copy", "somatic", "RUNX1 is a transcription factor", "2024-01-01")
        newest = self._classify("newer_copy", "somatic", "RUNX1 is a transcription factor", "2024-06-01")
        self._classify("rewritten", "somatic", "RUNX1 drives myeloid differentiation", "2024-03-01")

        groups = self._groups()
        self.assertEqual([group.record_count for group in groups], [2, 1])
        self.assertEqual(groups[0].representative.classification, newest)
        self.assertEqual(groups[0].other_count, 1)

    def test_the_most_used_content_comes_first_and_later_rows_name_what_differs(self):
        self._classify("copy_1", "somatic", "RUNX1 is a transcription factor", "2024-01-01")
        self._classify("copy_2", "somatic", "RUNX1 is a transcription factor", "2024-02-01")
        self._classify("reworded_newest", "somatic", "RUNX1 is a core transcription factor", "2024-06-01")

        groups = self._groups()
        self.assertEqual([group.record_count for group in groups], [2, 1])
        self.assertEqual([value.label for value in groups[0].content], ["Summary"])
        self.assertEqual(groups[0].content[0].value, "RUNX1 is a transcription factor")
        self.assertEqual(groups[0].differs_in, [])
        self.assertEqual(groups[1].differs_in, ["Summary"])

    def test_the_conditions_spread_is_distinct_and_most_used_first(self):
        self._classify("aml_1", "somatic", "RUNX1 is a transcription factor", "2024-01-01", condition="AML")
        self._classify("mds", "somatic", "RUNX1 is a transcription factor", "2024-02-01", condition="MDS")
        self._classify("aml_2", "somatic", "RUNX1 is a transcription factor", "2024-03-01", condition="AML")

        [group] = self._groups()
        self.assertEqual([condition.references[0].name for condition in group.conditions], ["AML", "MDS"])

    def test_the_other_bucket_is_not_offered(self):
        self._classify("germline_record", "germline", "RUNX1 is a transcription factor", "2024-01-01")
        self.assertEqual(self._groups(AlleleOriginBucket.SOMATIC), [])
        self.assertEqual(len(self._groups(AlleleOriginBucket.GERMLINE)), 1)

    def test_a_record_with_no_gene_content_is_not_a_candidate(self):
        self._classify("no_gene_content", "somatic", "", "2024-01-01")
        self.assertEqual(self._groups(), [])


class ConsensusComparisonTestCase(TestCase):
    """ The side by side view of the allele's previous classifications - what a pick would copy, with the fields the
        choices differ in separated from the ones they all carry """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _modification(self, lab_record_id: str, curation_date: str, **evidence) -> ClassificationModification:
        classification = Classification.create(
            user=self.user, lab=self.lab, lab_record_id=lab_record_id, source=SubmissionSource.VARIANT_GRID,
            data={SpecialEKeys.ALLELE_ORIGIN: "somatic", SpecialEKeys.CURATION_DATE: curation_date, **evidence})
        classification.publish_latest(self.user)
        return classification.last_published_version

    def _comparison(self, *modifications) -> ConsensusComparison:
        keys = ClassificationConsensus.copy_keys(COPY_SCOPES_ALL, AlleleOriginBucket.SOMATIC)
        return ConsensusComparison.of(ClassificationConsensus.consensus_groups(modifications, keys))

    def test_records_differing_only_in_uncopied_keys_are_one_choice_with_everything_shared(self):
        same_content = {"h_summary": "RUNX1 is a transcription factor", "condition": "MONDO:0018874"}
        comparison = self._comparison(
            self._modification("a", "2024-01-01", **same_content, **{"somatic:tmb_value": 12}),
            self._modification("b", "2024-02-01", **same_content, **{"somatic:tmb_value": 30}),
        )
        self.assertEqual(len(comparison.groups), 1)
        self.assertEqual(comparison.groups[0].record_count, 2)
        self.assertEqual(comparison.differing, [])
        shared_keys = {row.key for row in comparison.shared}
        self.assertTrue({"h_summary", "condition"} <= shared_keys)
        self.assertNotIn("somatic:tmb_value", shared_keys)

    def test_the_fields_the_choices_differ_in_come_apart_from_the_shared_ones(self):
        comparison = self._comparison(
            self._modification("a", "2024-01-01", h_summary="RUNX1 is a transcription factor", condition="MONDO:0018874"),
            self._modification("b", "2024-02-01", h_summary="RUNX1 drives myeloid differentiation", condition="MONDO:0018874"),
        )
        self.assertEqual(len(comparison.groups), 2)
        self.assertEqual([row.key for row in comparison.differing], ["h_summary"])
        self.assertEqual([cell.value for cell in comparison.differing[0].cells],
                         ["RUNX1 drives myeloid differentiation", "RUNX1 is a transcription factor"])
        self.assertIn("condition", {row.key for row in comparison.shared})

    def test_a_field_only_one_choice_has_is_a_difference_with_an_empty_cell(self):
        comparison = self._comparison(
            self._modification("a", "2024-01-01", h_summary="RUNX1 is a transcription factor"),
            self._modification("b", "2024-02-01", h_summary="RUNX1 is a transcription factor", condition="MONDO:0018874"),
        )
        [row] = comparison.differing
        self.assertEqual(row.key, "condition")
        self.assertEqual([cell is None for cell in row.cells], [False, True])


class ResolveDefaultValueTestCase(TestCase):
    """ The tokens a lab's classification_config can use for an evidence key's default_value """

    def test_tokens_resolve_and_anything_else_is_left_alone(self):
        user = User(username='defaulty')
        self.assertEqual(resolve_default_value("$user", user), "defaulty")

        user.first_name = "Dave"
        user.last_name = "Lawrence"
        self.assertEqual(resolve_default_value("$user", user), "Dave Lawrence")

        self.assertEqual(resolve_default_value("$today", user), date.today().isoformat())
        self.assertEqual(resolve_default_value("Labby curation", user), "Labby curation")


class CreateWithDefaultsTestCase(TestCase):
    """ What a web created record starts with when the lab's config gives a key a default_value """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.user.first_name = "Dave"
        self.user.last_name = "Lawrence"
        self.user.save()
        self.lab.classification_config = {
            SpecialEKeys.CURATED_BY: {"mandatory": True, "default_value": "$user"},
            SpecialEKeys.CURATION_DATE: {"mandatory": True, "default_value": "$today"}
        }
        self.lab.save()

    def _create(self, **kwargs) -> Classification:
        return Classification.create(
            user=self.user, lab=self.lab, lab_record_id="defaults", source=SubmissionSource.VARIANT_GRID,
            data={SpecialEKeys.ALLELE_ORIGIN: "somatic"}, **kwargs)

    def test_defaults_fill_the_mandatory_fields(self):
        classification = self._create(populate_with_defaults=True)
        self.assertEqual(classification.get(SpecialEKeys.CURATED_BY), "Dave Lawrence")
        self.assertEqual(classification.get(SpecialEKeys.CURATION_DATE), date.today().isoformat())

    def test_a_record_made_without_defaults_starts_empty(self):
        """ The API, file imports and sync create records without defaults """
        classification = self._create()
        self.assertIsNone(classification.get(SpecialEKeys.CURATED_BY))
        self.assertIsNone(classification.get(SpecialEKeys.CURATION_DATE))


class SomaticTierFromAmpLevelTestCase(TestCase):
    """ The tier a new somatic record starts at, from the AMP levels a copy brought with it """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def _classify(self, lab_record_id: str, **data) -> Classification:
        return Classification.create(
            user=self.user, lab=self.lab, lab_record_id=lab_record_id, source=SubmissionSource.VARIANT_GRID,
            data={SpecialEKeys.ALLELE_ORIGIN: "somatic", **data})

    def test_level_d_starts_at_tier_2(self):
        classification = self._classify("level_d", **{"amp:level_d": "therapeutic"})
        apply_somatic_tier_from_amp_level(classification, self.user)
        self.assertEqual(classification.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE), "tier_2")

    def test_the_highest_level_decides_the_tier(self):
        classification = self._classify("level_a_and_d", **{"amp:level_a": "therapeutic", "amp:level_d": "prognostic"})
        apply_somatic_tier_from_amp_level(classification, self.user)
        self.assertEqual(classification.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE), "tier_1")

    def test_a_tier_the_curator_already_has_is_kept(self):
        classification = self._classify("already_tiered", **{
            "amp:level_d": "therapeutic",
            SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE: "tier_3"
        })
        apply_somatic_tier_from_amp_level(classification, self.user)
        self.assertEqual(classification.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE), "tier_3")

    def test_a_record_with_no_levels_is_left_alone(self):
        classification = self._classify("no_levels")
        modifications = ClassificationModification.objects.filter(classification=classification).count()
        apply_somatic_tier_from_amp_level(classification, self.user)
        self.assertIsNone(classification.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE))
        self.assertEqual(ClassificationModification.objects.filter(classification=classification).count(),
                         modifications)
