from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.analysis_templates import populate_analysis_from_template_run
from analysis.enums import AnalysisTemplateType
from analysis.forms.forms_nodes import PhenotypeNodeForm
from analysis.models import (
    Analysis,
    AnalysisTemplate,
    AnalysisTemplateRun,
    AnalysisVariable,
    CohortNode,
    PhenotypeNode,
    TrioNode,
)
from analysis.tests.utils import AnalysisSetupMixin, create_annotated_variants, node_keeps
from annotation.fake_data import create_fake_variant_annotation, get_fake_annotation_version
from annotation.models.models import VariantAnnotation
from genes.enums import AnnotationConsortium, HGNCStatus
from genes.fake_data import create_fake_transcript_version
from genes.models.models_gene import HGNC, Gene, HGNCImport, UniProt
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Patient
from patients.models.models_phenotype import (
    PhenotypeDescription,
    TextPhenotype,
    TextPhenotypeMatch,
    TextPhenotypeSentence,
)
from snpdb.fake_data import create_fake_cohort, create_fake_trio
from snpdb.models import Cohort, GenomeBuild
from snpdb.models.models_cached_web_resource import CachedWebResource
from snpdb.models.models_variant import Variant
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

# In the test ontology this HPO term reaches these genes through OMIM (ontology/tests/test_ontology.py)
COHORT_TERM = "HP:0001507"
COHORT_TERM_GENES = {"BCKDHA", "MEN1"}


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPhenotypeNodeText(AnalysisSetupMixin, TestCase):
    """ Free text is matched against gene descriptions as well as the ontology terms """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        vav = cls.analysis.annotation_version.variant_annotation_version
        transcript_version = create_fake_transcript_version(cls.grch37)
        uniprot = UniProt.objects.create(accession="TESTPHENO1",
                                         cached_web_resource=CachedWebResource.objects.create(name="test_phenotype"),
                                         function="Channel whose loss causes epilepsy.")
        gene_version = transcript_version.gene_version
        gene_version.hgnc = HGNC.objects.create(gene_symbol=gene_version.gene_symbol,
                                                hgnc_import=HGNCImport.objects.create(), status=HGNCStatus.APPROVED,
                                                approved_name="RUNX family transcription factor 1",
                                                alias_symbols="", uniprot=uniprot)
        gene_version.save()

        rows = {}
        for i, (name, gene_summary, variant_transcript_version) in enumerate([
                ("summary_epilepsy", "Mutations cause childhood epilepsy.", None),
                ("uniprot_epilepsy", None, transcript_version),
                ("summary_ataxia", "Mutations cause spinocerebellar ataxia.", None),
                ("unrelated", "Required for kidney development.", None)]):
            gene = None
            if gene_summary:
                gene = Gene.objects.create(identifier=f"TEST_PHENOTYPE_{i}",
                                           annotation_consortium=AnnotationConsortium.REFSEQ, summary=gene_summary)
            rows[name] = {"gene": gene, "transcript_version": variant_transcript_version}
        cls.variants = create_annotated_variants(vav, rows)

    def _node(self, text_phenotype=None) -> PhenotypeNode:
        return PhenotypeNode.objects.create(analysis=self.analysis, text_phenotype=text_phenotype)

    def _filtered(self, node: PhenotypeNode) -> set[str]:
        return node_keeps(node, self.variants)

    def test_no_terms_or_text_does_not_modify_parents(self):
        self.assertFalse(self._node().modifies_parents())

    def test_text_modifies_parents(self):
        self.assertTrue(self._node("epilepsy").modifies_parents())

    def test_text_searches_every_description_column(self):
        self.assertEqual({"summary_epilepsy", "uniprot_epilepsy"}, self._filtered(self._node("epilepsy")))

    def test_each_word_is_searched_separately(self):
        self.assertEqual({"summary_epilepsy", "uniprot_epilepsy", "summary_ataxia"},
                         self._filtered(self._node("epilepsy ataxia")))

    def test_node_name_shows_the_text(self):
        self.assertIn("epilepsy", self._node("epilepsy").get_node_name())


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPhenotypeNodeSource(AnalysisSetupMixin, TestCase):
    """ The node's terms come from an ancestor cohort or patient - @see PhenotypeNode.get_phenotype_source """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = cls.analysis.user
        cls.cohort = create_fake_cohort(cls.user, cls.grch37)
        # Skip the NLP matching - these tests patch get_ontology_term_ids on the source
        cls.cohort.phenotype = "Cohort condition"
        cls.cohort.save(check_patient_text_phenotype=False)

    def _cohort_node(self, cohort=None) -> CohortNode:
        return CohortNode.objects.create(analysis=self.analysis, cohort=cohort or self.cohort)

    def _phenotype_node(self, parent, **kwargs) -> PhenotypeNode:
        node = PhenotypeNode.objects.create(analysis=self.analysis, **kwargs)
        node.add_parent(parent)
        node._cached_parents = None
        node.save()
        return PhenotypeNode.objects.get(pk=node.pk)

    def _trio_node_with_proband_patient(self) -> tuple[TrioNode, Patient]:
        """ A trio names its proband sample, so the node has a proband patient to prefer """
        trio = create_fake_trio(self.user, self.grch37)
        trio.cohort.phenotype = "Trio cohort condition"
        trio.cohort.save(check_patient_text_phenotype=False)

        patient = Patient.objects.create(patient_code="PHENO_PAT")
        proband_sample = trio.proband.sample
        proband_sample.patient = patient
        proband_sample.save()
        return TrioNode.objects.create(analysis=self.analysis, trio=trio), patient

    # ── The source decides whose terms are read ──────────────────────────────

    def test_terms_and_genes_come_from_the_cohort(self):
        description = PhenotypeDescription.objects.create(cohort=self.cohort, original_text=self.cohort.phenotype)
        text_phenotype = TextPhenotype.objects.create(text=self.cohort.phenotype)
        TextPhenotypeSentence.objects.create(phenotype_description=description, text_phenotype=text_phenotype,
                                             sentence_offset=0)
        TextPhenotypeMatch.objects.create(text_phenotype=text_phenotype, ontology_term_id=COHORT_TERM,
                                          offset_start=0, offset_end=len(self.cohort.phenotype))

        node = self._phenotype_node(self._cohort_node(), accordion_panel=PhenotypeNode.PANEL_SOURCE)
        self.assertEqual(node.cohort, self.cohort)
        self.assertEqual([COHORT_TERM], list(node.get_ontology_term_ids()))
        self.assertEqual(COHORT_TERM_GENES, set(node.get_gene_symbols_qs().values_list("symbol", flat=True)))

    def test_node_name_names_the_cohort(self):
        node = self._phenotype_node(self._cohort_node(), accordion_panel=PhenotypeNode.PANEL_SOURCE)
        node.cohort = self.cohort
        with patch.object(Cohort, "get_ontology_term_ids", return_value=[COHORT_TERM]), \
                patch.object(PhenotypeNode, "modifies_parents", return_value=True):
            self.assertEqual(node.get_node_name(), f"{self.cohort.name} cohort phenotypes")

    # ── Auto-configuration ───────────────────────────────────────────────────

    def test_a_lone_ancestor_cohort_with_phenotype_text_is_auto_selected(self):
        node = self._phenotype_node(self._cohort_node())
        self.assertEqual(node.cohort, self.cohort)
        self.assertEqual(node.accordion_panel, PhenotypeNode.PANEL_SOURCE)

    def test_a_cohort_without_phenotype_text_is_left_unset(self):
        Cohort.objects.filter(pk=self.cohort.pk).update(phenotype=None)
        node = self._phenotype_node(self._cohort_node(Cohort.objects.get(pk=self.cohort.pk)))
        self.assertIsNone(node.cohort)
        self.assertIsNone(node.patient)

    def test_a_proband_patient_wins_over_a_described_cohort(self):
        trio_node, patient = self._trio_node_with_proband_patient()
        node = self._phenotype_node(trio_node)
        self.assertEqual(node.patient, patient)
        self.assertIsNone(node.cohort)

    def test_a_cohort_no_longer_an_ancestor_is_cleared_and_re_picked(self):
        other_cohort = create_fake_cohort(self.user, self.grch37)
        other_cohort.phenotype = "Other condition"
        other_cohort.save(check_patient_text_phenotype=False)

        node = self._phenotype_node(self._cohort_node(other_cohort))
        self.assertEqual(node.cohort, other_cohort)

        node.remove_parent(node.get_parent_subclasses()[0])
        node.add_parent(self._cohort_node())
        node._cached_parents = None
        node.save()

        node = PhenotypeNode.objects.get(pk=node.pk)
        self.assertEqual(node.cohort, self.cohort)

    # ── Configuration errors ─────────────────────────────────────────────────

    def test_a_cohort_outside_the_ancestors_is_a_configuration_error(self):
        other_cohort = create_fake_cohort(self.user, self.grch37)
        node = self._phenotype_node(self._cohort_node())
        node._set_cohort(other_cohort)
        errors = node._get_configuration_errors()
        self.assertTrue(any("is not a cohort in any ancestors" in str(e) for e in errors), errors)

    # ── The picker ───────────────────────────────────────────────────────────

    def _visible_patient(self) -> Patient:
        patient = Patient.objects.create(patient_code="PICKER_PAT")
        assign_permission_to_user_and_groups(self.user, patient)
        sample = self.cohort.get_samples()[0]
        sample.patient = patient
        sample.save()
        return patient

    def _form(self, node, data=None) -> PhenotypeNodeForm:
        return PhenotypeNodeForm(data, instance=PhenotypeNode.objects.get(pk=node.pk))

    def test_the_picker_offers_ancestor_patients_and_cohorts(self):
        patient = self._visible_patient()
        node = self._phenotype_node(self._cohort_node())
        grouped = dict(self._form(node).fields["phenotype_source"].choices)

        self.assertIn((f"patient:{patient.pk}", str(patient)), grouped["Patients"])
        self.assertIn((f"cohort:{self.cohort.pk}", self.cohort.name), grouped["Cohorts"])

    def test_a_cohort_without_phenotype_text_is_offered_but_flagged(self):
        Cohort.objects.filter(pk=self.cohort.pk).update(phenotype=None)
        node = self._phenotype_node(self._cohort_node(Cohort.objects.get(pk=self.cohort.pk)))
        grouped = dict(self._form(node).fields["phenotype_source"].choices)

        self.assertIn((f"cohort:{self.cohort.pk}", f"{self.cohort.name} (no phenotype)"), grouped["Cohorts"])

    def test_saving_a_cohort_clears_the_patient(self):
        patient = self._visible_patient()
        node = self._phenotype_node(self._cohort_node())
        node._set_patient(patient)
        node.save()

        form = self._form(node, {"phenotype_source": f"cohort:{self.cohort.pk}",
                                 "accordion_panel": PhenotypeNode.PANEL_SOURCE})
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()

        self.assertEqual(node.cohort, self.cohort)
        self.assertIsNone(node.patient)
        self.assertEqual(self._form(node).get_analysis_variable_field("phenotype_source"), "cohort")

    def test_saving_a_patient_clears_the_cohort(self):
        patient = self._visible_patient()
        node = self._phenotype_node(self._cohort_node())
        node._set_cohort(self.cohort)
        node.save()

        form = self._form(node, {"phenotype_source": f"patient:{patient.pk}",
                                 "accordion_panel": PhenotypeNode.PANEL_SOURCE})
        self.assertTrue(form.is_valid(), form.errors)
        node = form.save()

        self.assertEqual(node.patient, patient)
        self.assertIsNone(node.cohort)
        self.assertEqual(self._form(node).get_analysis_variable_field("phenotype_source"), "patient")


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPhenotypeNodeTemplateRun(TestCase):
    """ A template's phenotype node picks up the run's cohort, as the sample pattern does """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_phenotype_template")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.cohort = create_fake_cohort(cls.user, cls.grch37)
        cls.cohort.phenotype = "Cohort condition"
        cls.cohort.save(check_patient_text_phenotype=False)

    def test_the_run_cohort_is_set_on_the_phenotype_node(self):
        analysis = Analysis(genome_build=self.grch37, template_type=AnalysisTemplateType.TEMPLATE)
        analysis.set_defaults_and_save(self.user)
        cohort_node = CohortNode.objects.create(analysis=analysis)
        AnalysisVariable.objects.create(node=cohort_node, field="cohort", class_name="snpdb.Cohort")
        phenotype_node = PhenotypeNode.objects.create(analysis=analysis,
                                                      accordion_panel=PhenotypeNode.PANEL_SOURCE)
        phenotype_node.add_parent(cohort_node)
        phenotype_node._cached_parents = None
        phenotype_node.save()
        AnalysisVariable.objects.create(node=phenotype_node, field="cohort", class_name="snpdb.Cohort")

        template = AnalysisTemplate.objects.create(name="phenotype template", user=self.user, analysis=analysis)
        version = template.new_version("%(template)s for %(input)s")
        template_run = AnalysisTemplateRun.create(template, self.grch37, user=self.user, template_version=version)
        template_run.populate_arguments({"cohort": self.cohort})
        populate_analysis_from_template_run(template_run)

        run_node = PhenotypeNode.objects.get(analysis=template_run.analysis)
        self.assertEqual(run_node.cohort, self.cohort)
