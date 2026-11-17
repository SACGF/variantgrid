from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone

from genes.enums import PanelAppConfidence
from genes.models import (
    ActiveSampleGeneList,
    GeneCoverageCanonicalTranscript,
    GeneCoverageCollection,
    GeneList,
    GeneListCategory,
    GeneListGeneSymbol,
    GeneSymbol,
    PanelAppPanel,
    PanelAppPanelLocalCache,
    PanelAppPanelLocalCacheGeneSymbol,
    PanelAppServer,
    SampleGeneList,
)
from genes.panel_app import create_gene_list_from_panel_app
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.enums import DataState, ImportStatus
from snpdb.models import VCF, Sample
from snpdb.models.models_genome import GenomeBuild


def _create_gene_list(user, symbols, **kwargs) -> GeneList:
    gene_list = GeneList.objects.create(name="goi", user=user, import_status=ImportStatus.SUCCESS, **kwargs)
    GeneListGeneSymbol.objects.bulk_create([GeneListGeneSymbol(gene_list=gene_list, original_name=symbol,
                                                               gene_symbol_id=symbol) for symbol in symbols])
    return gene_list


class GenesOfInterestTest(TestCase):
    SYMBOLS = ["GOIA", "GOIB", "GOIC"]

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="goi_user")
        cls.genome_build = GenomeBuild.grch37()
        for symbol in cls.SYMBOLS + ["GOID"]:
            GeneSymbol.objects.get_or_create(symbol=symbol)

        vcf = VCF.objects.create(name="goi_vcf", genotype_samples=1, genome_build=cls.genome_build,
                                 import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())
        assign_permission_to_user_and_groups(cls.user, vcf)
        cls.sample = Sample.objects.create(name="SMP5008", vcf=vcf, import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.user, cls.sample)

        server = PanelAppServer.objects.create(name="Test PanelApp GOI", url="https://panelapp.example.com/goi",
                                               icon_css_class="")
        cls.panel = PanelAppPanel.objects.create(server=server, panel_id=1, disease_group="", disease_sub_group="",
                                                 name="GOI Panel", status="public", current_version="1.0")
        cls.local_cache = PanelAppPanelLocalCache.objects.create(panel_app_panel=cls.panel, version="1.0")
        for symbol, confidence in [("GOIA", "3"), ("GOIB", "3"), ("GOIC", "1")]:
            PanelAppPanelLocalCacheGeneSymbol.objects.create(panel_app_local_cache=cls.local_cache,
                                                             gene_symbol_reported=symbol,
                                                             data={"confidence_level": confidence})

    def _active_sample_gene_list(self):
        return ActiveSampleGeneList.objects.get(sample=self.sample).sample_gene_list

    def test_new_sample_gene_list_becomes_active_unless_pipeline_replaces_user(self):
        pipeline_1 = SampleGeneList.objects.create(sample=self.sample, gene_list=_create_gene_list(self.user, ["GOIA"]))
        active_pk = ActiveSampleGeneList.objects.get(sample=self.sample).pk
        self.assertEqual(pipeline_1, self._active_sample_gene_list())

        pipeline_2 = SampleGeneList.objects.create(sample=self.sample, gene_list=_create_gene_list(self.user, ["GOIB"]))
        self.assertEqual(pipeline_2, self._active_sample_gene_list())

        user_applied = SampleGeneList.objects.create(sample=self.sample, user=self.user,
                                                     gene_list=_create_gene_list(self.user, ["GOIC"]))
        self.assertEqual(user_applied, self._active_sample_gene_list())

        SampleGeneList.objects.create(sample=self.sample, gene_list=_create_gene_list(self.user, ["GOID"]))
        self.assertEqual(user_applied, self._active_sample_gene_list())
        # Switched in place, so the 'created' signal that auto-runs templates only fired for the first
        self.assertEqual(active_pk, ActiveSampleGeneList.objects.get(sample=self.sample).pk)

    @patch("genes.panel_app.get_panel_app_local_cache")
    def test_panel_app_provenance(self, mock_get_local_cache):
        mock_get_local_cache.return_value = self.local_cache
        gene_list = create_gene_list_from_panel_app(self.panel, PanelAppConfidence.INTERMEDIATE, self.user)
        self.assertEqual({"GOIA", "GOIB"}, set(gene_list.get_gene_names()))

        provenance = gene_list.get_panel_app_provenance()
        self.assertFalse(provenance.is_edited)
        self.assertTrue(provenance.label.startswith("PanelApp: GOI Panel v1.0 (Intermediate), created by goi_user"))

        gene_list.add_and_remove_gene_symbols(["GOID"], ["GOIA"], {"GOID": "Added by goi_user"})
        copy = gene_list.copy_for_user(self.user)
        provenance = copy.get_panel_app_provenance()
        self.assertEqual("PanelApp: GOI Panel v1.0 (Intermediate) + manual changes: 1 added, 1 removed",
                         provenance.label)
        self.assertEqual([("GOID", "Added by goi_user")],
                         [(a.name, a.modification_info) for a in provenance.added])
        self.assertEqual(["GOIA"], provenance.removed)

        self.assertFalse(provenance.has_newer_version)
        PanelAppPanel.objects.filter(pk=self.panel.pk).update(current_version="1.1")
        self.assertTrue(GeneList.objects.get(pk=copy.pk).get_panel_app_provenance().has_newer_version)

        self.client.force_login(self.user)
        response = self.client.get(reverse("view_gene_list", kwargs={"gene_list_id": copy.pk}))
        self.assertContains(response, "PanelApp v1.0 (pinned)")
        self.assertContains(response, "Create from latest version")

    def test_coverage_splits_below_threshold_from_no_data(self):
        gcc = GeneCoverageCollection.objects.create(path="/tmp/goi.cov.tsv", data_state=DataState.COMPLETE,
                                                    genome_build=self.genome_build)
        for symbol, min_depth, mean, size in [("GOIA", 30, 100.0, 100), ("GOIB", 5, 40.0, 300)]:
            GeneCoverageCanonicalTranscript.objects.create(gene_coverage_collection=gcc, gene_symbol_id=symbol,
                                                           original_gene_symbol=symbol, original_transcript="",
                                                           size=size, min=min_depth, mean=mean, std_dev=1.0)
        gene_list = _create_gene_list(self.user, self.SYMBOLS)
        GeneListGeneSymbol.objects.create(gene_list=gene_list, original_name="NOTAGENE")

        gene_symbols = GeneSymbol.objects.filter(pk__in=self.SYMBOLS)
        self.assertEqual(["GOIB"], list(gcc.get_below_threshold_gene_symbols(gene_symbols, 20)
                                        .values_list("symbol", flat=True)))
        self.assertEqual(["GOIC"], list(gcc.get_no_coverage_gene_symbols(gene_symbols)
                                        .values_list("symbol", flat=True)))

        summary = gcc.get_gene_lists_coverage_summary([gene_list], 20)
        self.assertEqual((4, 2, 5, 1), (summary.num_genes, summary.num_genes_with_coverage,
                                        summary.lowest_min, summary.num_below_threshold))
        self.assertAlmostEqual((100 * 100 + 40 * 300) / 400, summary.mean_depth)
        self.assertEqual((["GOIC"], ["NOTAGENE"]), (summary.no_coverage_symbols, summary.unmatched_names))

    def test_apply_existing_gene_list_to_sample(self):
        source = _create_gene_list(self.user, ["GOIA", "GOIB"])
        self.client.force_login(self.user)
        url = reverse("sample_gene_lists_tab", kwargs={"sample_id": self.sample.pk})
        response = self.client.post(url, {"action": "apply", "apply-gene_list": source.pk})
        self.assertEqual(200, response.status_code)

        sample_gene_list = self._active_sample_gene_list()
        self.assertEqual(self.user, sample_gene_list.user)
        applied = sample_gene_list.gene_list
        self.assertNotEqual(source, applied)  # a copy, so editing the source doesn't change the GOI
        self.assertEqual(("goi", GeneListCategory.SAMPLE_GENE_LIST), (applied.name, applied.category.name))
        self.assertEqual({"GOIA", "GOIB"}, set(applied.get_gene_names()))
        self.assertContains(response, "Applied by goi_user")
