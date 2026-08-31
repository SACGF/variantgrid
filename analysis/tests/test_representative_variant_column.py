"""
Tests for the representative Variant column and the Classifications column (variantgrid_private#686).

The renderers are client side, so what's tested here is what the server has to get right for them:
the hidden row fields riding along, the colmodels, the genomic sort and the expanded row.
"""
from django.test.client import Client
from django.urls.base import reverse
from django.utils.html import escape

from analysis.grids import VariantGrid
from analysis.models.nodes.sources.cohort_node import CohortNode
from analysis.tests.test_grid_export import GridExportTestCase
from library.django_utils import FakeRequest
from library.django_utils.jqgrid_datatable_adapter import datatable_definition
from snpdb.grid_columns.custom_columns import (
    CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS,
    CLASSIFICATIONS_COLUMN_ROW_FIELDS,
    COMPOSITE_COLUMN_ROW_FIELDS,
    VARIANT_COLUMN_ROW_FIELDS,
)
from snpdb.grid_columns.grid_sample_columns import get_available_format_columns
from snpdb.models import CustomColumnsCollection, UserSettings


class RepresentativeVariantColumnTest(GridExportTestCase):
    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.grid = VariantGrid(self.user, self.node)

    def _colmodels_by_name(self) -> dict:
        return {cm["name"]: cm for cm in self.grid.get_colmodels()}

    def test_variant_column_sorts_in_genome_build_order(self):
        qs = self.grid._sort_items(self.node.get_queryset(), "id", "asc")
        coords = [(v.locus.contig.name, v.locus.position) for v in qs]
        self.assertEqual(coords, [("1", 1000), ("1", 2000), ("2", 5000), ("2", 9000),
                                  ("3", 3000), ("10", 500), ("10", 1000)])

    def test_variant_column_sort_descending(self):
        qs = self.grid._sort_items(self.node.get_queryset(), "id", "desc")
        self.assertEqual([(v.locus.contig.name, v.locus.position) for v in qs][:2],
                         [("10", 1000), ("10", 500)])

    def test_renderer_fields_ride_along_hidden(self):
        colmodels = self._colmodels_by_name()
        self.assertEqual(colmodels["id"]["formatter"], "representativeVariant")
        for field in (VARIANT_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_FIELDS
                      + CLASSIFICATIONS_COLUMN_ROW_ANNOTATIONS):
            self.assertIn(field, colmodels)
        for visible, partners in COMPOSITE_COLUMN_ROW_FIELDS.items():
            # A partner only rides along while the collection shows the column that reads it
            for partner in partners:
                self.assertEqual(visible in colmodels, partner in colmodels, partner)
        # Not in the default collection, so hidden, and labelled from the catalogue for the CSV header
        self.assertTrue(colmodels["locus__ref__seq"]["hidden"])
        self.assertEqual(colmodels["locus__ref__seq"]["label"], "Reference")
        self.assertEqual(colmodels["alt__seq"]["label"], "Alt")

    def test_rows_carry_renderer_fields(self):
        request = FakeRequest(user=self.user)
        request.GET = {"rows": "1", "page": "1"}
        row = self.grid.get_data(request)["rows"][0]
        composite_partners = [p for visible, partners in COMPOSITE_COLUMN_ROW_FIELDS.items()
                              for p in partners if visible in row]
        for field in VARIANT_COLUMN_ROW_FIELDS + CLASSIFICATIONS_COLUMN_ROW_FIELDS + composite_partners:
            self.assertIn(field, row)
        self.assertNotIn("classifications", row)  # renderer-only column, like tags

    def test_classifications_column_is_not_sortable(self):
        cm = self._colmodels_by_name()["classifications"]
        self.assertEqual(cm["formatter"], "classificationsFormatter")
        self.assertFalse(cm["sortable"])

    def test_composite_columns_carry_their_renderers(self):
        colmodels = self._colmodels_by_name()
        for name, formatter in [
                ("variantannotation__consequence", "impactConsequenceFormatter"),
                ("variantannotation__spliceai_max_ds", "spliceaiFormatter"),
                ("variantannotation__maxentscan_percent_diff_ref", "maxentscanFormatter"),
                ("variantannotation__mastermind_count_1_cdna", "mastermindFormatter"),
                ("variantannotation__aloft_pred", "aloftFormatter"),
                ("variantannotation__predictions_num_pathogenic", "predictionsFormatter"),
                ("global_variant_zygosity__het_count", "dbZygosityCountsFormatter")]:
            self.assertEqual(colmodels[name]["formatter"], formatter, name)

    def test_gnomad_columns_keep_their_server_side_formatting(self):
        """ The gnomAD AFs are formatted server side (unit -> percent) so the CSV matches the grid -
            the client renderers only add the population / the Pass-Fail link beside them """
        colmodels = self._colmodels_by_name()
        for name, formatter in [("variantannotation__gnomad_popmax_af", "gnomadPopmaxFormatter"),
                                ("variantannotation__gnomad_af", "gnomadAfFormatter")]:
            cm = colmodels[name]
            self.assertEqual(cm["formatter"], formatter)
            self.assertIn("server_side_formatter", cm)

    def test_sample_zygosity_cell_carries_its_partners(self):
        """ The whole call is drawn inside the zygosity cell, so the rest ride along hidden """
        colmodels = self._colmodels_by_name()
        sample_pk = self.sample.pk
        # A VCF without a GQ/PL/FT field has no such column to draw, sort on or ride along
        available = get_available_format_columns(self.grid.cohorts)
        partners = [c for c in ["allele_frequency", "allele_depth", "read_depth",
                                "genotype_quality", "phred_likelihood", "filters"]
                    if available[f"samples_{c}"]]
        zygosity = colmodels[f"sample_{sample_pk}_samples_zygosity"]
        self.assertEqual(zygosity["formatter"], "sampleZygosityFormatter")
        self.assertEqual(zygosity["formatter_kwargs"], {"samplePrefix": f"sample_{sample_pk}_"})
        # One sort key per value the cell shows, each naming the column carrying that key's sort index
        self.assertEqual([entry["column"] for entry in zygosity["sort_menu"]],
                         [f"sample_{sample_pk}_samples_{c}" for c in ["zygosity"] + partners])
        for partner in partners:
            self.assertTrue(colmodels[f"sample_{sample_pk}_samples_{partner}"]["hidden"], partner)

    def test_two_line_rows_is_a_user_setting(self):
        self.assertEqual([], self.grid.get_extra_table_classes())
        user_settings_override = UserSettings.get_settings_overrides(user=self.user)[-1]
        user_settings_override.variant_grid_two_line_rows = True
        user_settings_override.save()
        self.assertEqual(["two-line-rows"], VariantGrid(self.user, self.node).get_extra_table_classes())

    def test_definition_declares_row_expansion(self):
        definition = datatable_definition(self.grid)
        annotation_version_id = self.node.analysis.annotation_version_id
        self.assertEqual(definition["expandClientRenderer"],
                         f"variantGridRowDetail.bind(null, {annotation_version_id})")
        self.assertFalse(definition["expandPrefetch"])

    def test_system_default_collection(self):
        columns = list(CustomColumnsCollection.get_system_default().customcolumn_set
                       .order_by("sort_order").values_list("column_id", flat=True))
        self.assertEqual(columns[:4], ["variant", "classifications", "tags", "tags_global"])
        # Coordinates now live in the Variant cell; the rest inside their composite cells
        for removed in ["chrom", "position", "ref", "alt", "svlen", "hgvs_g", "impact", "gnomad_popmax",
                        "gnomad_filtered", "spliceai_pred_ds_ag", "spliceai_pred_dp_dl", "maxentscan_ref",
                        "maxentscan_alt", "maxentscan_diff", "mastermind_mmid3", "predictions_num_benign",
                        "total_db_hom", "total_db_ref", "total_db_unk", "aloft_prob_dominant",
                        "aloft_high_confidence", "aloft_ensembl_transcript"]:
            self.assertNotIn(removed, columns)
        for merged in ["spliceai_max_ds", "maxentscan_percent_diff_ref", "mastermind_count_1_cdna",
                       "predictions_num_pathogenic", "total_db_het", "aloft_pred"]:
            self.assertIn(merged, columns)


class CohortNodeCompositeColumnsTest(GridExportTestCase):
    """ The cohort node's own columns - the counts cell and the record level FILTER """

    def setUp(self):
        super().setUp()
        self.node = CohortNode.objects.create(analysis=self.analysis, cohort=self.cohort,
                                              accordion_panel=CohortNode.COUNT)
        self.grid = VariantGrid(self.user, self.node)

    def _colmodels_by_name(self) -> dict:
        return {cm["name"]: cm for cm in self.grid.get_colmodels()}

    def test_counts_cell_carries_its_partners(self):
        colmodels = self._colmodels_by_name()
        het = colmodels[self.node.het_count_column]
        self.assertEqual(het["formatter"], "dbZygosityCountsFormatter")
        self.assertEqual(het["formatter_kwargs"], {"countPrefix": self.node.count_column_prefix})
        self.assertEqual([entry["column"] for entry in het["sort_menu"]],
                         [self.node.het_count_column, self.node.hom_count_column,
                          self.node.ref_count_column])
        for partner in [self.node.hom_count_column, self.node.ref_count_column]:
            self.assertTrue(colmodels[partner]["hidden"], partner)

    def test_rows_carry_every_count_the_cell_draws(self):
        request = FakeRequest(user=self.user)
        request.GET = {"rows": "1", "page": "1"}
        row = self.grid.get_data(request)["rows"][0]
        for column in [self.node.het_count_column, self.node.hom_count_column, self.node.ref_count_column]:
            self.assertIn(column, row)

    def test_record_filters_column_is_drawn_client_side(self):
        cgc = self.node.cohort_genotype_collection
        filters = self._colmodels_by_name()[f"{cgc.cohortgenotype_alias}__filters"]
        self.assertEqual(filters["formatter"], "vcfFiltersFormatter")
        # Still expanded to the VCF's own filter descriptions server side, so the CSV matches
        self.assertIn("server_side_formatter", filters)


class VariantGridRowDetailViewTest(GridExportTestCase):
    def test_row_detail_renders(self):
        client = Client()
        client.force_login(self.user)
        variant = self.variants[0]
        url = reverse("variant_grid_row_detail",
                      kwargs={"variant_id": variant.pk,
                              "annotation_version_id": self.analysis.annotation_version_id})
        response = client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, f"v {variant.pk}")
        self.assertContains(response, escape(variant.full_string))  # HGVS/coordinate '>' is escaped
