"""
Tests for the representative Variant column and the Classifications column (variantgrid_private#686).

The renderers are client side, so what's tested here is what the server has to get right for them:
the members riding along hidden, the columns, the genomic sort and the expanded row.
"""
from django.test.client import Client
from django.urls.base import reverse
from django.utils.html import escape

from analysis.grids import VariantGrid
from analysis.models.nodes.sources.cohort_node import CohortNode
from analysis.tests.test_grid_export import GridExportTestCase
from annotation.models import ClinVar
from annotation.models.models_enums import ClinVarReviewStatus
from library.django_utils import FakeRequest
from snpdb.grid_columns.grid_sample_columns import get_available_format_columns
from snpdb.models import CompositeColumnMember, CustomColumnsCollection, UserSettings
from snpdb.views.datatable_view import datatable_definition, datatable_response


class RepresentativeVariantColumnTest(GridExportTestCase):
    def setUp(self):
        super().setUp()
        self.node = self._sample_node()
        self.grid = VariantGrid(FakeRequest(user=self.user), self.node)

    def _columns_by_name(self) -> dict:
        return {rc.name: rc for rc in self.grid.enabled_columns}

    def _order_by_column(self, name: str, direction: str = "asc"):
        index = next(i for i, rc in enumerate(self.grid.enabled_columns) if rc.name == name)
        self.grid.request.GET = {"order[0][column]": str(index), "order[0][dir]": direction}
        return self.grid.ordering(self.grid.get_initial_queryset())

    def test_variant_column_sorts_in_genome_build_order(self):
        qs = self._order_by_column("id")
        coords = [(v.locus.contig.name, v.locus.position) for v in qs]
        self.assertEqual(coords, [("1", 1000), ("1", 2000), ("2", 5000), ("2", 9000),
                                  ("3", 3000), ("10", 500), ("10", 1000)])

    def test_variant_column_sort_descending(self):
        qs = self._order_by_column("id", "desc")
        self.assertEqual([(v.locus.contig.name, v.locus.position) for v in qs][:2],
                         [("10", 1000), ("10", 500)])

    def _member_paths(self, composite_id: str) -> list[str]:
        return list(CompositeColumnMember.objects.filter(composite_id=composite_id)
                    .values_list("column__variant_column", flat=True))

    def test_renderer_fields_ride_along_hidden(self):
        """ A member rides along hidden while the collection shows the cell that draws it - unless
            it isn't annotated for this columns version, where there would be nothing to draw """
        columns = self._columns_by_name()
        self.assertEqual(columns["id"].client_renderer, "VariantGridFormat.representativeVariant")
        for composite_id, path in [("variant", "id"), ("classifications", "classifications"),
                                   ("gnomad", "gnomad"), ("spliceai", "spliceai")]:
            self.assertIn(path, columns, composite_id)
            drawn = [m for m in self._member_paths(composite_id) if m in columns]
            self.assertIn(self._member_paths(composite_id)[0], drawn, composite_id)  # the headline
            for member in drawn:
                self.assertEqual(member == path, columns[member].visible, member)
        for composite_id in ["annotsv_omim", "denovo_db"]:  # not in the default collection
            for member in self._member_paths(composite_id):
                self.assertNotIn(member, columns, member)
        # Not in the default collection, so hidden, and labelled from the catalogue for the CSV header
        self.assertFalse(columns["locus__ref__seq"].visible)
        self.assertEqual(columns["locus__ref__seq"].label, "Reference")
        self.assertEqual(columns["alt__seq"].label, "Alt")

    def test_rows_carry_renderer_fields(self):
        columns = self._columns_by_name()
        self.grid.request.GET = {"length": "1"}
        row = datatable_response(self.grid)["data"][0]
        for composite_id in ["variant", "classifications", "gnomad"]:
            for member in self._member_paths(composite_id):
                if member in columns:
                    self.assertIn(member, row, member)
        self.assertIsNone(row["classifications"])  # renderer-only column, like tags

    def test_composite_columns_carry_their_renderers(self):
        columns = self._columns_by_name()
        for name, client_renderer in [
                ("classifications", "VariantGridFormat.classifications"),
                ("consequence_impact", "VariantGridFormat.impactConsequence"),
                ("gnomad", "VariantGridFormat.gnomad"),
                ("spliceai", "VariantGridFormat.spliceai"),
                ("maxentscan", "VariantGridFormat.maxentscan"),
                ("mastermind", "VariantGridFormat.mastermind"),
                ("aloft", "VariantGridFormat.aloft"),
                ("predictions", "VariantGridFormat.predictions"),
                ("db_zygosity", "VariantGridFormat.dbZygosityCounts")]:
            self.assertEqual(columns[name].client_renderer, client_renderer, name)

    def test_composite_sort_menus_name_columns_in_the_grid(self):
        """ An entry naming a column this grid doesn't carry is dropped client side, so the menu has
            to name the composite's (hidden) members """
        columns = self._columns_by_name()
        menus = {name: rc.sort_menu for name, rc in columns.items() if rc.sort_menu}
        self.assertIn("consequence_impact", menus)
        for name, sort_menu in menus.items():
            # An entry for a member this columns version doesn't annotate is dropped client side
            offered = [entry["column"] for entry in sort_menu if entry["column"] in columns]
            self.assertTrue(offered, name)
            for column in offered:
                self.assertTrue(columns[column].orderable, column)
                # Picking the entry sorts on that column - which has to resolve
                list(self._order_by_column(column)[:1])

    def test_gnomad_members_keep_their_server_side_formatting(self):
        """ The gnomAD AFs are formatted server side (unit -> percent) so the CSV matches the grid -
            the client renderer only lays them out and adds the Pass/Fail link """
        columns = self._columns_by_name()
        for name in ["variantannotation__gnomad_popmax_af", "variantannotation__gnomad_af"]:
            rc = columns[name]
            self.assertFalse(rc.visible, "Drawn inside the gnomAD cell")
            self.assertIsNotNone(rc.renderer)
            self.assertTrue(rc.csv_rendered, "Server rendered value is what goes in the CSV")

    def test_sample_zygosity_cell_carries_its_partners(self):
        """ The whole call is drawn inside the zygosity cell, so the rest ride along hidden """
        columns = self._columns_by_name()
        sample_pk = self.sample.pk
        # A VCF without a GQ/PL/FT field has no such column to draw, sort on or ride along
        available = get_available_format_columns(self.grid.cohorts)
        partners = [c for c in ["allele_frequency", "allele_depth", "read_depth",
                                "genotype_quality", "phred_likelihood", "filters"]
                    if available[f"samples_{c}"]]
        zygosity = columns[f"sample_{sample_pk}_samples_zygosity"]
        self.assertEqual(zygosity.client_renderer, "VariantGridFormat.sampleZygosity")
        self.assertEqual(zygosity.client_renderer_kwargs, {"samplePrefix": f"sample_{sample_pk}_"})
        # One sort key per value the cell shows, each naming the column carrying that key's sort index
        self.assertEqual([entry["column"] for entry in zygosity.sort_menu],
                         [f"sample_{sample_pk}_samples_{c}" for c in ["zygosity"] + partners])
        for partner in partners:
            self.assertFalse(columns[f"sample_{sample_pk}_samples_{partner}"].visible, partner)

    def test_packed_genotype_columns_sort_on_their_own_value(self):
        """ A sample column's value is packed into the cohort's array/string, so sorting annotates
            this sample's value out of it """
        sample_pk = self.sample.pk
        for column in ["samples_zygosity", "samples_read_depth"]:
            qs = self._order_by_column(f"sample_{sample_pk}_{column}")
            self.assertIn(VariantGrid.GENOTYPE_SORT_ALIAS_PREFIX, str(qs.query), column)
            list(qs[:1])  # ...and it has to be a query the database will run

    def test_two_line_rows_is_a_user_setting(self):
        self.assertEqual(["variantgrid-datatable"], self.grid.get_table_classes())
        user_settings_override = UserSettings.get_settings_overrides(user=self.user)[-1]
        user_settings_override.variant_grid_two_line_rows = True
        user_settings_override.save()
        grid = VariantGrid(FakeRequest(user=self.user), self.node)
        self.assertEqual(["variantgrid-datatable", "two-line-rows"], grid.get_table_classes())

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
        for removed in ["chrom", "position", "ref", "alt", "svlen", "hgvs_g", "consequence", "impact",
                        "gnomad_af", "gnomad_popmax", "gnomad_filtered", "spliceai_max_ds",
                        "spliceai_pred_ds_ag", "maxentscan_ref", "mastermind_count_1_cdna",
                        "predictions_num_benign", "total_db_het", "total_db_hom", "aloft_pred"]:
            self.assertNotIn(removed, columns)
        for composite in ["consequence_impact", "gnomad", "spliceai", "maxentscan", "mastermind",
                          "predictions", "db_zygosity", "aloft"]:
            self.assertIn(composite, columns)


class CohortNodeCompositeColumnsTest(GridExportTestCase):
    """ The cohort node's own columns - the counts cell and the record level FILTER """

    def setUp(self):
        super().setUp()
        self.node = CohortNode.objects.create(analysis=self.analysis, cohort=self.cohort,
                                              accordion_panel=CohortNode.COUNT)
        self.grid = VariantGrid(FakeRequest(user=self.user), self.node)

    def _columns_by_name(self) -> dict:
        return {rc.name: rc for rc in self.grid.enabled_columns}

    def test_counts_cell_carries_its_partners(self):
        columns = self._columns_by_name()
        het = columns[self.node.het_count_column]
        self.assertEqual(het.client_renderer, "VariantGridFormat.dbZygosityCounts")
        self.assertEqual(het.client_renderer_kwargs, {"countPrefix": self.node.count_column_prefix})
        self.assertEqual([entry["column"] for entry in het.sort_menu],
                         [self.node.het_count_column, self.node.hom_count_column,
                          self.node.ref_count_column])
        for partner in [self.node.hom_count_column, self.node.ref_count_column]:
            self.assertFalse(columns[partner].visible, partner)

    def test_rows_carry_every_count_the_cell_draws(self):
        self.grid.request.GET = {"length": "1"}
        row = datatable_response(self.grid)["data"][0]
        for column in [self.node.het_count_column, self.node.hom_count_column, self.node.ref_count_column]:
            self.assertIn(column, row)

    def test_record_filters_column_is_drawn_client_side(self):
        cgc = self.node.cohort_genotype_collection
        filters = self._columns_by_name()[f"{cgc.cohortgenotype_alias}__filters"]
        self.assertEqual(filters.client_renderer, "VariantGridFormat.vcfFilters")
        # Still expanded to the VCF's own filter descriptions server side, so the CSV matches
        self.assertIsNotNone(filters.renderer)
        self.assertTrue(filters.csv_rendered)


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

    def test_row_detail_shows_clinvar(self):
        """ The Classifications chips are a summary - the expanded row is where the ClinVar detail is """
        variant = self.variants[0]
        annotation_version = self.analysis.annotation_version
        ClinVar.objects.create(version=annotation_version.clinvar_version, variant=variant,
                               clinvar_variation_id=12345, clinvar_allele_id=678,
                               preferred_disease_name="Cystic fibrosis",
                               review_status=ClinVarReviewStatus.CRITERIA_PROVIDED_SINGLE_SUBMITTER,
                               clinical_significance="Likely_pathogenic", highest_pathogenicity=4,
                               origin=1)
        client = Client()
        client.force_login(self.user)
        url = reverse("variant_grid_row_detail",
                      kwargs={"variant_id": variant.pk,
                              "annotation_version_id": annotation_version.pk})
        response = client.get(url)
        self.assertContains(response, "clinvar/variation/12345")
        self.assertContains(response, "Likely Pathogenic")
        self.assertContains(response, "Cystic fibrosis")
        self.assertContains(response, "germline")
