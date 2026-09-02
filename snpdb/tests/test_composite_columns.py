"""
Composite grid columns - the cell, the members riding along hidden, and the migration helper that
collapses a collection onto them.

@see snpdb.models.models_columns.CompositeColumnMember, snpdb/migrations/0238_composite_columns.py
"""
from django.apps import apps
from django.contrib.auth.models import User
from django.test import TestCase
from django.test.client import Client
from django.urls.base import reverse

from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils.composite_columns import collapse_into_composite
from snpdb.grid_columns.custom_columns import get_variant_grid_columns
from snpdb.models import CompositeColumnMember, CustomColumn, CustomColumnsCollection, VariantGridColumn
from snpdb.models.models_genome import GenomeBuild

SPLICEAI_MEMBERS = ["spliceai_max_ds", "spliceai_pred_ds_ag", "spliceai_pred_ds_al",
                    "spliceai_pred_ds_dg", "spliceai_pred_ds_dl", "spliceai_pred_dp_ag",
                    "spliceai_pred_dp_al", "spliceai_pred_dp_dg", "spliceai_pred_dp_dl",
                    "spliceai_gene_symbol"]


class CompositeColumnGridTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_composite_columns")[0]
        cls.annotation_version = get_fake_annotation_version(GenomeBuild.get_name_or_alias("GRCh38"))

    def _collection(self, column_ids: list[str]) -> CustomColumnsCollection:
        ccc = CustomColumnsCollection.objects.create(name=f"test {column_ids}", user=self.user)
        ccc.customcolumn_set.all().delete()  # incl. the mandatory Variant column
        for sort_order, column_id in enumerate(column_ids):
            CustomColumn.objects.create(custom_columns_collection=ccc, column_id=column_id,
                                        sort_order=sort_order)
        return ccc

    def _columns(self, column_ids: list[str]) -> list:
        rich_columns, _ = get_variant_grid_columns(self._collection(column_ids),
                                                   self.annotation_version, {})
        return rich_columns

    def test_composite_draws_its_members_hidden(self):
        columns = self._columns(["spliceai"])
        composite = columns[0]
        self.assertEqual("spliceai", composite.name)
        self.assertIsNone(composite.key, "Display only - it has no value of its own")
        self.assertEqual(["variantannotation__spliceai_max_ds"], composite.sort_keys)
        self.assertTrue(composite.orderable)

        members = VariantGridColumn.objects.filter(pk__in=SPLICEAI_MEMBERS)
        member_paths = {vgc.pk: vgc.variant_column for vgc in members}
        self.assertEqual([member_paths[pk] for pk in SPLICEAI_MEMBERS], [rc.name for rc in columns[1:]])
        for rc in columns[1:]:
            self.assertFalse(rc.visible, rc.name)
            self.assertTrue(rc.include_in_csv, rc.name)

    def test_sort_menu_and_render_kwargs_come_from_the_members(self):
        composite = self._columns(["spliceai"])[0]
        in_menu = list(CompositeColumnMember.objects.filter(composite_id="spliceai", in_sort_menu=True)
                       .values_list("column__variant_column", flat=True))
        self.assertEqual(in_menu, [entry["column"] for entry in composite.sort_menu])
        self.assertTrue(in_menu, "SpliceAI offers its four delta scores")
        self.assertEqual([VariantGridColumn.objects.get(pk=pk).variant_column for pk in SPLICEAI_MEMBERS],
                         [m["path"] for m in composite.client_renderer_kwargs["members"]])

    def test_generic_renderer_unless_an_override_names_one(self):
        self.assertEqual("VariantGridFormat.composite", self._columns(["cadd"])[0].client_renderer)

    def test_member_renderer_travels_with_the_member_entry(self):
        overrides = {"variantannotation__cosmic_id": {"client_renderer": "VariantGridFormat.cosmicLink"}}
        rich_columns, _ = get_variant_grid_columns(self._collection(["cosmic"]),
                                                   self.annotation_version, overrides)
        members = rich_columns[0].client_renderer_kwargs["members"]
        self.assertEqual("VariantGridFormat.cosmicLink", members[0]["renderer"])
        self.assertNotIn("renderer", members[1])

    def test_keyed_composite_keeps_its_key(self):
        columns = self._columns(["variant"])
        self.assertEqual("id", columns[0].key)
        self.assertIsNone(columns[0].sort_keys)
        self.assertIn("locus__contig__name", [rc.name for rc in columns[1:]])

    def test_member_shown_standalone_is_emitted_once_visible(self):
        columns = self._columns(["spliceai_pred_ds_ag", "spliceai"])
        names = [rc.name for rc in columns]
        self.assertEqual(1, names.count("variantannotation__spliceai_pred_ds_ag"))
        self.assertEqual(0, names.index("variantannotation__spliceai_pred_ds_ag"))
        self.assertTrue(columns[0].visible)

    def test_every_composite_has_members_that_are_real_columns(self):
        """ Guards the migration's table - a composite with no members would draw an empty cell """
        catalogue = set(VariantGridColumn.objects.values_list("pk", flat=True))
        # Display only and not drawn inside anything else - it has to be a composite or it draws nothing
        composites = VariantGridColumn.objects.filter(queryset_field=False, composite_membership__isnull=True) \
                                              .exclude(pk__in=["tags", "tags_global", "Sample"])
        for composite in composites.prefetch_related("composite_members__column"):
            self.assertTrue(composite.is_composite, composite.pk)
            for member in composite.member_columns:
                self.assertIn(member.pk, catalogue, member.pk)


class CollapseIntoCompositeTest(TestCase):
    """ The migration helper - every collection, user and system, comes through it """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_collapse_into_composite")[0]

    def _collection(self, column_ids: list[str]) -> CustomColumnsCollection:
        ccc = CustomColumnsCollection.objects.create(name=f"test {column_ids}", user=self.user)
        ccc.customcolumn_set.all().delete()
        for sort_order, column_id in enumerate(column_ids):
            CustomColumn.objects.create(custom_columns_collection=ccc, column_id=column_id,
                                        sort_order=sort_order)
        ccc.refresh_from_db()
        return ccc

    def _collapse(self, ccc, composite_id="spliceai") -> list[str]:
        collapse_into_composite(apps, composite_id)
        ccc.refresh_from_db()
        return list(ccc.customcolumn_set.order_by("sort_order").values_list("column_id", flat=True))

    def test_run_of_members_becomes_the_composite_at_the_runs_position(self):
        ccc = self._collection(["gene_symbol", "spliceai_max_ds", "spliceai_pred_ds_ag",
                                "spliceai_pred_dp_ag", "hgvs_c"])
        version = ccc.version_id
        self.assertEqual(["gene_symbol", "spliceai", "hgvs_c"], self._collapse(ccc))
        self.assertGreater(ccc.version_id, version, "The node grid definition cache is keyed on it")

    def test_lone_member_becomes_the_composite(self):
        ccc = self._collection(["gene_symbol", "spliceai_pred_dp_dl"])
        self.assertEqual(["gene_symbol", "spliceai"], self._collapse(ccc))

    def test_composite_already_present_loses_the_stray_members(self):
        ccc = self._collection(["spliceai", "gene_symbol", "spliceai_pred_ds_al"])
        self.assertEqual(["spliceai", "gene_symbol"], self._collapse(ccc))

    def test_collection_without_the_group_is_untouched(self):
        ccc = self._collection(["gene_symbol", "hgvs_c"])
        version = ccc.version_id
        self.assertEqual(["gene_symbol", "hgvs_c"], self._collapse(ccc))
        self.assertEqual(version, ccc.version_id, "Nothing changed, so nothing to invalidate")

    def test_is_idempotent(self):
        ccc = self._collection(["spliceai_max_ds", "gene_symbol"])
        self.assertEqual(["spliceai", "gene_symbol"], self._collapse(ccc))
        version = ccc.version_id
        self.assertEqual(["spliceai", "gene_symbol"], self._collapse(ccc))
        self.assertEqual(version, ccc.version_id)


class CustomColumnsArrangePageTest(TestCase):
    """ The arrange page offers composites; the columns they draw are there but out of the way """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_custom_columns_arrange")[0]
        cls.ccc = CustomColumnsCollection.objects.create(name="test arrange page", user=cls.user)

    def test_composites_and_members_are_marked_up(self):
        client = Client()
        client.force_login(self.user)
        url = reverse("view_custom_columns", kwargs={"custom_columns_collection_id": self.ccc.pk})
        response = client.get(url)
        self.assertEqual(200, response.status_code)
        content = response.content.decode()
        self.assertIn('composite-column cursor-move" column_id="spliceai"', content)
        self.assertIn('composite-member cursor-move" column_id="spliceai_max_ds"', content)
        self.assertIn("Show columns already inside a composite", content)
