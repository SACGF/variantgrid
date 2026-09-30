import json
from collections import defaultdict

from django.contrib.auth.models import User
from django.urls import reverse

from library.django_utils.unittest_utils import URLTestCase


class AnnotationDescriptionsTest(URLTestCase):
    """ The composite cells card - @see annotation.views_descriptions.view_annotation_descriptions """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username=f"test_user_{__file__}")[0]

    def _get_context(self):
        self.client.force_login(self.user)
        response = self.client.get(reverse("view_annotation_descriptions"))
        self.assertEqual(response.status_code, 200)
        return response.context

    def _get_rows(self) -> dict[str, list[dict]]:
        """ column -> its rows, one per columns version range, from the composite sections and level cards """
        context = self._get_context()
        rows = defaultdict(list)
        for section in context["composite_sections"]:
            for row in section["rows"]:
                rows[row["column"].pk].append(row)
        for level_rows in context["columns_by_annotation_level"].values():
            for row in level_rows:
                rows[row["column"].pk].append(row)
        return rows

    def test_a_column_lists_a_row_per_columns_version(self):
        """ CADD phred moved from dbNSFP v1 to v4 - the reader picking a version should find the source
            that version reads, not whichever definition was declared first """
        versions = {(row["min_columns_version"], row["max_columns_version"]) for row in self._get_rows()["cadd_phred"]}
        self.assertEqual(versions, {(None, 1), (4, None)})

    def test_sources_for_columns_no_tool_writes(self):
        rows = self._get_rows()
        self.assertEqual(rows["chrom"][0]["source"], "Variant")
        self.assertEqual(rows["hgvs_g"][0]["source"], "Variant")
        self.assertEqual(rows["clinvar_review_status"][0]["source"], "ClinVar")
        self.assertEqual(rows["max_internal_classification"][0]["source"], "VariantGrid")
        self.assertEqual(rows["max_internal_classification"][0]["category"].name, "CLASSIFICATIONS")
        self.assertEqual(rows["clinvar_review_status"][0]["category"].name, "CLASSIFICATIONS")

    def test_example_rows_carry_every_member_the_cell_draws(self):
        """ The column definition and the example row are built separately, and a member this build
            doesn't annotate has to drop out of both - otherwise the cell draws detail the grid can't """
        for section in self._get_context()["composite_sections"]:
            column, row = json.loads(section["column_json"]), json.loads(section["row_json"])
            with self.subTest(composite=section["column"].pk):
                for member in column.get("renderKwargs", {}).get("members", []):
                    self.assertIn(member["path"], row)
                for entry in column.get("sortMenu") or []:
                    self.assertIn(entry["column"], row)
