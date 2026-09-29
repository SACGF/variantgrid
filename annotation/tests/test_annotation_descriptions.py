import html
import json
import re
from collections import defaultdict

from django.contrib.auth.models import User
from django.urls import reverse

from library.django_utils.unittest_utils import URLTestCase


class AnnotationDescriptionsTest(URLTestCase):
    """ The composite cells card - @see annotation.views_descriptions.view_annotation_descriptions """

    SECTION_RE = re.compile(
        r'id="composite-(?P<pk>\w+)".*?'
        r'data-column-json="(?P<column>[^"]*)"\s*'
        r'data-row-json="(?P<row>[^"]*)"', re.DOTALL)
    ROW_RE = re.compile(
        r"""<tr class='annotation-row' data-min-cv="(?P<min>\d*)" data-max-cv="(?P<max>\d*)">"""
        r'.*?annotation-grid-column-name">(?P<pk>\w+)<'
        r'(?P<body>.*?)</tr>', re.DOTALL)

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username=f"test_user_{__file__}")[0]

    def _get_sections(self) -> dict[str, tuple[dict, dict]]:
        """ composite -> (the column definition, the example row) the page hands the client renderer """
        self.client.force_login(self.user)
        response = self.client.get(reverse("view_annotation_descriptions"))
        self.assertEqual(response.status_code, 200)
        sections = {m["pk"]: (json.loads(html.unescape(m["column"])), json.loads(html.unescape(m["row"])))
                    for m in self.SECTION_RE.finditer(response.content.decode())}
        self.assertTrue(sections, "No composite cell sections on the page")
        return sections

    def _get_rows(self) -> dict[str, list[dict]]:
        """ column -> its rows, one per columns version range, as the page's version toggle reads them """
        self.client.force_login(self.user)
        response = self.client.get(reverse("view_annotation_descriptions"))
        self.assertEqual(response.status_code, 200)
        rows = defaultdict(list)
        for m in self.ROW_RE.finditer(response.content.decode()):
            rows[m["pk"]].append({"min": m["min"], "max": m["max"], "body": m["body"]})
        self.assertTrue(rows, "No annotation rows on the page")
        return rows

    def test_a_column_lists_a_row_per_columns_version(self):
        """ CADD phred moved from dbNSFP v1 to v4 - the reader picking a version should find the source
            that version reads, not whichever definition was declared first """
        versions = {(row["min"], row["max"]) for row in self._get_rows()["cadd_phred"]}
        self.assertEqual(versions, {("", "1"), ("4", "")})

    def test_sources_for_columns_no_tool_writes(self):
        rows = self._get_rows()
        self.assertIn("annotation-source-variant", rows["chrom"][0]["body"])
        self.assertIn("annotation-source-variant", rows["hgvs_g"][0]["body"])
        self.assertIn("annotation-source-clinvar", rows["clinvar_review_status"][0]["body"])
        self.assertIn("annotation-source-variantgrid", rows["max_internal_classification"][0]["body"])
        self.assertIn("annotation-category-classifications", rows["max_internal_classification"][0]["body"])
        self.assertIn("annotation-category-classifications", rows["clinvar_review_status"][0]["body"])

    def test_example_rows_carry_every_member_the_cell_draws(self):
        """ The column definition and the example row are built separately, and a member this build
            doesn't annotate has to drop out of both - otherwise the cell draws detail the grid can't """
        for composite_id, (column, row) in self._get_sections().items():
            with self.subTest(composite=composite_id):
                for member in column.get("renderKwargs", {}).get("members", []):
                    self.assertIn(member["path"], row)
                for entry in column.get("sortMenu") or []:
                    self.assertIn(entry["column"], row)

    def test_gnomad_example_cell(self):
        column, row = self._get_sections()["gnomad"]
        self.assertEqual(column["render"], "VariantGridFormat.gnomad")
        self.assertIn("variantannotation__gnomad_popmax_af",
                      {entry["column"] for entry in column["sortMenu"]})
        self.assertTrue(row["variantannotation__gnomad_af"])
