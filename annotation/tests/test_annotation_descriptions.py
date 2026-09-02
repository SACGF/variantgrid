import html
import json
import re

from django.contrib.auth.models import User
from django.urls import reverse

from library.django_utils.unittest_utils import URLTestCase


class AnnotationDescriptionsTest(URLTestCase):
    """ The composite cells card - @see annotation.views_descriptions.view_annotation_descriptions """

    # The response goes through htmlmin, which picks the attribute quoting - so match either
    SECTION_RE = re.compile(
        r"""id="composite-(?P<pk>\w+)".*?"""
        r"""data-column-json=(?P<cq>["'])(?P<column>.*?)(?P=cq)\s*"""
        r"""data-row-json=(?P<rq>["'])(?P<row>.*?)(?P=rq)""", re.DOTALL)

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
