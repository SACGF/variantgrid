import pandas as pd
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.test import RequestFactory, TestCase
from django.utils.html import format_html
from django.utils.safestring import mark_safe

from snpdb.models import Cohort, GenomeBuild, ImportStatus
from snpdb.views.datatable_dataframe import (
    INDEX_COLUMN_KEY,
    DataFrameDatatableConfig,
    DataFrameTableView,
)
from snpdb.views.datatable_view import DatabaseTableView, DatatableConfig, RichColumn

PAYLOAD = "<img src=x onerror=alert(1)>"
ESCAPED_PAYLOAD = "&lt;img src=x onerror=alert(1)&gt;"


class EscapingColumns(DatatableConfig[Cohort]):
    server_csv_download = True

    def __init__(self, request):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="name", name="plain", label=PAYLOAD),
            RichColumn(key="name", name="client_rendered", label=mark_safe("Client<br/>rendered"),
                       client_renderer="TableFormat.text"),
            RichColumn(key="name", name="safe_html", renderer=self.render_safe_html, csv_rendered=True),
        ]

    @staticmethod
    def render_safe_html(cell):
        return format_html("<b>{}</b>", cell.value)

    def get_initial_queryset(self) -> QuerySet[Cohort]:
        return Cohort.objects.all()


class EscapingDataFrameConfig(DataFrameDatatableConfig):
    def get_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame({PAYLOAD: [PAYLOAD]}, index=[PAYLOAD])


class DatatableEscapingTests(TestCase):
    """ A grid cell or label goes into the page as HTML - escaped unless a client renderer takes the raw
        value or the server renderer built HTML on purpose. The CSV keeps the raw value """

    def setUp(self):
        self.user = User.objects.create(username='datatable_escaping_user')
        Cohort.objects.create(name=PAYLOAD, user=self.user, genome_build=GenomeBuild.get_name_or_alias("GRCh38"),
                              import_status=ImportStatus.SUCCESS, vcf=None)

    def _request(self):
        request = RequestFactory().get("/fake/")
        request.user = self.user
        return request

    def _view(self) -> DatabaseTableView:
        view = DatabaseTableView(column_class=EscapingColumns)
        view.request = self._request()
        view.config = EscapingColumns(view.request)
        return view

    def test_cells(self):
        [row] = self._view().get_context_data()["data"]
        self.assertEqual(row["plain"], ESCAPED_PAYLOAD)
        self.assertEqual(row["client_rendered"], PAYLOAD)
        self.assertEqual(row["safe_html"], f"<b>{ESCAPED_PAYLOAD}</b>")

    def test_labels(self):
        labels = [c["label"] for c in self._view().json_definition()["columns"]]
        self.assertEqual(labels, [ESCAPED_PAYLOAD, "Client<br/>rendered", "Safe Html"])

    def test_csv_is_raw(self):
        response = self._view().download_csv()
        lines = b"".join(response.streaming_content).decode().splitlines()
        self.assertEqual(lines[1].split(",")[0], PAYLOAD)

    def test_dataframe(self):
        view = DataFrameTableView(column_class=EscapingDataFrameConfig)
        view.request = self._request()
        view.config = EscapingDataFrameConfig(view.request)
        [row] = view.get_context_data()["data"]
        self.assertEqual(row[INDEX_COLUMN_KEY], ESCAPED_PAYLOAD)
        self.assertEqual(row["c1"], ESCAPED_PAYLOAD)
        self.assertEqual(view.json_definition()["columns"][1]["label"], ESCAPED_PAYLOAD)
