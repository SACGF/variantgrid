from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db.models import QuerySet
from django.test import RequestFactory, TestCase

from snpdb.enums import ImportStatus
from snpdb.models import Cohort, GenomeBuild
from snpdb.views.datatable_view import CellData, DatabaseTableView, DatatableConfig, RichColumn, SortOrder


class CohortCsvColumns(DatatableConfig[Cohort]):
    server_csv_download = True
    csv_name = "cohorts"

    def __init__(self, request):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="id", label="ID", visible=False),
            RichColumn(key="name", label="Name", orderable=True, default_sort=SortOrder.ASC),
            RichColumn(key="import_status", label="Status", orderable=True),
            # An action column shares its key with the column it acts on, and has nothing to export
            RichColumn(key="id", name="delete", label="", renderer=self.render_delete),
        ]

    def get_initial_queryset(self) -> QuerySet[Cohort]:
        return Cohort.objects.all()

    def filter_queryset(self, qs: QuerySet[Cohort]) -> QuerySet[Cohort]:
        if name := self.get_query_param("name"):
            qs = qs.filter(name=name)
        return qs


class NoCsvColumns(CohortCsvColumns):
    server_csv_download = False


class CohortRenderedCsvColumns(DatatableConfig[Cohort]):
    """ What the user sees: a choice's label, and a server renderer that looks its value up per page """
    server_csv_download = True

    def __init__(self, request):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="name", label="Name", orderable=True, default_sort=SortOrder.ASC),
            RichColumn(key="import_status", label="Status",
                       **RichColumn.choices_client_renderer(ImportStatus.choices)),
            RichColumn(key="id", name="upper_name", label="Upper Name", renderer=self.render_upper_name,
                       csv_rendered=True),
        ]

    def pre_render(self, qs, rows):
        super().pre_render(qs, rows)
        self._page_upper_names = {row["id"]: row["name"].upper() for row in rows}

    def render_upper_name(self, cell: CellData):
        return self._page_upper_names[cell.value]

    def get_initial_queryset(self) -> QuerySet[Cohort]:
        return Cohort.objects.all()


class DatatableServerCsvTests(TestCase):
    def setUp(self):
        self.user = User.objects.create(username='datatable_csv_user')
        genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        for name in ["zebra", "alpha"]:
            Cohort.objects.create(name=name, user=self.user, genome_build=genome_build,
                                  import_status=ImportStatus.SUCCESS, vcf=None)

    def _view(self, column_class=CohortCsvColumns, **params) -> DatabaseTableView:
        request = RequestFactory().get("/fake/cohorts/", params)
        request.user = self.user
        view = DatabaseTableView(column_class=column_class)
        view.request = request
        view.config = column_class(request)
        return view

    def _csv_lines(self, **params) -> list[str]:
        response = self._view(**params).download_csv()
        return b"".join(response.streaming_content).decode().strip().splitlines()

    def test_csv_columns_deduplicated(self):
        # 'id' is both the hidden column and the delete column's key - one column, under its own label.
        # Excel reads a file starting with the letters ID as SYLK, so the header goes out quoted
        self.assertEqual(self._csv_lines()[0], '"ID","Name","Status"')

    def test_csv_applies_default_sort(self):
        self.assertEqual([line.split(",")[1] for line in self._csv_lines()[1:]], ["alpha", "zebra"])

    def test_csv_applies_filters(self):
        lines = self._csv_lines(name="zebra")
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[1].split(",")[1], "zebra")

    def test_definition_offers_download_url(self):
        definition = self._view().json_definition()
        self.assertIn("dataTableCsv=1", definition["downloadUrl"])

        self.assertNotIn("downloadUrl", self._view(column_class=NoCsvColumns).json_definition())

    def test_csv_writes_what_the_grid_shows(self):
        lines = self._csv_lines(column_class=CohortRenderedCsvColumns)
        self.assertEqual(lines, ["Name,Status,Upper Name", "alpha,success,ALPHA", "zebra,success,ZEBRA"])

    def test_csv_refused_where_not_enabled(self):
        with self.assertRaises(PermissionDenied):
            self._view(column_class=NoCsvColumns).download_csv()
