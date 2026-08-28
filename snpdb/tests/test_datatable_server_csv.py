from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db.models import QuerySet
from django.test import RequestFactory, TestCase

from snpdb.models import Cohort, GenomeBuild, ImportStatus
from snpdb.views.datatable_view import (
    DatatableConfig,
    RichColumn,
    SortOrder,
    datatable_csv_response,
    datatable_definition,
)


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


class DatatableServerCsvTests(TestCase):
    def setUp(self):
        self.user = User.objects.create(username='datatable_csv_user')
        genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        for name in ["zebra", "alpha"]:
            Cohort.objects.create(name=name, user=self.user, genome_build=genome_build,
                                  import_status=ImportStatus.SUCCESS, vcf=None)

    def _config(self, column_class=CohortCsvColumns, **params) -> DatatableConfig:
        request = RequestFactory().get("/fake/cohorts/", params)
        request.user = self.user
        return column_class(request)

    def _csv_lines(self, **params) -> list[str]:
        response = datatable_csv_response(self._config(**params))
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
        definition = datatable_definition(self._config())
        self.assertIn("dataTableCsv=1", definition["downloadUrl"])
        self.assertEqual(definition["csvName"], "cohorts")

        self.assertNotIn("downloadUrl", datatable_definition(self._config(column_class=NoCsvColumns)))

    def test_csv_refused_where_not_enabled(self):
        with self.assertRaises(PermissionDenied):
            datatable_csv_response(self._config(column_class=NoCsvColumns))
