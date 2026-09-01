from typing import Optional

from django.db.models import QuerySet
from django.http import HttpRequest
from django.urls import reverse

from library.utils import JsonDataType
from pathtests.models import Case, PathologyTest, PathologyTestOrder
from pathtests.models_enums import CaseState, CaseWorkflowStatus, InvestigationType
from snpdb.views.datatable_view import CellData, DatatableConfig, RichColumn, SortOrder


class PathologyTestOrdersColumns(DatatableConfig[PathologyTestOrder]):
    grid_name = "Pathology Test Orders"
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC,
                       renderer=self.view_primary_key, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="external_pk__code", label="External Code", orderable=True),
            RichColumn(key="case__external_pk__code", label="Case", orderable=True),
            RichColumn(key="pathology_test_version__pathology_test__name", label="Pathology Test", orderable=True),
            RichColumn(key="user__username", label="User", orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="modified", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="started_library", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="finished_library", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="started_sequencing", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="finished_sequencing", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="order_completed", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="experiment__name", label="Experiment", orderable=True),
            RichColumn(key="sequencing_run__name", label="Sequencing Run", orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[PathologyTestOrder]:
        return PathologyTestOrder.objects.all()


class CasesColumns(DatatableConfig[Case]):
    grid_name = "Cases"
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC,
                       renderer=self.view_primary_key, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="external_pk__code", label="External Code", orderable=True),
            RichColumn(key="name", label="Name", orderable=True),
            RichColumn(key="lead_scientist__username", label="Lead Scientist", orderable=True,
                       extra_columns=["lead_scientist__id"], renderer=self.render_user),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="result_required_date", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="modified", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="patient__external_pk__code", label="Patient", orderable=True),
            RichColumn(key="report_date", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="details", label="Details", orderable=True),
            RichColumn(key="status", label="Status", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(CaseState.choices)),
            RichColumn(key="workflow_status", label="Workflow Status", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(CaseWorkflowStatus.choices)),
            RichColumn(key="investigation_type", label="Investigation Type", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(InvestigationType.choices)),
        ]

    def get_initial_queryset(self) -> QuerySet[Case]:
        return Case.objects.all()


class PathologyTestsColumns(DatatableConfig[PathologyTest]):
    grid_name = "Pathology Tests"
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="name", label="Name", orderable=True,
                       renderer=self.view_primary_key, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="curator__username", label="Curator", orderable=True,
                       extra_columns=["curator__id"], renderer=self.render_user),
            RichColumn(key="activepathologytestversion__pathology_test_version__version", label="Active Version",
                       orderable=True, renderer=self._render_active_version,
                       extra_columns=["activepathologytestversion__pathology_test_version__id"],
                       client_renderer='renderOptionalLink'),
            RichColumn(key="modified", client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
        ]

    @staticmethod
    def _render_active_version(cell: CellData) -> Optional[JsonDataType]:
        if version := cell.value:
            pk = cell["activepathologytestversion__pathology_test_version__id"]
            return {
                "text": f"Version {version}",
                "url": reverse("view_pathology_test_version", kwargs={"pk": pk}),
            }
        return None

    def get_initial_queryset(self) -> QuerySet[PathologyTest]:
        return PathologyTest.objects.filter(deleted=False)
