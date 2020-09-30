from library.jqgrid_user_row_config import JqGridUserRowConfig
from pathtests.models import PathologyTestOrder, Case, PathologyTest


class PathologyTestOrdersGrid(JqGridUserRowConfig):
    model = PathologyTestOrder
    caption = "Pathology Test Orders"
    fields = ["id", "external_pk__code", "case__external_pk__code", "pathology_test_version__pathology_test__name",
              "user", "created", "modified", "started_library", "finished_library", "started_sequencing",
              "finished_sequencing", "order_completed", "experiment", "sequencing_run"]
    colmodel_overrides = {'id': {"width": 30, 'formatter': 'viewPathologyOrderLink'}}

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())


class CasesGrid(JqGridUserRowConfig):
    model = Case
    caption = "Cases"
    fields = ["id", "external_pk__code", "name", "lead_scientist",
              "created", "result_required_date", "modified", "patient__external_pk__code",
              "report_date", "details", "status", "workflow_status", "investigation_type"]
    colmodel_overrides = {'id': {"width": 30, 'formatter': 'viewCaseLink'}}

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        self.queryset = queryset.values(*self.get_field_names())


class PathologyTestsGrid(JqGridUserRowConfig):
    model = PathologyTest
    caption = "Pathology Tests"
    fields = ["name", "curator__username",
              'activepathologytestversion__pathology_test_version',
              'activepathologytestversion__pathology_test_version__version',
              'modified']
    colmodel_overrides = {'name': {"width": 150, 'formatter': 'viewPathologyTestLink'},
                          'curator__username': {"label": 'Curator'},
                          'activepathologytestversion__pathology_test_version': {"hidden": True},
                          'activepathologytestversion__pathology_test_version__version': {"label": "Active Version",
                                                                                          "formatter": "viewPathologyTestVersionLink"}}

    def __init__(self, user, **kwargs):
        super().__init__(user)
        queryset = self.model.objects.all()
        queryset = queryset.filter(deleted=False)
        self.queryset = queryset.values(*self.get_field_names())

        self.extra_config.update({'sortname': 'modified',
                                  'sortorder': 'desc'})
