from typing import TypeVar
from unittest import mock

from django.contrib.auth.models import User
from django.db.models import Model, QuerySet
from django.test import RequestFactory, SimpleTestCase

from snpdb.models import Cohort, Trio
from snpdb.views.datatable_view import DatatableConfig, RichColumn, datatable_definition

M = TypeVar("M", bound=Model)


class DirectColumns(DatatableConfig[Cohort]):
    rich_columns = [RichColumn(key="id", label="ID")]


class IntermediateColumns(DatatableConfig[M]):
    rich_columns = [RichColumn(key="id", label="ID")]


class ThroughIntermediateColumns(IntermediateColumns[Trio]):
    pass


class ThroughConcreteBaseColumns(DirectColumns):
    pass


class BareColumns(DatatableConfig):
    rich_columns = [RichColumn(key="id", label="ID")]

    def get_initial_queryset(self) -> QuerySet[Cohort]:
        return Cohort.objects.all()


class DatatableModelTests(SimpleTestCase):
    """ The model behind a grid comes from the DatatableConfig[Model] declaration, so naming the CSV in
        the definition doesn't build the initial queryset (#1913) """

    def _config(self, column_class):
        request = RequestFactory().get("/fake/")
        request.user = User(username="datatable_model_user")
        return column_class(request)

    def test_model_from_declaration(self):
        self.assertIs(DirectColumns._declared_model(), Cohort)
        self.assertIs(ThroughIntermediateColumns._declared_model(), Trio)
        self.assertIs(ThroughConcreteBaseColumns._declared_model(), Cohort)
        self.assertIsNone(IntermediateColumns._declared_model())
        self.assertIsNone(BareColumns._declared_model())

    def test_definition_does_not_build_queryset(self):
        config = self._config(DirectColumns)
        with mock.patch.object(DirectColumns, "get_initial_queryset", side_effect=AssertionError("built")):
            self.assertEqual(datatable_definition(config)["csvName"], "Cohort")

    def test_bare_config_falls_back_with_warning(self):
        config = self._config(BareColumns)
        with self.assertLogs("snpdb.views.datatable_view", level="WARNING") as logs:
            self.assertIs(config._model, Cohort)
        self.assertIn("BareColumns", logs.output[0])
