from django.core.exceptions import EmptyResultSet
from django.test import TestCase

from library.django_utils.database_utils import (
    get_pg_setting,
    get_queryset_row_estimate,
    get_table_row_estimates,
    pg_settings,
    queryset_to_sql,
    render_empty_result_set_sql,
)
from snpdb.models import GenomeBuild


class RenderEmptyResultSetSQLTest(TestCase):
    def test_empty_in_renders_sql(self):
        qs = GenomeBuild.objects.filter(pk__in=[])
        with render_empty_result_set_sql():
            sql = queryset_to_sql(qs)
        self.assertIn("IN (NULL)", sql)

    def test_none_renders_sql(self):
        with render_empty_result_set_sql():
            sql = queryset_to_sql(GenomeBuild.objects.none())
        self.assertIn("0 = 1", sql)

    def test_patches_restored(self):
        with render_empty_result_set_sql():
            pass
        with self.assertRaises(EmptyResultSet):
            queryset_to_sql(GenomeBuild.objects.filter(pk__in=[]))


class PgSettingsTest(TestCase):
    def test_nested_restores_the_outer_value(self):
        before = get_pg_setting("statement_timeout")
        with pg_settings(statement_timeout="7s"):
            with pg_settings(statement_timeout=5000):
                self.assertEqual("5s", get_pg_setting("statement_timeout"))
            self.assertEqual("7s", get_pg_setting("statement_timeout"))
        self.assertEqual(before, get_pg_setting("statement_timeout"))


class RowEstimatesTest(TestCase):
    def test_table_row_estimates_only_named_tables(self):
        estimates = get_table_row_estimates(["snpdb_genomebuild", "not_a_table"])
        self.assertEqual({"snpdb_genomebuild"}, set(estimates))

    def test_queryset_row_estimate(self):
        self.assertIsInstance(get_queryset_row_estimate(GenomeBuild.objects.all()), int)
