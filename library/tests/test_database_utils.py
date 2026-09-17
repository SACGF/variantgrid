from django.core.exceptions import EmptyResultSet
from django.test import TestCase

from library.utils.database_utils import queryset_to_sql, render_empty_result_set_sql
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
