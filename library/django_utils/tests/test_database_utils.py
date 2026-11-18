from django.contrib.auth.models import Group, User
from django.core.exceptions import EmptyResultSet
from django.db import connection
from django.test import TestCase

from library.django_utils.database_utils import (
    get_pg_setting,
    get_queryset_row_estimate,
    get_table_row_estimates,
    iter_ordered,
    iter_ordered_batches,
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


class QuerysetToSQLTest(TestCase):
    def test_params_are_quoted_not_run(self):
        statement_timeout = get_pg_setting("statement_timeout")
        name = "x'\n; SET statement_timeout = 1234 --"
        sql = queryset_to_sql(GenomeBuild.objects.filter(name=name).values_list("name"))
        self.assertEqual(statement_timeout, get_pg_setting("statement_timeout"))
        with connection.cursor() as cursor:
            cursor.execute(sql)
            self.assertEqual([], cursor.fetchall())
        self.assertEqual(statement_timeout, get_pg_setting("statement_timeout"))


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


class IterOrderedTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        # last names tie, so the order across batches depends on the pk tiebreaker
        cls.users = [User.objects.create_user(username=f"iter_ordered_{i}", last_name=f"name_{i % 3}")
                     for i in range(7)]
        group_a, group_b = Group.objects.create(name="iter_ordered_a"), Group.objects.create(name="iter_ordered_b")
        for user in cls.users:
            user.groups.add(group_a, group_b)

    def _qs(self):
        return User.objects.filter(pk__in=[u.pk for u in self.users])

    def test_keeps_order_across_batches(self):
        qs = self._qs().order_by("-last_name")
        self.assertEqual(list(qs.order_by("-last_name", "pk")), list(iter_ordered(qs, batch_size=2)))
        self.assertEqual([2, 2, 2, 1], [len(batch_qs) for batch_qs in iter_ordered_batches(qs, batch_size=2)])

    def test_multi_valued_join_rows_as_the_queryset_gives_them(self):
        """ A repeated pk is fetched in one batch, so its rows aren't doubled across a batch boundary """
        qs = self._qs().filter(groups__name__startswith="iter_ordered").order_by("last_name")
        self.assertEqual(7 * 2, qs.count())
        self.assertEqual(list(qs.order_by("last_name", "pk")), list(iter_ordered(qs, batch_size=2)))
