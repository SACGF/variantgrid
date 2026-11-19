import math

from django.contrib.auth.models import Group, Permission, User
from django.core.exceptions import EmptyResultSet
from django.db import connection
from django.db.models import F, TextField
from django.db.models.functions import Cast
from django.test import TestCase
from django.test.utils import CaptureQueriesContext

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

    def test_objects_filters_run_once(self):
        """ The pk query applies the filters; batches fetch by pk, so a subquery in them isn't repeated """
        in_group_a = Group.objects.filter(name="iter_ordered_a").values("user")
        qs = User.objects.filter(pk__in=in_group_a).order_by("last_name")
        with CaptureQueriesContext(connection) as queries:
            users = list(iter_ordered(qs, batch_size=3))
        self.assertEqual(list(qs.order_by("last_name", "pk")), users)
        self.assertEqual(1 + math.ceil(len(users) / 3), len(queries))
        self.assertFalse(any("auth_group" in q["sql"] for q in queries[1:]))

    def test_objects_annotated_on_own_columns(self):
        """ As Guardian's get_objects_for_user: an annotation of the pk as text, filtered on """
        qs = self._qs().annotate(pk_text=Cast("pk", TextField())).filter(
            pk_text__in=[str(u.pk) for u in self.users[:5]]).order_by("last_name")
        with CaptureQueriesContext(connection) as queries:
            rows = [(u.pk, u.pk_text) for u in iter_ordered(qs, batch_size=2)]
        self.assertEqual([(u.pk, str(u.pk)) for u in qs.order_by("last_name", "pk")], rows)
        self.assertFalse(any(" IN ('" in q["sql"] for q in queries[1:]))

    def test_objects_annotated_through_a_join_keep_the_filters(self):
        qs = self._qs().filter(groups__name="iter_ordered_a").annotate(group_name=F("groups__name")).order_by("last_name")
        rows = [(u.pk, u.group_name) for u in iter_ordered(qs, batch_size=2)]
        self.assertEqual([(u.pk, "iter_ordered_a") for u in qs.order_by("last_name", "pk")], rows)

    def test_objects_keep_select_related_and_only(self):
        qs = Permission.objects.filter(content_type__app_label="auth").select_related("content_type").only(
            "codename", "content_type__app_label").order_by("codename")
        expected = [(p.codename, "auth") for p in qs.order_by("codename", "pk")]
        with self.assertNumQueries(1 + math.ceil(len(expected) / 10)):
            permissions = [(p.codename, p.content_type.app_label) for p in iter_ordered(qs, batch_size=10)]
        self.assertEqual(expected, permissions)

    def test_objects_through_multi_valued_filter_come_back_once(self):
        qs = self._qs().filter(groups__name__startswith="iter_ordered").order_by("last_name")
        self.assertEqual(7 * 2, qs.count())
        self.assertEqual(list(qs.distinct().order_by("last_name", "pk")), list(iter_ordered(qs, batch_size=2)))

    def test_values_keep_the_filters_on_a_multi_valued_join(self):
        """ The filter on groups decides which group each row has, so values() batches are fetched with it """
        qs = self._qs().filter(groups__name="iter_ordered_a").values("pk", "groups__name").order_by("last_name")
        rows = [row for batch_qs in iter_ordered_batches(qs, batch_size=2) for row in batch_qs]
        self.assertEqual(list(qs.order_by("last_name", "pk")), rows)
        self.assertEqual({"iter_ordered_a"}, {row["groups__name"] for row in rows})
