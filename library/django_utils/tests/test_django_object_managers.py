"""
The caching querysets directly: under UNIT_TEST the managers hand out plain querysets, so model tests never
reach this code.
"""
from django.contrib.auth.models import Group
from django.test import TestCase

from library.django_utils.django_object_managers import _CACHED_QUERIES, QuerySetImmutable


class QuerySetImmutableTests(TestCase):
    def setUp(self):
        _CACHED_QUERIES.clear()
        self.group_a = Group.objects.create(name="caching_test_a")
        self.group_b = Group.objects.create(name="caching_test_b")

    def tearDown(self):
        _CACHED_QUERIES.clear()

    def test_bare_get_is_cached(self):
        self.assertEqual(self.group_a, QuerySetImmutable(model=Group).get(name="caching_test_a"))
        with self.assertNumQueries(0):
            self.assertEqual(self.group_a, QuerySetImmutable(model=Group).get(name="caching_test_a"))

    def test_filters_are_not_dropped(self):
        QuerySetImmutable(model=Group).get(name="caching_test_a")
        with self.assertRaises(Group.DoesNotExist):
            QuerySetImmutable(model=Group).filter(name="caching_test_b").get(name="caching_test_a")

    def test_filtered_no_arg_gets_do_not_share_a_key(self):
        qs = QuerySetImmutable(model=Group)
        self.assertEqual(self.group_a, qs.filter(pk=self.group_a.pk).get())
        self.assertEqual(self.group_b, qs.filter(pk=self.group_b.pk).get())

    def test_values_is_not_served_a_cached_instance(self):
        QuerySetImmutable(model=Group).get(name="caching_test_a")
        values = QuerySetImmutable(model=Group).values("name").get(name="caching_test_a")
        self.assertEqual({"name": "caching_test_a"}, values)
