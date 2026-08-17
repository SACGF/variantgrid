from django.contrib.auth.models import User
from django.test import TestCase

from library.jqgrid.jqgrid import KnownCountPaginator


class KnownCountPaginatorTest(TestCase):
    """ Grids that already know their row count hand it to the paginator instead of having it
        re-derived with a COUNT(*) """

    def test_supplied_count_is_used_without_querying(self):
        paginator = KnownCountPaginator(User.objects.order_by("pk"), 10, count=1234)
        with self.assertNumQueries(0):
            self.assertEqual(1234, paginator.count)
            self.assertEqual(124, paginator.num_pages)
