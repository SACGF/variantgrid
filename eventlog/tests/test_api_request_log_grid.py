from datetime import date, datetime

from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.utils import timezone

from eventlog.grids import ApiRequestLogColumns
from eventlog.models import ApiRequestLog


class ApiRequestLogGridDateFilterTest(TestCase):
    """ date_from / date_to on the API Log grid are inclusive local (TIME_ZONE) days """

    @classmethod
    def setUpTestData(cls):
        cls.admin = User.objects.create_superuser("api_log_grid_admin")
        tz = timezone.get_current_timezone()
        cls.rows_by_day = {}
        for day, local_time in [(date(2026, 3, 1), datetime(2026, 3, 1, 23, 30)),
                                (date(2026, 3, 2), datetime(2026, 3, 2, 0, 30)),
                                (date(2026, 3, 3), datetime(2026, 3, 3, 12, 0))]:
            row = ApiRequestLog.objects.create(method="GET", path="/x", view_name="v", status_code=200, duration_ms=1)
            ApiRequestLog.objects.filter(pk=row.pk).update(created=timezone.make_aware(local_time, tz))
            cls.rows_by_day[day] = row

    def _filtered_days(self, **params) -> set[date]:
        request = RequestFactory().get("/fake/api_request_log/datatable", params)
        request.user = self.admin
        config = ApiRequestLogColumns(request)
        pks = set(config.filter_queryset(config.get_initial_queryset()).values_list("pk", flat=True))
        return {day for day, row in self.rows_by_day.items() if row.pk in pks}

    def test_one_day(self):
        self.assertEqual({date(2026, 3, 2)}, self._filtered_days(date_from="2026-03-02", date_to="2026-03-02"))

    def test_open_ended(self):
        self.assertEqual({date(2026, 3, 2), date(2026, 3, 3)}, self._filtered_days(date_from="2026-03-02"))
        self.assertEqual({date(2026, 3, 1), date(2026, 3, 2)}, self._filtered_days(date_to="2026-03-02"))

    def test_bad_or_empty_date_is_no_filter(self):
        all_days = set(self.rows_by_day)
        self.assertEqual(all_days, self._filtered_days(date_from="", date_to=""))
        self.assertEqual(all_days, self._filtered_days(date_from="yesterday"))
