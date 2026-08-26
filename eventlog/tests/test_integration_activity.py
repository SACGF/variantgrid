from datetime import timedelta

from django.test import TestCase
from django.utils.timezone import now

from eventlog.models import IntegrationActivity, IntegrationActivityStatus
from eventlog.signals.integration_activity_status import integration_activity_status


class IntegrationActivityTrackTest(TestCase):

    def test_clean_run_stamps_success_without_change(self):
        with IntegrationActivity.track("test-integration", name="Test Integration"):
            pass

        activity = IntegrationActivity.objects.get(key="test-integration")
        self.assertIsNotNone(activity.last_attempt)
        self.assertIsNotNone(activity.last_success)
        self.assertIsNone(activity.last_change)
        self.assertIsNone(activity.last_error)

    def test_changed_body_stamps_last_change(self):
        with IntegrationActivity.track("test-integration", name="Test Integration") as activity:
            activity.changed = True

        activity = IntegrationActivity.objects.get(key="test-integration")
        self.assertIsNotNone(activity.last_change)

    def test_exception_stamps_error_and_re_raises(self):
        with self.assertRaises(ValueError):
            with IntegrationActivity.track("test-integration", name="Test Integration"):
                raise ValueError("upstream is down")

        activity = IntegrationActivity.objects.get(key="test-integration")
        self.assertIsNotNone(activity.last_attempt)
        self.assertIsNone(activity.last_success)
        self.assertIsNotNone(activity.last_error)
        self.assertEqual("upstream is down", activity.last_error_message)

    def test_counters_accumulate_over_repeated_calls(self):
        for _ in range(2):
            with IntegrationActivity.track("test-integration", name="Test Integration"):
                pass
        with self.assertRaises(ValueError):
            with IntegrationActivity.track("test-integration", name="Test Integration"):
                raise ValueError("upstream is down")

        activity = IntegrationActivity.objects.get(key="test-integration")
        self.assertEqual(3, activity.attempt_count)
        self.assertEqual(2, activity.success_count)
        self.assertEqual(1, activity.error_count)


class IntegrationActivityStatusReceiverTest(TestCase):

    def _last_failure_status(self, last_success, last_error) -> str:
        IntegrationActivity.objects.create(
            key="test-integration",
            name="Test Integration",
            last_attempt=last_error,
            last_success=last_success,
            last_error=last_error,
            last_error_message="upstream is down"
        )
        statuses = integration_activity_status(sender=None)
        details = {detail.label: detail for detail in statuses[0].details}
        return details["Last Failure"].status

    def test_failure_newer_than_success_is_danger(self):
        right_now = now()
        self.assertEqual("danger", self._last_failure_status(
            last_success=right_now - timedelta(hours=2), last_error=right_now))

    def test_failure_older_than_success_is_warning(self):
        right_now = now()
        self.assertEqual("warning", self._last_failure_status(
            last_success=right_now, last_error=right_now - timedelta(hours=2)))

    def test_recorded_error_renders_a_failing_row(self):
        IntegrationActivity.record("test-integration", name="Test Integration",
                                   status=IntegrationActivityStatus.ERROR, message="HTTP 500")
        status = integration_activity_status(sender=None)[0]
        self.assertEqual("Test Integration", status.name)
        self.assertFalse(status.ok)
        self.assertIsNotNone(status.dismiss)

    def test_dismissed_error_goes_quiet(self):
        IntegrationActivity.record("test-integration", name="Test Integration",
                                   status=IntegrationActivityStatus.ERROR, message="HTTP 500")
        IntegrationActivity.acknowledge_error("test-integration")

        status = integration_activity_status(sender=None)[0]
        details = {detail.label: detail for detail in status.details}
        failure = details["Last Failure"]
        self.assertEqual("secondary", failure.status)
        self.assertIsNone(failure.extra)
        self.assertTrue(status.ok)
        self.assertIsNone(status.dismiss)

    def test_new_error_after_dismissal_shows_again(self):
        IntegrationActivity.record("test-integration", name="Test Integration",
                                   status=IntegrationActivityStatus.ERROR, message="HTTP 500")
        IntegrationActivity.acknowledge_error("test-integration")
        IntegrationActivity.record("test-integration", name="Test Integration",
                                   status=IntegrationActivityStatus.ERROR, message="HTTP 400")

        status = integration_activity_status(sender=None)[0]
        self.assertFalse(status.ok)
        self.assertIsNotNone(status.dismiss)
