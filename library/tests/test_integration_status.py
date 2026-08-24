from contextlib import contextmanager
from datetime import timedelta

from django.test import TestCase
from django.utils.timezone import now

from library.integration_status import (
    IntegrationDetail,
    IntegrationStatus,
    IntegrationTrigger,
    get_integration_statuses,
    integration_status_signal,
    run_integration_trigger,
)

TEST_NAMES = {"Alpha", "Beta", "Gamma"}


@contextmanager
def temporarily_connected(*receivers):
    for receiver in receivers:
        integration_status_signal.connect(receiver)
    try:
        yield
    finally:
        for receiver in receivers:
            integration_status_signal.disconnect(receiver)


def _test_status_names(statuses: list[IntegrationStatus]) -> list[str]:
    """ Deployment providers are registered too - only look at what this test contributed """
    return [status.name for status in statuses if status.name in TEST_NAMES]


class IntegrationStatusCollectionTest(TestCase):

    def test_flattens_lists_and_drops_none(self):
        def single(sender, **kwargs):
            return IntegrationStatus(name="Gamma")

        def several(sender, **kwargs):
            return [IntegrationStatus(name="Alpha"), IntegrationStatus(name="Beta")]

        def nothing_configured(sender, **kwargs):
            return None

        with temporarily_connected(single, several, nothing_configured):
            statuses, messages = get_integration_statuses()

        self.assertEqual(["Alpha", "Beta", "Gamma"], _test_status_names(statuses))
        self.assertEqual([], messages)

    def test_a_registered_status_supersedes_the_fallback_for_its_key(self):
        def shipped_fallback(sender, **kwargs):
            return IntegrationStatus(name="Alpha", key="alpha", fallback=True)

        def richer_provider(sender, **kwargs):
            return IntegrationStatus(name="Beta", key="alpha", record_count=7)

        with temporarily_connected(shipped_fallback, richer_provider):
            statuses, _ = get_integration_statuses()

        self.assertEqual(["Beta"], _test_status_names(statuses))

    def test_fallback_kept_when_nothing_claims_its_key(self):
        def shipped_fallback(sender, **kwargs):
            return IntegrationStatus(name="Alpha", key="alpha", fallback=True)

        with temporarily_connected(shipped_fallback):
            statuses, _ = get_integration_statuses()

        self.assertEqual(["Alpha"], _test_status_names(statuses))

    def test_raising_receiver_becomes_a_message(self):
        def broken(sender, **kwargs):
            raise ValueError("provider blew up")

        def working(sender, **kwargs):
            return IntegrationStatus(name="Alpha")

        with temporarily_connected(broken, working):
            statuses, messages = get_integration_statuses()

        self.assertEqual(["Alpha"], _test_status_names(statuses))
        self.assertEqual(1, len(messages))
        self.assertIn("provider blew up", messages[0])


class IntegrationStatusTriggerTest(TestCase):

    def setUp(self):
        self.runs = []
        self.status = IntegrationStatus(
            name="Alpha",
            trigger=IntegrationTrigger(action_id="alpha-sync", run=lambda: self.runs.append("ran"))
        )

    def _provider(self, sender, **kwargs):
        return self.status

    def test_runs_a_registered_action_id(self):
        with temporarily_connected(self._provider):
            triggered = run_integration_trigger("alpha-sync")
        self.assertEqual(self.status.name, triggered.name)
        self.assertEqual(["ran"], self.runs)

    def test_unregistered_action_id_runs_nothing(self):
        with temporarily_connected(self._provider):
            self.assertIsNone(run_integration_trigger("not-registered"))
        self.assertEqual([], self.runs)


class IntegrationStatusLastActivityTest(TestCase):

    def test_last_activity_is_the_newest_timestamp(self):
        right_now = now()
        status = IntegrationStatus(name="Alpha", details=[
            IntegrationDetail(label="Last Run", timestamp=right_now),
            IntegrationDetail(label="Last Success", timestamp=right_now - timedelta(days=3)),
            IntegrationDetail(label="Last Change", timestamp=None),
        ])
        self.assertEqual(right_now, status.last_activity)

    def test_no_timestamps_means_never_run(self):
        status = IntegrationStatus(name="Alpha", details=[IntegrationDetail(label="Last Run")])
        self.assertIsNone(status.last_activity)
        self.assertFalse(status.has_run)
