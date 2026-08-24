from unittest.mock import patch

from django.test import TestCase, override_settings

from snpdb.models import JobsControl
from snpdb.signals.jobs_autopause import autopause_jobs_after_reboot

BOOT = 1_700_000_000


@override_settings(JOBS_AUTOPAUSE_ON_REBOOT=True, JOBS_AUTOPAUSE_ON_REBOOT_UPTIME_SECS=600,
                   JOBS_AUTOPAUSE_ON_REBOOT_WINDOW_SECS=3600)
class JobsAutopauseTest(TestCase):
    @staticmethod
    def _worker_ready(boot_time: int, uptime: float = 30):
        with patch("snpdb.signals.jobs_autopause._read_proc_boot_time", return_value=boot_time), \
                patch("snpdb.signals.jobs_autopause._read_proc_uptime_seconds", return_value=uptime):
            autopause_jobs_after_reboot()

    def test_single_reboot_keeps_running(self):
        self._worker_ready(BOOT)
        self.assertFalse(JobsControl.is_paused())

    def test_second_reboot_within_window_pauses(self):
        self._worker_ready(BOOT)
        self._worker_ready(BOOT + 900)
        self.assertTrue(JobsControl.is_paused())

    def test_second_reboot_after_window_keeps_running(self):
        self._worker_ready(BOOT)
        self._worker_ready(BOOT + 86400)
        self.assertFalse(JobsControl.is_paused())

    def test_worker_restart_within_boot_does_not_re_pause(self):
        self._worker_ready(BOOT)
        self._worker_ready(BOOT + 900)
        JobsControl.resume(by="test")
        self._worker_ready(BOOT + 900)
        self.assertFalse(JobsControl.is_paused())

    def test_long_uptime_records_boot_without_pausing(self):
        self._worker_ready(BOOT)
        self._worker_ready(BOOT + 900, uptime=5000)
        self.assertFalse(JobsControl.is_paused())
        self._worker_ready(BOOT + 1800)
        self.assertTrue(JobsControl.is_paused())
