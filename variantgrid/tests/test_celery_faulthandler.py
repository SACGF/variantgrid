import os
import subprocess
import tempfile
from unittest.mock import patch

from django.test import SimpleTestCase

from variantgrid.celery import report_dead_worker_fault_dumps


class ReportDeadWorkerFaultDumpsTest(SimpleTestCase):
    def test_reports_dead_process_dumps_and_leaves_live_ones(self):
        dead_pid = subprocess.Popen(["true"]).pid
        os.waitpid(dead_pid, 0)
        empty_dead_pid = subprocess.Popen(["true"]).pid
        os.waitpid(empty_dead_pid, 0)

        with tempfile.TemporaryDirectory() as fault_dir:
            def dump_file(pid, text):
                with open(os.path.join(fault_dir, f"faulthandler_{pid}.log"), "w") as f:
                    f.write(text)

            dump_file(dead_pid, "Fatal Python error: Aborted\n\nCurrent thread 0x1 (most recent call first):\n")
            dump_file(empty_dead_pid, "")
            dump_file(os.getpid(), "")

            with patch("variantgrid.celery.rollbar.report_message") as report_message:
                report_dead_worker_fault_dumps(fault_dir)

            report_message.assert_called_once()
            self.assertEqual("Celery pool process died: Fatal Python error: Aborted", report_message.call_args.args[0])
            self.assertEqual(dead_pid, report_message.call_args.kwargs["extra_data"]["pid"])
            self.assertEqual([f"faulthandler_{os.getpid()}.log"], os.listdir(fault_dir))
