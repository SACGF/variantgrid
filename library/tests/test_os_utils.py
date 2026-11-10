import sys
import time

from django.test import SimpleTestCase

from library.utils.os_utils import execute_cmd


def _process_is_running(pid: int) -> bool:
    try:
        with open(f"/proc/{pid}/stat", encoding="utf-8") as f:
            return f.read().rsplit(")", 1)[1].split()[0] != "Z"  # a zombie has exited
    except FileNotFoundError:
        return False


class ExecuteCmdTest(SimpleTestCase):
    def test_max_rss_is_the_largest_descendant(self):
        """ #2749: the memory that matters is a tool's child (AnnotSV's bedtools), not the tool itself """
        megabytes = 150
        allocate = f"data = b'x' * ({megabytes} * 2 ** 20)"
        cmd_output = execute_cmd(["sh", "-c", f"{sys.executable} -c \"{allocate}\"; true"])
        self.assertEqual(cmd_output.return_code, 0)
        self.assertGreaterEqual(cmd_output.max_rss_mb, megabytes)

    def test_timeout_kills_the_processes_the_command_started(self):
        return_code, std_out, _std_err = execute_cmd(["sh", "-c", "sleep 60 & echo $!; wait"], timeout=0.5)
        self.assertNotEqual(return_code, 0)
        child_pid = int(std_out)
        deadline = time.monotonic() + 5
        while _process_is_running(child_pid) and time.monotonic() < deadline:
            time.sleep(0.05)
        self.assertFalse(_process_is_running(child_pid), "the command's child outlived the timeout")
