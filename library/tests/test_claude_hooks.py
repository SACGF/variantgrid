"""
The Bash PreToolUse hook (`.claude/hooks/pre_bash.py`): the regexes that decide which commands a shared box
prompts for, which are refused everywhere, and the per-host exceptions. Loaded by path because `.claude` is not a
package.
"""
import importlib.util
from pathlib import Path
from unittest import TestCase

HOOK = Path(__file__).resolve().parents[2] / ".claude" / "hooks" / "pre_bash.py"
_spec = importlib.util.spec_from_file_location("pre_bash", HOOK)
pre_bash = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pre_bash)


class PreBashHookTest(TestCase):
    def _decide(self, command, host="vg-test2", unguarded=False):
        return pre_bash.decision(command, host, unguarded)

    def test_plain_commands_pass(self):
        for command in ["git status", "python3 manage.py vg status", "python3 manage.py test --keepdb snpdb",
                        "scripts/upgrade.sh --quick < /dev/null", "pkill runserver", "git pull-request"]:
            self.assertIsNone(self._decide(command), command)

    def test_shared_box_writes_ask(self):
        for command in ["python3 manage.py migrate", "sudo systemctl restart nginx",
                        "python3 manage.py vg page /x --js-writes-to-database",
                        "psql -c 'DELETE FROM snpdb_variant WHERE id = 1'",
                        "cd /opt/variantgrid && git pull", "python3 manage.py test snpdb.tests"]:
            permission, reason = self._decide(command)
            self.assertEqual(permission, "ask", command)
            self.assertTrue(reason.startswith("vg-test2 - this "), reason)

    def test_sanctioned_restart_passes_on_vg_test2_only(self):
        command = "sudo -n systemctl restart gunicorn\nsudo -n systemctl restart celeryd_web_workers"
        self.assertIsNone(self._decide(command, host="vg-test2"))
        self.assertEqual(self._decide(command, host="sharianttest")[0], "ask")
        self.assertEqual(self._decide("sudo -n systemctl restart nginx", host="vg-test2")[0], "ask")

    def test_unguarded_box_skips_asks_but_not_denials(self):
        self.assertIsNone(self._decide("python3 manage.py migrate", unguarded=True))
        self.assertEqual(self._decide("pkill -f runserver", unguarded=True)[0], "deny")

    def test_denied_everywhere(self):
        for command in ["git commit -m 'x\n\nCo-Authored-By: Claude <noreply@anthropic.com>'",
                        "git commit -F msg.txt && echo Claude-Session: https://x",
                        "pkill -f runserver", "pkill -9 -f 'manage.py test'", "pkill -if celery",
                        "gh search issues tso --owner SACGF", "echo hi > CLAUDE.local.md"]:
            self.assertEqual(self._decide(command)[0], "deny", command)
        self.assertIsNone(self._decide("gh search issues tso --repo SACGF/variantgrid"))
