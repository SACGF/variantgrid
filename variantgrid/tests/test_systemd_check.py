"""
The predicate the systemd env files manual migration and deployment_check share: when the installed units
would fail to start because /etc/variantgrid/variantgrid.env or a unit's own .env is missing, or
VG_INSTALL_DIR is another checkout. Also a unit's own .env left over from an older unit (#2482).
"""
import os
import tempfile
from unittest.mock import patch

from django.conf import settings
from django.test import SimpleTestCase

from variantgrid.deployment_validation.systemd_check import (
    SHARED_ENV_FILE,
    check_systemd_units,
    detached_simple_units,
    systemd_env_files_problem,
)

MODULE = "variantgrid.deployment_validation.systemd_check"
GUNICORN_ENV = "/etc/variantgrid/gunicorn/gunicorn.env"


def _units(**env_files):
    return patch(f"{MODULE}._installed_units_env_files", return_value=env_files)


class SystemdEnvFilesProblemTest(SimpleTestCase):

    def test_no_units_installed_has_nothing_to_do(self):
        with _units():
            self.assertIsNone(systemd_env_files_problem())
            self.assertEqual(check_systemd_units(), {})

    def test_missing_env_file_is_named(self):
        with _units(gunicorn=[SHARED_ENV_FILE, GUNICORN_ENV]), patch(f"{MODULE}.os.path.exists", return_value=False):
            self.assertEqual(systemd_env_files_problem(), f"missing {GUNICORN_ENV}, {SHARED_ENV_FILE}")

    def test_install_dir_of_another_checkout(self):
        with _units(gunicorn=[SHARED_ENV_FILE]), patch(f"{MODULE}.os.path.exists", return_value=True), \
                patch(f"{MODULE}.shared_env_install_dir", return_value="/mnt/elsewhere"):
            self.assertIn("/mnt/elsewhere", systemd_env_files_problem())

    def test_this_checkout_is_fine(self):
        with _units(gunicorn=[SHARED_ENV_FILE]), patch(f"{MODULE}.os.path.exists", return_value=True), \
                patch(f"{MODULE}.shared_env_install_dir", return_value=settings.BASE_DIR):
            self.assertIsNone(systemd_env_files_problem())
            self.assertTrue(check_systemd_units()["env files"]["valid"])

    def test_hardcoded_units_warn_but_pass(self):
        with _units(gunicorn=[], celeryd_beat=[]):
            self.assertIsNone(systemd_env_files_problem())
            checks = check_systemd_units()
            self.assertTrue(checks["env files"]["valid"])
            self.assertIn("celeryd_beat, gunicorn", checks["env-driven units"]["warning"])


class DetachedSimpleUnitsTest(SimpleTestCase):

    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp_dir.cleanup)

    def _env_file(self, name: str, celeryd_opts: str) -> str:
        path = os.path.join(self.tmp_dir.name, name)
        with open(path, "w", encoding="utf-8") as f:
            f.write(f'CELERY_BIN=".venv/bin/celery"\nCELERYD_OPTS="{celeryd_opts}"\n')
        return path

    def _installed(self, units: dict[str, tuple[str, str]]):
        """ units: {unit: (Type, its own env file)} """
        def show(unit, *_properties):
            if unit not in units:
                return None
            unit_type, env_file = units[unit]
            return {"LoadState": ["loaded"], "Type": [unit_type],
                    "EnvironmentFiles": [f"{SHARED_ENV_FILE} (ignore_errors=yes)", f"{env_file} (ignore_errors=yes)"]}
        return patch(f"{MODULE}.systemctl_show", side_effect=show)

    def test_old_beat_env_under_simple_unit(self):
        old_beat_env = self._env_file("celeryd_beat.env", "--detach --pidfile=/var/run/beat.pid --schedule=x")
        worker_env = self._env_file("celeryd_db_workers.env", "--detach -n db_workers")
        with self._installed({"celeryd_beat": ("simple", old_beat_env),
                              "celeryd_db_workers": ("forking", worker_env)}):
            self.assertEqual({"celeryd_beat": old_beat_env}, detached_simple_units())

    def test_current_beat_env_is_fine(self):
        beat_env = self._env_file("celeryd_beat.env", "--schedule=data/celery_beat_schedule")
        with self._installed({"celeryd_beat": ("simple", beat_env)}):
            self.assertEqual({}, detached_simple_units())
