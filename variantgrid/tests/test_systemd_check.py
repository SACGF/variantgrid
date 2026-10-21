"""
The predicate the systemd env files manual migration and deployment_check share: when the installed units
would fail to start because /etc/variantgrid/variantgrid.env or a unit's own .env is missing, or
VG_INSTALL_DIR is another checkout.
"""
from unittest.mock import patch

from django.conf import settings
from django.test import SimpleTestCase

from variantgrid.deployment_validation.systemd_check import (
    SHARED_ENV_FILE,
    check_systemd_units,
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
