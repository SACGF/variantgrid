import os
from unittest import TestCase
from unittest.mock import patch

from variantgrid.settings.components import secret_settings
from variantgrid.settings.components.secret_settings import (
    _get_env_variable,
    _settings_file,
    get_secret,
)


class GetEnvVariableTest(TestCase):
    """ Plain unittest: the module runs at settings load, so it has to work without Django """

    def test_missing_and_empty_are_both_unset(self):
        with patch.dict(os.environ, {"VG_TEST_SECRET": ""}):
            self.assertEqual((None, False), _get_env_variable("VG_TEST_SECRET"))
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual((None, False), _get_env_variable("VG_TEST_SECRET"))

    def test_conversions(self):
        with patch.dict(os.environ, {"A": "true", "B": "false", "C": '"quoted"', "D": "plain"}):
            self.assertEqual((True, True), _get_env_variable("A"))
            self.assertEqual((False, True), _get_env_variable("B"))
            self.assertEqual(("quoted", True), _get_env_variable("C"))
            self.assertEqual(("plain", True), _get_env_variable("D"))

    def test_empty_settings_config_uses_default_path(self):
        with patch.dict(os.environ, {"SETTINGS_CONFIG": ""}):
            self.assertEqual("/etc/variantgrid/settings_config.json", _settings_file())


class GetSecretTest(TestCase):
    def test_empty_env_var_falls_through_to_file_then_default(self):
        with patch.dict(os.environ, {"DB.host": ""}), \
                patch.object(secret_settings, "_settings_json", {"DB": {"host": "from-file"}}):
            self.assertEqual("from-file", get_secret("DB.host"))
        with patch.dict(os.environ, {"DB.host": ""}), \
                patch.object(secret_settings, "_settings_json", {}):
            with self.assertLogs(level="WARNING") as logs:
                self.assertEqual("localhost", get_secret("DB.host"))
            self.assertIn("using the built-in default", logs.output[0])
            self.assertNotIn("localhost", logs.output[0])

    def test_env_var_wins(self):
        with patch.dict(os.environ, {"DB.host": "from-env"}), \
                patch.object(secret_settings, "_settings_json", {"DB": {"host": "from-file"}}):
            self.assertEqual("from-env", get_secret("DB.host"))
