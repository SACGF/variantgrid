from unittest.mock import patch

from django.test import SimpleTestCase, override_settings

from library.log_utils import report_message


@patch("library.log_utils.rollbar.report_message")
class RollbarMinLevelTest(SimpleTestCase):

    @override_settings(ROLLBAR={"min_level": "warning"})
    def test_min_level(self, mock_report_message):
        report_message("below", level="info")
        mock_report_message.assert_not_called()

        for level in ["warning", "error"]:
            report_message("at or above", level=level)
        self.assertEqual(mock_report_message.call_count, 2)

    @override_settings(ROLLBAR={})
    def test_no_min_level_sends_everything(self, mock_report_message):
        report_message("info", level="info")
        mock_report_message.assert_called_once()
