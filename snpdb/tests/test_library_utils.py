from django.test import TestCase

from library.log_utils import NotificationBuilder
from library.utils import format_significant_digits


class TestUtils(TestCase):
    """ Test for library.utils - Django unit tests are run in apps so this needs to be here """
    def test_sig_digits(self):
        self.assertEqual("0", format_significant_digits(0))
        self.assertEqual("1", format_significant_digits(1))
        self.assertEqual("10000", format_significant_digits(10000))
        self.assertEqual("1.23", format_significant_digits(1.234567))
        self.assertEqual("-1.23", format_significant_digits(-1.234567))
        self.assertEqual("456", format_significant_digits(456.12))
        self.assertEqual("1.1", format_significant_digits(1.10004))
        self.assertEqual("1.11", format_significant_digits(1.114))
        self.assertEqual("1.12", format_significant_digits(1.116))
        self.assertEqual("0.0000015", format_significant_digits(0.00000150002))
        self.assertEqual("-0.0000015", format_significant_digits(-0.00000150002))

    def test_markdown(self):
        markdown_text = ":hospital: This is **bold**"
        html = NotificationBuilder.slack_markdown_to_html(markdown_text)
        self.assertEqual("This is <strong>bold</strong>", html)
