""" NotificationBuilder's Slack markdown -> email HTML conversion """
from django.test import SimpleTestCase

from eventlog.notifications import NotificationBuilder


class SlackMarkdownToHtmlTest(SimpleTestCase):
    def test_emoji_stripped_and_bold_converted(self):
        self.assertEqual(NotificationBuilder.slack_markdown_to_html(":hospital: This is **bold**"),
                         "This is <strong>bold</strong>")
