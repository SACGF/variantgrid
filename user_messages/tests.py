from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from user_messages.models import Message, inbox_count_for


class MessageModelTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.alice = User.objects.create_user("alice", password="x")
        cls.bob = User.objects.create_user("bob", password="x")

    def test_body_html_renders_markdown_link(self):
        msg = Message.objects.create(subject="s", body="see [vcf #1](/view/1)", sender=self.alice, recipient=self.bob)
        self.assertIn('<a href="/view/1">vcf #1</a>', msg.body_html)

    def test_body_html_neutralizes_script(self):
        msg = Message.objects.create(subject="s", body="<script>alert('x')</script>", sender=self.alice, recipient=self.bob)
        html = msg.body_html
        self.assertNotIn("<script>", html)
        self.assertIn("&lt;script&gt;", html)

    def test_body_html_strips_javascript_href(self):
        msg = Message.objects.create(subject="s", body="[click](javascript:alert(1))", sender=self.alice, recipient=self.bob)
        self.assertNotIn("javascript:", msg.body_html)

    def test_sent_at_stamped_on_create(self):
        msg = Message.objects.create(subject="s", body="b", sender=self.alice, recipient=self.bob)
        self.assertIsNotNone(msg.sent_at)

    def test_manager_partitions_and_count(self):
        Message.objects.create(subject="to bob", body="b", sender=self.alice, recipient=self.bob)
        self.assertEqual(Message.objects.inbox_for(self.bob).count(), 1)
        self.assertEqual(Message.objects.outbox_for(self.alice).count(), 1)
        self.assertEqual(Message.objects.inbox_for(self.alice).count(), 0)
        self.assertEqual(inbox_count_for(self.bob), 1)


class MessageViewTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.alice = User.objects.create_user("alice", password="x")
        cls.bob = User.objects.create_user("bob", password="x")
        cls.carol = User.objects.create_user("carol", password="x")

    def setUp(self):
        self.client.force_login(self.bob)

    def _msg(self):
        return Message.objects.create(subject="hi", body="body", sender=self.alice, recipient=self.bob)

    def test_view_marks_read(self):
        msg = self._msg()
        self.assertTrue(msg.new())
        self.client.get(reverse('messages_detail', args=[msg.id]))
        msg.refresh_from_db()
        self.assertIsNotNone(msg.read_at)

    def test_view_non_party_404(self):
        msg = self._msg()
        self.client.force_login(self.carol)
        response = self.client.get(reverse('messages_detail', args=[msg.id]))
        self.assertEqual(response.status_code, 404)

    def test_compose_creates_message(self):
        response = self.client.post(reverse('messages_compose'),
                                    {"recipient": "alice", "subject": "hello", "body": "hi there"})
        self.assertEqual(response.status_code, 302)
        self.assertTrue(Message.objects.filter(sender=self.bob, recipient=self.alice, subject="hello").exists())

    def test_delete_marks_recipient_deleted(self):
        msg = self._msg()
        self.client.get(reverse('messages_delete', args=[msg.id]))
        msg.refresh_from_db()
        self.assertIsNotNone(msg.recipient_deleted_at)
        self.assertEqual(Message.objects.trash_for(self.bob).count(), 1)

    def test_list_and_compose_pages_render(self):
        self._msg()
        for name in ['messages_inbox', 'messages_outbox', 'messages_trash', 'messages_compose']:
            response = self.client.get(reverse(name))
            self.assertEqual(response.status_code, 200, name)
