from unittest.mock import PropertyMock, patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone

from classification.enums.overlaps_enums import ClassificationResultValue, OverlapType
from classification.models.overlaps_model import Overlap
from classification.tests.models.test_utils import ClassificationTestUtils
from review.models import Review, ReviewTopic


class OverlapReviewPermissionTestCase(TestCase):
    """ Only members of an Overlap's reviewing labs can start, edit or action a Review of it """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.lab_user = ClassificationTestUtils.lab_and_user()
        self.outsider = User.objects.get_or_create(username="overlap_review_outsider")[0]
        self.topic = ReviewTopic.objects.get_or_create(key="discordance_report", defaults={"name": "Discordance"})[0]
        self.overlap = Overlap.objects.create(overlap_type=OverlapType.SINGLE_CONTEXT,
                                              value_type=ClassificationResultValue.ONC_PATH)
        self.reviewed_object = self.overlap.reviews_safe
        reviewing_labs = patch.object(Overlap, "reviewing_labs", new_callable=PropertyMock,
                                      return_value={self.lab})
        reviewing_labs.start()
        self.addCleanup(reviewing_labs.stop)

    def _review(self) -> Review:
        return Review.objects.create(reviewing=self.reviewed_object, topic=self.topic, user=self.lab_user,
                                     review_date=timezone.now().date(), meeting_meta={})

    def test_new_review_requires_reviewing_lab_membership(self):
        url = reverse("start_review", kwargs={"reviewed_object_id": self.reviewed_object.pk,
                                              "topic_id": self.topic.pk})
        self.client.force_login(self.outsider)
        self.assertEqual(self.client.post(url, {}).status_code, 403)
        self.assertFalse(Review.objects.filter(reviewing=self.reviewed_object).exists())

        self.client.force_login(self.lab_user)
        self.assertEqual(self.client.get(url).status_code, 200)

    def test_edit_review_is_read_only_for_non_members(self):
        url = reverse("edit_review", kwargs={"review_id": self._review().pk})
        self.client.force_login(self.outsider)
        self.assertTemplateUsed(self.client.get(url), "review/review_detail.html")

        self.client.force_login(self.lab_user)
        self.assertTemplateUsed(self.client.get(url), "review/review.html")

    def test_overlap_review_action_requires_reviewing_lab_membership(self):
        review = self._review()
        self.client.force_login(self.outsider)
        response = self.client.post(reverse("action_overlap_review", kwargs={"review_id": review.pk}),
                                    {"action": "postpone"})
        self.assertEqual(response.status_code, 403)
        review.refresh_from_db()
        self.assertFalse(review.is_complete)
