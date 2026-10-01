from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from snpdb.models import AvatarDetails, UserAward
from snpdb.models.models_enums import UserAwardKind, UserAwardLevel
from snpdb.user_award_updates import update_badge, update_user_awards
from snpdb.user_awards import AwardDefinition


class _Counts:
    """ A counter whose results the test sets directly """
    def __init__(self):
        self.counts = {}

    def __call__(self):
        return self.counts


def _badge(counter, key="test_badge") -> AwardDefinition:
    return AwardDefinition(key=key, title="Test badge", description="", icon="fa-tag",
                           counter=counter, tiers=(5, 10, 20))


class UserAwardUpdateTest(TestCase):
    """ Badge recompute rules - #1819 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.alice = User.objects.create(username="award_alice")

    def test_badge_tier_only_moves_up(self):
        counter = _Counts()
        definition = _badge(counter)
        counter.counts = {self.alice.pk: 3}
        update_badge(definition)
        award = UserAward.objects.get(user=self.alice, definition_key="test_badge")
        self.assertFalse(award.active)  # below bronze - progress only
        self.assertEqual(award.count, 3)

        counter.counts = {self.alice.pk: 12}
        update_badge(definition)
        award.refresh_from_db()
        self.assertTrue(award.active)
        self.assertEqual(award.award_level, UserAwardLevel.SILVER)

        counter.counts = {self.alice.pk: 1}  # data deleted - never revoked, tier never drops
        update_badge(definition)
        award.refresh_from_db()
        self.assertTrue(award.active)
        self.assertEqual(award.award_level, UserAwardLevel.SILVER)
        self.assertEqual(award.count, 1)

    def test_badge_progress(self):
        counter = _Counts()
        counter.counts = {self.alice.pk: 7}
        update_badge(_badge(counter))
        progress = AvatarDetails.avatar_for(self.alice).awards.progress(_badge(counter))
        self.assertTrue(progress.earned)
        self.assertEqual(progress.next_threshold, 10)
        self.assertEqual(progress.next_tier_name, "silver")
        self.assertEqual(progress.percent, 70)

    @override_settings(USER_AWARDS_DISABLED_KEYS={"test_badge"})
    def test_disabled_definition_deactivated(self):
        counter = _Counts()
        counter.counts = {self.alice.pk: 50}
        UserAward.objects.create(user=self.alice, kind=UserAwardKind.BADGE, definition_key="test_badge",
                                 award_text="old", active=True)
        update_user_awards(definitions=[_badge(counter)])
        award = UserAward.objects.get(user=self.alice, definition_key="test_badge")
        self.assertFalse(award.active)

    @override_settings(USER_AWARDS_ENABLED=False)
    def test_disabled_deployment_short_circuits(self):
        counter = _Counts()
        counter.counts = {self.alice.pk: 5}
        self.assertFalse(update_user_awards(definitions=[_badge(counter)]))
        self.assertFalse(UserAward.objects.filter(definition_key="test_badge").exists())
