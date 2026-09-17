from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from snpdb.models import AvatarDetails, GlobalSettings, UserAward, UserSettings
from snpdb.models.models_enums import AwardPeriod, UserAwardKind, UserAwardLevel
from snpdb.models.models_user_settings import UserSettingsOverride
from snpdb.user_award_updates import update_badge, update_title, update_user_awards
from snpdb.user_awards import AwardDefinition


class _Counts:
    """ A counter whose results the test sets directly """
    def __init__(self):
        self.counts = {}

    def __call__(self, _since):
        return self.counts


def _title(counter, key="test_title") -> AwardDefinition:
    return AwardDefinition(key=key, kind=UserAwardKind.TITLE, title="Test title", description="", icon="fa-tag",
                           counter=counter, periods=(AwardPeriod.ALL_TIME, AwardPeriod.DAY))


def _badge(counter, key="test_badge", hidden=False) -> AwardDefinition:
    return AwardDefinition(key=key, kind=UserAwardKind.BADGE, title="Test badge", description="", icon="fa-tag",
                           counter=counter, tiers=(5, 10, 20), hidden=hidden)


class UserAwardUpdateTest(TestCase):
    """ Title/badge recompute rules - #1819 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.alice = User.objects.create(username="award_alice")
        cls.bob = User.objects.create(username="award_bob")

    def _titles(self, key="test_title", period=AwardPeriod.ALL_TIME, active=True):
        return set(UserAward.objects.filter(definition_key=key, period=period, active=active).values_list("user_id", flat=True))

    def test_title_changes_hands(self):
        counter = _Counts()
        definition = _title(counter)
        counter.counts = {None: {self.alice.pk: 5, self.bob.pk: 3}}
        update_title(definition, AwardPeriod.ALL_TIME)
        self.assertEqual(self._titles(), {self.alice.pk})
        alice_award = UserAward.objects.get(user=self.alice, definition_key="test_title", period=AwardPeriod.ALL_TIME)
        self.assertEqual(alice_award.count, 5)
        self.assertEqual(alice_award.award_text, "Test title (all time)")

        counter.counts = {None: {self.alice.pk: 5, self.bob.pk: 8}}
        update_title(definition, AwardPeriod.ALL_TIME)
        self.assertEqual(self._titles(), {self.bob.pk})
        self.assertEqual(self._titles(active=False), {self.alice.pk})

        # Regained - the same row flips back, no duplicate
        counter.counts = {None: {self.alice.pk: 9, self.bob.pk: 8}}
        update_title(definition, AwardPeriod.ALL_TIME)
        self.assertEqual(self._titles(), {self.alice.pk})
        self.assertEqual(UserAward.objects.filter(user=self.alice, definition_key="test_title").count(), 1)

    def test_ties_all_hold_and_zero_never_holds(self):
        counter = _Counts()
        counter.counts = {None: {self.alice.pk: 4, self.bob.pk: 4}}
        update_title(_title(counter), AwardPeriod.ALL_TIME)
        self.assertEqual(self._titles(), {self.alice.pk, self.bob.pk})

        counter.counts = {None: {self.alice.pk: 0, self.bob.pk: 0}}
        update_title(_title(counter), AwardPeriod.ALL_TIME)
        self.assertEqual(self._titles(), set())

    def test_title_subjects_are_independent(self):
        counter = _Counts()
        counter.counts = {None: {self.alice.pk: 5, self.bob.pk: 3}, "Pathogenic": {self.bob.pk: 3}}
        update_title(_title(counter), AwardPeriod.DAY)
        bob_award = UserAward.objects.get(user=self.bob, definition_key="test_title", subject="Pathogenic")
        self.assertEqual(bob_award.award_text, "Test title of 'Pathogenic' (today)")
        self.assertEqual(self._titles(period=AwardPeriod.DAY), {self.alice.pk, self.bob.pk})

        # Subject disappears (tag deleted / nobody tagged today) - its holder loses it
        counter.counts = {None: {self.alice.pk: 5}}
        update_title(_title(counter), AwardPeriod.DAY)
        self.assertEqual(self._titles(period=AwardPeriod.DAY), {self.alice.pk})

    def test_badge_tier_only_moves_up(self):
        counter = _Counts()
        definition = _badge(counter)
        counter.counts = {None: {self.alice.pk: 3}}
        update_badge(definition)
        award = UserAward.objects.get(user=self.alice, definition_key="test_badge")
        self.assertFalse(award.active)  # below bronze - progress only
        self.assertEqual(award.count, 3)

        counter.counts = {None: {self.alice.pk: 12}}
        update_badge(definition)
        award.refresh_from_db()
        self.assertTrue(award.active)
        self.assertEqual(award.award_level, UserAwardLevel.SILVER)

        counter.counts = {None: {self.alice.pk: 1}}  # data deleted - never revoked, tier never drops
        update_badge(definition)
        award.refresh_from_db()
        self.assertTrue(award.active)
        self.assertEqual(award.award_level, UserAwardLevel.SILVER)
        self.assertEqual(award.count, 1)

    def test_badge_progress(self):
        counter = _Counts()
        counter.counts = {None: {self.alice.pk: 7}}
        update_badge(_badge(counter))
        progress = AvatarDetails.avatar_for(self.alice).awards.progress(_badge(counter))
        self.assertTrue(progress.earned)
        self.assertEqual(progress.next_threshold, 10)
        self.assertEqual(progress.next_tier_name, "silver")
        self.assertEqual(progress.percent, 70)

        hidden = _badge(counter, key="hidden_badge", hidden=True)
        self.assertFalse(AvatarDetails.avatar_for(self.bob).awards.progress(hidden).visible)

    @override_settings(USER_AWARDS_DISABLED_KEYS={"test_title"})
    def test_disabled_definition_deactivated(self):
        counter = _Counts()
        counter.counts = {None: {self.alice.pk: 5}}
        UserAward.objects.create(user=self.alice, kind=UserAwardKind.TITLE, definition_key="test_title",
                                 period=AwardPeriod.ALL_TIME, award_text="old", active=True)
        update_user_awards([AwardPeriod.ALL_TIME], badges=True, definitions=[_title(counter)])
        self.assertEqual(self._titles(), set())

    @override_settings(USER_AWARDS_ENABLED=False)
    def test_disabled_deployment_short_circuits(self):
        counter = _Counts()
        counter.counts = {None: {self.alice.pk: 5}}
        self.assertFalse(update_user_awards([AwardPeriod.ALL_TIME], badges=True, definitions=[_title(counter)]))
        self.assertFalse(UserAward.objects.filter(definition_key="test_title").exists())


class GridLabelTest(TestCase):
    """ AvatarDetails.grid_label_html - the single decoration rule - #1819 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create(username="label_user", first_name="Evelyn")
        cls.viewer = User.objects.create(username="label_viewer")
        global_settings = GlobalSettings.objects.get()
        global_settings.show_user_awards = True
        global_settings.save()

    def _label(self, viewer=None) -> str:
        return str(AvatarDetails.avatar_for(self.user).grid_label_html(UserSettings.get_for_user(viewer or self.viewer)))

    def _give_title(self):
        UserAward.objects.create(user=self.user, kind=UserAwardKind.TITLE, definition_key="top_tagger",
                                 period=AwardPeriod.ALL_TIME, award_text="Top tagger (all time)")

    def test_plain_without_title(self):
        self.assertEqual(self._label(), "Evelyn")

    def test_decorated_with_title(self):
        self._give_title()
        label = self._label()
        self.assertIn("fa-crown", label)
        self.assertTrue(label.endswith("Evelyn"), label)

    def test_viewer_opts_out(self):
        self._give_title()
        UserSettingsOverride.objects.update_or_create(user=self.viewer, defaults={"show_user_awards": False})
        self.assertEqual(self._label(), "Evelyn")

    @override_settings(USER_AWARDS_ENABLED=False)
    def test_deployment_off(self):
        self._give_title()
        self.assertEqual(self._label(), "Evelyn")

    def test_highest_title_icon_wins(self):
        UserAward.objects.create(user=self.user, kind=UserAwardKind.TITLE, definition_key="top_tagger",
                                 period=AwardPeriod.DAY, award_text="Top tagger (today)")
        UserAward.objects.create(user=self.user, kind=UserAwardKind.TITLE, definition_key="top_tagger",
                                 period=AwardPeriod.MONTH, award_text="Top tagger (this month)")
        avatar = AvatarDetails.avatar_for(self.user)
        self.assertEqual([t.period for t in avatar.titles], [AwardPeriod.MONTH, AwardPeriod.DAY])
        self.assertIn("fa-medal", avatar.title_icon_html)
        self.assertIn("Top tagger (today)", avatar.title_icon_html)
