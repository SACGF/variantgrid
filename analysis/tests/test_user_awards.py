from datetime import datetime, timezone as dt_timezone

from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import Analysis, VariantTag
from analysis.user_awards import _tags_in_hours, _tags_per_tag
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from snpdb.models import GenomeBuild, Tag, Variant
from snpdb.models.models_user_settings import UserSettingsOverride


class TagAwardCountersTest(TestCase):
    """ The tag based award counters - #1819 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create(username="award_tagger")
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        create_fake_variants(cls.grch37)
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)
        cls.variant = Variant.objects.filter(Variant.get_no_reference_q()).order_by("pk").first()

    def _tag(self, name: str, created: datetime) -> VariantTag:
        tag = Tag.objects.get_or_create(pk=name)[0]
        variant_tag = VariantTag.objects.create(genome_build=self.grch37, analysis=self.analysis,
                                                variant=self.variant, tag=tag, user=self.user)
        VariantTag.objects.filter(pk=variant_tag.pk).update(created=created)
        return variant_tag

    def test_top_tagger_per_tag_and_overall(self):
        when = datetime(2026, 3, 1, 12, 0, tzinfo=dt_timezone.utc)
        self._tag("artefact", when)
        self._tag("artefact", when)
        self._tag("pathogenic", when)
        self.assertEqual(_tags_per_tag(None), {
            "artefact": {self.user.pk: 2},
            "pathogenic": {self.user.pk: 1},
            None: {self.user.pk: 3},
        })
        self.assertEqual(_tags_per_tag(datetime(2026, 4, 1, tzinfo=dt_timezone.utc)), {})

    def test_night_owl_uses_the_users_timezone(self):
        night_owl = _tags_in_hours({22, 23, 0, 1, 2, 3, 4})
        self._tag("artefact", datetime(2026, 3, 1, 23, 30, tzinfo=dt_timezone.utc))

        UserSettingsOverride.objects.update_or_create(user=self.user, defaults={"timezone": "UTC"})
        self.assertEqual(night_owl(None), {None: {self.user.pk: 1}})

        # 23:30 UTC is 10:00 the next morning in Adelaide
        UserSettingsOverride.objects.update_or_create(user=self.user, defaults={"timezone": "Australia/Adelaide"})
        self.assertEqual(night_owl(None), {None: {}})
