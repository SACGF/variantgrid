from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import Analysis, VariantTag
from analysis.user_awards import _analyses_worked_on
from annotation.fake_data import create_fake_variants, get_fake_annotation_version
from snpdb.models import GenomeBuild, Tag, Variant


class AnalystAwardCounterTest(TestCase):
    """ The analyst award counter - #1819 """

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

    def test_analyses_worked_on_unions_sources(self):
        """ Created, tagged in, and audit-logged analyses count once each """
        other = Analysis(genome_build=self.grch37)
        other.set_defaults_and_save(User.objects.create(username="other_analyst"))
        # Tagging in someone else's analysis counts it; tagging in our own doesn't double count
        VariantTag.objects.create(genome_build=self.grch37, analysis=other, variant=self.variant,
                                  tag=Tag.objects.get_or_create(pk="artefact")[0], user=self.user)
        self.assertEqual(_analyses_worked_on()[self.user.pk], 2)
