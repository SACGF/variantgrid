from django.contrib.auth.models import User

from analysis.models import Analysis
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild


class AnalysisSetupMixin:
    """setUpTestData mixin that provides cls.analysis and cls.grch37.

    Shared by multiple analysis test classes to avoid repeating the same
    user + genome build + annotation version + analysis setup boilerplate.
    """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username=f"test_{cls.__name__}")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)
