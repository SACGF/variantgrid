"""
'manage.py create_fake_data all' on the test database: the steps find each other, the pages they are for render,
a second run adds nothing, and '--delete' takes it all away again.
"""
from io import StringIO

from django.contrib.auth.models import User
from django.core.management import call_command

from analysis.fake_data import FAKE_ANALYSIS_NAME
from analysis.models import Analysis, VariantTag
from annotation.fake_data import FAKE_VARIANTS_PIPELINE_COMMAND
from annotation.models import VariantAnnotation
from classification.models import Classification
from library.django_utils.unittest_utils import URLTestCase
from snpdb.fake_data import FAKE_ORGANIZATION_GROUP, FAKE_USERS
from snpdb.models import Lab, Trio

SMALL = {
    "variants": 60,
    "max_genotypes": 60,
    "somatic_events": 100,
    "germline_events": 100,
    "events_per_variant": 2,
    "alleles": 12,
    "classifications": 20,
    "somatic_classifications": 10,
}


def fake_row_counts() -> dict[str, int]:
    return {
        "users": User.objects.filter(username__in=FAKE_USERS).count(),
        "labs": Lab.objects.filter(group_name__startswith=FAKE_ORGANIZATION_GROUP).count(),
        "variant annotations": VariantAnnotation.objects.filter(
            annotation_run__pipeline_command=FAKE_VARIANTS_PIPELINE_COMMAND).count(),
        "trios": Trio.objects.filter(user__username__in=FAKE_USERS).count(),
        "analyses": Analysis.objects.filter(name=FAKE_ANALYSIS_NAME).count(),
        "variant tags": VariantTag.objects.filter(tag__pk__startswith="fake-").count(),
        "classifications": Classification.objects.filter(lab_record_id__startswith="fake-").count(),
    }


class CreateFakeDataTest(URLTestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        call_command("create_fake_data", "all", genome_build="GRCh37", stdout=StringIO(), **SMALL)
        cls.admin = User.objects.create_superuser("create_fake_data_admin", "admin@fake.com", "password")

    def test_all_creates_every_step(self):
        for name, count in fake_row_counts().items():
            self.assertGreater(count, 0, f"No fake {name}")

    def test_pages_render(self):
        trio = Trio.objects.get(user__username__in=FAKE_USERS)
        analysis = Analysis.objects.get(name=FAKE_ANALYSIS_NAME)
        allele = Classification.objects.filter(lab_record_id__startswith="fake-class-").first().allele
        self._test_urls([
            ("analysis", {"analysis_id": analysis.pk}, 200),
            ("view_trio", {"pk": trio.pk}, 200),
        ], trio.user)
        self._test_urls([
            ("classifications", {}, 200),
            ("view_allele", {"allele_id": allele.pk}, 200),
            ("tag_stats", {}, 200),
            ("classification_reclassification_analytics", {}, 200),
        ], self.admin)

    def test_running_again_creates_nothing(self):
        counts = fake_row_counts()
        call_command("create_fake_data", "all", genome_build="GRCh37", stdout=StringIO(), **SMALL)
        self.assertEqual(fake_row_counts(), counts)

    def test_delete_leaves_no_fake_rows(self):
        call_command("create_fake_data", "all", "--delete", genome_build="GRCh37", stdout=StringIO())
        self.assertEqual(set(fake_row_counts().values()), {0})
