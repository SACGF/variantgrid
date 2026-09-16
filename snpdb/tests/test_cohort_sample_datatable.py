from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.urls.base import resolve, reverse

from snpdb.grids import CohortSampleListColumns
from snpdb.models import Cohort, GenomeBuild, ImportStatus, Sample
from snpdb.tests.utils.fake_cohort_data import create_fake_trio


class CohortSampleDatatableTests(TestCase):
    """ The cohort page runs the one endpoint twice - the samples it holds, and the ones it could add """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="cohort_sample_datatable_user")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.trio = create_fake_trio(cls.user, cls.grch37)
        cls.samples = list(cls.trio.get_samples())
        Sample.objects.filter(pk__in=[s.pk for s in cls.samples]).update(import_status=ImportStatus.SUCCESS)
        cls.cohort = Cohort.objects.create(name="picker_cohort", user=cls.user, genome_build=cls.grch37)
        cls.cohort.add_sample(cls.samples[0].pk)

    def _sample_names(self, **params) -> set[str]:
        url = reverse('cohort_sample_datatable', kwargs={"cohort_id": self.cohort.pk})
        request = RequestFactory().get(url, params)
        request.resolver_match = resolve(url)
        request.user = self.user
        config = CohortSampleListColumns(request)
        qs = config.filter_queryset(config.get_initial_queryset())
        return set(qs.values_list("name", flat=True))

    def test_show_cohort_samples(self):
        self.assertEqual(self._sample_names(cohort_op="show_cohort"), {self.samples[0].name})

    def test_available_samples_exclude_cohort(self):
        available = self._sample_names(cohort_op="exclude_cohort")
        self.assertNotIn(self.samples[0].name, available)
        self.assertIn(self.samples[1].name, available)

    def test_defaults_to_available_samples(self):
        self.assertEqual(self._sample_names(), self._sample_names(cohort_op="exclude_cohort"))

    def test_unknown_cohort_op(self):
        with self.assertRaises(ValueError):
            self._sample_names(cohort_op="nonsense")
