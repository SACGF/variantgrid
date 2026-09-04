from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.urls.base import resolve, reverse
from django.utils.timezone import now

from analysis.grids import AnalysesListColumns
from analysis.models import Analysis
from seqauto.grids.sequencing_data_grids import SequencingRunColumns
from seqauto.models import DataGeneration, Sequencer, SequencerModel, SequencingRun
from snpdb.grids import SamplesListColumns, VCFListColumns
from snpdb.models import VCF, GenomeBuild, ImportStatus, Sample
from snpdb.views.datatable_view import DatabaseTableView


class GridSearchTests(TestCase):
    """ The list grids' search boxes OR icontains over every searchable column, so a numeric/date
        column left in the search set is a 500 on the first thing anyone types """

    def setUp(self):
        self.user = User.objects.create(username='grid_search_user', is_superuser=True)
        genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        self.vcf = VCF.objects.create(name="findme_vcf", genome_build=genome_build, user=self.user,
                                      date=now(), genotype_samples=1, import_status=ImportStatus.SUCCESS)
        self.sample = Sample.objects.create(name="findme_sample", vcf=self.vcf,
                                            import_status=ImportStatus.SUCCESS)
        self.analysis = Analysis.objects.create(name="findme_analysis", user=self.user,
                                                genome_build=genome_build, visible=True)
        sequencer_model, _ = SequencerModel.objects.get_or_create(model="MiSeq",
                                                                  data_naming_convention=DataGeneration.MISEQ)
        sequencer = Sequencer.objects.create(name="M0_search", sequencer_model=sequencer_model)
        self.sequencing_run = SequencingRun.objects.create(name="findme_run", sequencer=sequencer)

    def _rows(self, url_name, column_class, search):
        url = reverse(url_name)
        request = RequestFactory().get(url, {"search[value]": search})
        request.resolver_match = resolve(url)
        request.user = self.user
        view = DatabaseTableView(column_class=column_class)
        view.request = request
        view.config = column_class(request)
        qs = view.config.get_initial_queryset()
        qs = view.config.filter_queryset(qs)
        qs = view.config.power_search(qs, search)
        return view.prepare_results(view.config.ordering(qs))

    def test_vcf(self):
        self.assertEqual(len(self._rows('vcfs_datatable', VCFListColumns, "findme")), 1)
        self.assertEqual(len(self._rows('vcfs_datatable', VCFListColumns, "nope")), 0)

    def test_sample(self):
        self.assertEqual(len(self._rows('samples_list_datatable', SamplesListColumns, "findme")), 1)
        self.assertEqual(len(self._rows('samples_list_datatable', SamplesListColumns, "nope")), 0)

    def test_analysis(self):
        self.assertEqual(len(self._rows('analyses_list_datatable', AnalysesListColumns, "findme")), 1)
        self.assertEqual(len(self._rows('analyses_list_datatable', AnalysesListColumns, "nope")), 0)

    def test_sequencing_run(self):
        self.assertEqual(len(self._rows('sequencing_run_datatable', SequencingRunColumns, "findme")), 1)
        self.assertEqual(len(self._rows('sequencing_run_datatable', SequencingRunColumns, "nope")), 0)

    def test_search_by_pk(self):
        """ IDs are written down and typed back in - they're integers, so they can't go through the
            icontains the text columns use """
        for url_name, column_class, obj in [('vcfs_datatable', VCFListColumns, self.vcf),
                                            ('samples_list_datatable', SamplesListColumns, self.sample),
                                            ('analyses_list_datatable', AnalysesListColumns, self.analysis)]:
            with self.subTest(url_name):
                self.assertEqual(len(self._rows(url_name, column_class, str(obj.pk))), 1)
                self.assertEqual(len(self._rows(url_name, column_class, str(obj.pk + 5000))), 0)
                # A number too wide for the column is a miss, not a database error
                self.assertEqual(len(self._rows(url_name, column_class, "9" * 30)), 0)
