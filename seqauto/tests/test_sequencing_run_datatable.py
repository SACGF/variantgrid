from django.contrib.auth.models import User
from django.test import RequestFactory, TestCase
from django.urls.base import resolve, reverse
from django.utils.timezone import now

from seqauto.grids.sequencing_data_grids import SequencingRunColumns
from seqauto.models import (
    DataGeneration,
    EnrichmentKit,
    SampleSheet,
    Sequencer,
    SequencerModel,
    SequencingRun,
    SequencingRunCurrentSampleSheet,
    SequencingSample,
    VariantCaller,
    VCFFromSequencingRun,
)
from snpdb.models import VCF, GenomeBuild, ImportStatus, UserGridConfig
from snpdb.views.datatable_view import prepare_rows


class SequencingRunDatatableTests(TestCase):
    """ The row aggregates join through both sequencing samples and VCFs, so each multiplies the other out """

    def setUp(self):
        self.user = User.objects.create(username='sequencing_run_datatable_user')
        sequencer_model, _ = SequencerModel.objects.get_or_create(model="MiSeq",
                                                                  data_naming_convention=DataGeneration.MISEQ)
        sequencer = Sequencer.objects.create(name="M02027_datatable", sequencer_model=sequencer_model)
        self.sequencing_run = SequencingRun.objects.create(name="RUN_DT", sequencer=sequencer)
        sample_sheet = SampleSheet.objects.create(sequencing_run=self.sequencing_run,
                                                  path="/data/RUN_DT/SampleSheet.csv", hash="HASH_DT")
        SequencingRunCurrentSampleSheet.objects.create(sequencing_run=self.sequencing_run,
                                                       sample_sheet=sample_sheet)
        for i, sample_name in enumerate(["s1", "s2", "s3"], start=1):
            SequencingSample.objects.create(sample_sheet=sample_sheet, sample_id=sample_name,
                                            sample_name=sample_name, sample_number=i, barcode="ACGT")

        genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        for i, (caller_name, import_status) in enumerate([("gatk", ImportStatus.SUCCESS),
                                                          ("freebayes", ImportStatus.ERROR)]):
            variant_caller = VariantCaller.objects.create(name=caller_name, version="1")
            vcf = VCF.objects.create(name=f"vcf_{i}", genome_build=genome_build, user=self.user,
                                     date=now(), genotype_samples=0, import_status=import_status)
            VCFFromSequencingRun.objects.create(vcf=vcf, sequencing_run=self.sequencing_run,
                                                variant_caller=variant_caller)

    def _rows(self, **params) -> list[dict]:
        url = reverse('sequencing_run_datatable')
        request = RequestFactory().get(url, params)
        request.resolver_match = resolve(url)
        request.user = self.user
        config = SequencingRunColumns(request)
        qs = config.filter_queryset(config.get_initial_queryset())
        return prepare_rows(config, config.ordering(qs))

    def test_sample_count_and_vcfs(self):
        row = self._rows()[0]
        self.assertEqual(row["sample_count"], 3)
        self.assertEqual(row["vcf_ids"], [
            {"id": str(v.pk), "url": v.get_absolute_url(), "variant_caller": caller, "import_status": status}
            for v, caller, status in [
                (VCF.objects.get(name="vcf_0"), "gatk", ImportStatus.SUCCESS),
                (VCF.objects.get(name="vcf_1"), "freebayes", ImportStatus.ERROR),
            ]])

    def test_hidden_filter(self):
        SequencingRun.objects.filter(pk=self.sequencing_run.pk).update(hidden=True)
        self.assertEqual(self._rows(), [])

        user_grid_config = UserGridConfig.get(self.user, 'SequencingRuns')
        user_grid_config.show_hidden_data = True
        user_grid_config.save()
        self.assertEqual(len(self._rows()), 1)

    def test_enrichment_kit_filter(self):
        enrichment_kit = EnrichmentKit.objects.create(name="datatable_kit", version=1)
        self.assertEqual(self._rows(enrichment_kit_id=enrichment_kit.pk), [])

        SequencingRun.objects.filter(pk=self.sequencing_run.pk).update(enrichment_kit=enrichment_kit)
        self.assertEqual(len(self._rows(enrichment_kit_id=enrichment_kit.pk)), 1)
