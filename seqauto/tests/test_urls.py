import unittest

from django.contrib.auth.models import User

from library.django_utils.unittest_utils import URLTestCase
from seqauto.models import (
    DataGeneration,
    EnrichmentKit,
    QCColumn,
    Sequencer,
    SequencerModel,
    SequencingRun,
)
from snpdb.models import Manufacturer


class Test(URLTestCase):
    enrichment_kit = None
    qc_column = None
    sequencing_run = None
    user_owner = None

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user_non_owner = User.objects.get_or_create(username='different_user')[0]

        illumina = Manufacturer.objects.get_or_create(name="Illumina")[0]
        sequencer_model = SequencerModel.objects.get_or_create(model='NextSeq 500',
                                                               data_naming_convention=DataGeneration.MISEQ,
                                                               manufacturer=illumina)[0]
        sequencer = Sequencer.objects.get_or_create(name='NB501008', sequencer_model=sequencer_model)[0]
        cls.sequencing_run = SequencingRun.objects.get_or_create(name="200626_NB501009_0391_AHFHLJBGXG",
                                                                 sequencer=sequencer)[0]

        cls.enrichment_kit = EnrichmentKit.objects.get_or_create(name="fake_kit", version=1, manufacturer=illumina)[0]
        cls.qc_column = QCColumn.objects.first()

    def testDataGridUrls(self):
        """ Grids w/o permissions """

        GRID_LIST_URLS = [
            ("experiments_datatable", {}, 200),
            ("sequencing_run_datatable", {}, 200),
            ("unaligned_reads_datatable", {}, 200),
            ("alignment_file_datatable", {}, 200),
            ("vcf_file_datatable", {}, 200),
            ("qc_datatable", {}, 200),
            ("enrichment_kit_datatable", {}, 200),
            ("illumina_flowcell_qc_datatable", {}, 200),
            ("fastqc_datatable", {}, 200),
            ("flagstats_datatable", {}, 200),
            ("qc_exec_summary_datatable", {}, 200),
            ("sequencing_samples_datatable", {}, 200),
            ("sequencing_samples_historical_datatable", {"time_frame": "year"}, 200),
        ]
        self._test_datatable_urls(GRID_LIST_URLS, self.user_non_owner)

    def testSoftwareVersionsStaffOnly(self):
        staff_user = User.objects.get_or_create(username='staff_user', is_staff=True)[0]
        names = ["library_datatable", "sequencer_datatable", "assay_datatable", "aligner_datatable",
                 "variant_caller_datatable", "variant_calling_pipeline_datatable"]
        self._test_urls([("sequencing_software_versions", {}, 200)], staff_user)
        self._test_datatable_urls([(name, {}, 200) for name in names], staff_user)
        self._test_urls([("sequencing_software_versions", {}, 302)] + [(name, {}, 302) for name in names],
                        self.user_non_owner)

    def testAutocompleteUrls(self):
        # panel_app_forward = json.dumps({"server_id": self.panel_app_panel.server_id})
        AUTOCOMPLETE_URLS = [
            ('qc_column_autocomplete', self.qc_column, {"q": self.qc_column.name}),
            ('enrichment_kit_autocomplete', self.enrichment_kit, {"q": self.enrichment_kit.name}),
            ('sequencing_run_autocomplete', self.sequencing_run, {"q": self.sequencing_run.name}),
        ]
        self._test_autocomplete_urls(AUTOCOMPLETE_URLS, self.user_non_owner, True)


if __name__ == "__main__":
    unittest.main()
