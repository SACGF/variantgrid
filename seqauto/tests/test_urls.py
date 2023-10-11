import unittest

from django.contrib.auth.models import User

from library.django_utils.unittest_utils import URLTestCase
from seqauto.models import QCColumn, EnrichmentKit, SequencingRun, SequencerModel, DataGeneration, Sequencer
from snpdb.models import Manufacturer, DataState


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
                                                                 sequencer=sequencer,
                                                                 data_state=DataState.COMPLETE)[0]

        cls.enrichment_kit = EnrichmentKit.objects.get_or_create(name="fake_kit", version=1, manufacturer=illumina)[0]
        cls.qc_column = QCColumn.objects.first()

    def testDataGridUrls(self):
        """ Grids w/o permissions """

        GRID_LIST_URLS = [
            ("experiments_datatable", {}, 200),
            ("enrichment_kit_datatable", {}, 200),
        ]
        self._test_urls(GRID_LIST_URLS, self.user_non_owner)

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
