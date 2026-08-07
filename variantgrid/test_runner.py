from django.test.runner import DiscoverRunner

from genes.tests.utils.mock_transcript_sequence_retrieval import MockTranscriptSequenceFetcher
from genes.transcript_sequence_retrieval import TranscriptSequenceFetcher
from snpdb.clingen_allele_api import ClinGenAlleleRegistryAPI
from snpdb.tests.utils.mock_clingen_api import MockClinGenAlleleRegistryAPI


class VariantGridTestRunner(DiscoverRunner):
    """ Points the external service clients at implementations serving recorded data, so tests neither
        depend on those services being up nor take the latency of calling them.

        Each mock raises when asked for something its recordings don't cover, naming the fixture to add. """

    def setup_test_environment(self, **kwargs):
        super().setup_test_environment(**kwargs)
        ClinGenAlleleRegistryAPI.override_class = MockClinGenAlleleRegistryAPI
        TranscriptSequenceFetcher.override_class = MockTranscriptSequenceFetcher
