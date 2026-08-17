import atexit
import json
import os

import cdot.hgvs.dataproviders.fasta_seqfetcher as fasta_seqfetcher
from django.test.runner import DiscoverRunner

import library.genomics.fasta_wrapper as fasta_wrapper
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


class FastaRecordingRunner(VariantGridTestRunner):
    """ Records every genome-fasta region the suite fetches, for regenerating the sparse test
        fastas CI runs against (see variantgrid/data/reference/sparse_test_fastas/README.md).

        Run on a machine with the real reference fastas:
            python3 manage.py test --keepdb --testrunner=variantgrid.test_runner.FastaRecordingRunner
        then feed the regions file to scripts/generate_sparse_test_fastas.py """

    REGIONS_FILE = os.environ.get("VG_FASTA_REGIONS_FILE", "/tmp/vg_fasta_regions.jsonl")

    def setup_test_environment(self, **kwargs):
        super().setup_test_environment(**kwargs)
        regions = []
        atexit.register(self._dump, regions)

        def _log(filename, contig, start, end):
            if isinstance(filename, bytes):
                filename = filename.decode()
            regions.append((str(filename), str(contig), start, end))

        orig_fetch_seq = fasta_seqfetcher.GenomeFastaSeqFetcher.fetch_seq

        def fetch_seq(seqfetcher, ac, start_i=None, end_i=None):
            result = orig_fetch_seq(seqfetcher, ac, start_i=start_i, end_i=end_i)
            if fasta_file := seqfetcher.contig_fastas.get(ac):
                _log(fasta_file.filename, ac, start_i, end_i)
            return result

        fasta_seqfetcher.GenomeFastaSeqFetcher.fetch_seq = fetch_seq

        orig_fetch_from_fasta = fasta_seqfetcher.ExonsFromGenomeFastaSeqFetcher._fetch_seq_from_fasta

        def _fetch_seq_from_fasta(seqfetcher, ac, alt_ac, alt_aln_method):
            result = orig_fetch_from_fasta(seqfetcher, ac, alt_ac, alt_aln_method)
            fasta_file = seqfetcher.contig_fastas[alt_ac]
            for exon in seqfetcher.hdp.get_tx_exons(ac, alt_ac, alt_aln_method):
                _log(fasta_file.filename, alt_ac, exon["alt_start_i"], exon["alt_end_i"])
            return result

        fasta_seqfetcher.ExonsFromGenomeFastaSeqFetcher._fetch_seq_from_fasta = _fetch_seq_from_fasta

        orig_getitem = fasta_wrapper.FastaFileContigWrapper.__getitem__

        def getitem(contig_wrapper, _slice):
            result = orig_getitem(contig_wrapper, _slice)
            _log(contig_wrapper.fasta_file.filename, contig_wrapper.contig, _slice.start, _slice.stop)
            return result

        fasta_wrapper.FastaFileContigWrapper.__getitem__ = getitem

    @classmethod
    def _dump(cls, regions):
        with open(cls.REGIONS_FILE, "w") as f:
            for region in regions:
                f.write(json.dumps(region) + "\n")
