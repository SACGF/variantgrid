import atexit
import json
import os

import cdot.hgvs.dataproviders.fasta_seqfetcher as fasta_seqfetcher
from django.db import connections
from django.db.migrations.loader import MigrationLoader
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

    def setup_databases(self, **kwargs):
        if self.keepdb and self.parallel > 1:
            self._drop_test_db_clones()
        old_config = super().setup_databases(**kwargs)
        if self.keepdb:
            self._check_kept_test_db_matches_disk()
        return old_config

    def _check_kept_test_db_matches_disk(self):
        """ --keepdb only migrates forwards: a migration applied while another branch was checked out keeps
            its schema and its django_migrations row after switching back, and the failures that causes
            (IntegrityError on a column the model no longer has) look nothing like the cause. """
        for connection in connections.all():
            loader = MigrationLoader(connection)
            orphans = sorted(loader.applied_migrations.keys() - loader.disk_migrations.keys())
            if orphans:
                test_db_name = connection.settings_dict["NAME"]
                orphan_list = "\n".join(f"  {app}.{name}" for app, name in orphans)
                raise SystemExit(f"Test database '{test_db_name}' has migrations applied that are not on disk "
                                 f"(applied on another branch?):\n{orphan_list}\n"
                                 f"Recreate it by running once without --keepdb.")

    def _drop_test_db_clones(self):
        """ --keepdb migrates the kept main test database but reuses an existing per-worker clone untouched,
            so the clones fall behind as migrations land. Cloning is a few seconds, so start them fresh. """
        for connection in connections.all():
            creation = connection.creation
            test_db_name = creation._get_test_db_name()
            with creation._nodb_cursor() as cursor:
                for index in range(self.parallel):
                    clone_name = connection.ops.quote_name(f"{test_db_name}_{index + 1}")
                    cursor.execute(f"DROP DATABASE IF EXISTS {clone_name}")


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
