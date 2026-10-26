import atexit
import json
import os
import shutil

import cdot.hgvs.dataproviders.fasta_seqfetcher as fasta_seqfetcher
from django.conf import settings
from django.db import DEFAULT_DB_ALIAS, connections
from django.db.migrations.loader import MigrationLoader
from django.test.runner import DiscoverRunner

import library.genomics.fasta_wrapper as fasta_wrapper
from annotation.fake_data import get_fake_annotation_version
from snpdb.models.models_genome import GenomeBuild


class VariantGridTestRunner(DiscoverRunner):
    """ Seeds the kept test database before cloning it for parallel workers, and cleans up the per-suite
        import processing dir. External service clients are swapped for recorded-data mocks in settings
        (CLINGEN_ALLELE_REGISTRY_API_CLASS etc.), not here - a spawned worker never runs this runner's
        setup_test_environment """

    def teardown_test_environment(self, **kwargs):
        super().teardown_test_environment(**kwargs)
        # settings.IMPORT_PROCESSING_DIR is a per-suite temp dir (see default_settings UNIT_TEST block),
        # so whatever the run left in it - a pipeline dir, fake annotation scratch - goes with it
        shutil.rmtree(settings.IMPORT_PROCESSING_DIR, ignore_errors=True)

    def setup_databases(self, **kwargs):
        """ Django creates the main test database and clones it for the workers in the one call; the fake
            annotation versions are seeded between the two, so every worker starts with them (below). """
        if self.keepdb and self.parallel > 1:
            self._drop_test_db_clones()
        keepdb = self.keepdb
        if keepdb and self._kept_test_dbs_not_matching_disk(kwargs.get("aliases", {DEFAULT_DB_ALIAS})):
            self.keepdb = False  # Django's own drop and recreate, as a run without --keepdb does
        parallel = self.parallel
        self.parallel = 1  # create the main test database only - we clone it ourselves, after seeding
        try:
            old_config = super().setup_databases(**kwargs)
        finally:
            self.parallel = parallel
            self.keepdb = keepdb  # teardown still keeps the recreated database
        test_connections = [connection for connection, _old_name, _created in old_config]
        # A run of SimpleTestCase-only modules sets up no database at all - seeding then would write the
        # fixture rows into the real one
        if any(connection.alias == DEFAULT_DB_ALIAS for connection in test_connections):
            self._seed_fake_annotation_versions()
        if parallel > 1:
            for connection, _old_name, created in old_config:
                if created:
                    for index in range(parallel):
                        connection.creation.clone_test_db(suffix=str(index + 1), verbosity=self.verbosity,
                                                          keepdb=self.keepdb)
        return old_config

    def _seed_fake_annotation_versions(self):
        """ 240-odd test classes call get_fake_annotation_version in setUpTestData, each re-importing the
            ontology OWL and rebuilding an identical version chain inside a transaction that is then rolled
            back. Seeded here (outside any transaction, before the workers are cloned) the per-class calls
            become a handful of get_or_create lookups.

            A test that counts AnnotationVersion / VariantAnnotationVersion / OntologyTerm rows, or creates
            its own ACTIVE VariantAnnotationVersion for a build, now sees these - use the fixture. """
        for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
            get_fake_annotation_version(genome_build)

    @staticmethod
    def _kept_test_dbs_not_matching_disk(aliases) -> bool:
        """ --keepdb only migrates forwards, so a kept test database holding migrations that aren't on disk - the
            pre-squash history (#736), or one applied while another branch was checked out - fails in ways that
            look nothing like the cause (InconsistentMigrationHistory, an IntegrityError on a column the model no
            longer has). Such a database is rebuilt rather than migrated: it is only a cache of the migrations. """
        disk_migrations = MigrationLoader(None, ignore_no_migrations=True).disk_migrations.keys()
        stale = False
        for alias in aliases:
            connection = connections[alias]
            test_db_name = connection.creation._get_test_db_name()
            with connection.creation._nodb_cursor() as cursor:
                cursor.execute("SELECT 1 FROM pg_database WHERE datname = %s", [test_db_name])
                if not cursor.fetchone():
                    continue
            # Raw SQL: MigrationRecorder queries through the alias, ie the real database until setup swaps it
            test_connection = connection.__class__({**connection.settings_dict, "NAME": test_db_name}, alias)
            try:
                with test_connection.cursor() as cursor:
                    if "django_migrations" not in test_connection.introspection.table_names(cursor):
                        continue
                    cursor.execute("SELECT app, name FROM django_migrations")
                    orphans = sorted(set(cursor.fetchall()) - disk_migrations)
            finally:
                test_connection.close()
            if orphans:
                example = ", ".join(f"{app}.{name}" for app, name in orphans[:3])
                print(f"Test database '{test_db_name}' has {len(orphans)} migrations applied that aren't on disk "
                      f"({example}, ...) - recreating it")
                stale = True
        return stale

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
        with open(cls.REGIONS_FILE, "w", encoding="utf-8") as f:
            for region in regions:
                f.write(json.dumps(region) + "\n")
