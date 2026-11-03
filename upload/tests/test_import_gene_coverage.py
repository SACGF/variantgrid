"""ImportGeneCoverageTask: the collection's data state on failure, and the record count it reports (#2617)."""
import os
import tempfile

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from genes.enums import AnnotationConsortium
from genes.models import CanonicalTranscript, CanonicalTranscriptCollection, GeneCoverageCollection
from snpdb.enums import DataState, ImportSource
from snpdb.models import GenomeBuild
from upload.models import FileUpload
from upload.tasks.import_gene_coverage_task import ImportGeneCoverageTask

GENE_COVERAGE_FILE = os.path.join(settings.BASE_DIR, "seqauto", "test_data", "clinical_hg38", "idt_haem",
                                  "Haem_20_999_201231_M02027_0112_000000000_JFT79", "4_QC", "bam_stats",
                                  "samples", "fake_sample_1.per_gene_coverage.tsv.gz")


class TestImportGeneCoverage(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='testuser')[0]
        cls.canonical_collection = CanonicalTranscriptCollection.objects.create(
            genome_build=GenomeBuild.legacy_build(), annotation_consortium=AnnotationConsortium.REFSEQ)
        for symbol, transcript in [("PLCXD1", "NM_018390.4"), ("KIF1B", "NM_183416.4")]:
            CanonicalTranscript.objects.create(collection=cls.canonical_collection,
                                               original_gene_symbol=symbol, original_transcript=transcript)

    def setUp(self):
        self.file_upload = FileUpload.objects.create(path=GENE_COVERAGE_FILE, import_source=ImportSource.COMMAND_LINE,
                                                     name="fake_sample_1", user=self.user)
        processing_dir = tempfile.TemporaryDirectory()
        self.addCleanup(processing_dir.cleanup)
        settings_override = override_settings(IMPORT_PROCESSING_DIR=processing_dir.name,
                                              GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID=self.canonical_collection.pk,
                                              SEQAUTO_QC_GENE_COVERAGE_STORE_ALL=False,
                                              SEQAUTO_QC_GENE_COVERAGE_STORE_CANONICAL=True)
        settings_override.enable()
        self.addCleanup(settings_override.disable)

    def _collections(self):
        return GeneCoverageCollection.objects.filter(path=GENE_COVERAGE_FILE)

    def test_counts_canonical_rows_when_not_storing_all(self):
        num_records = ImportGeneCoverageTask.process_items(self.file_upload)
        self.assertEqual(2, num_records)
        self.assertEqual(DataState.COMPLETE, self._collections().get().data_state)

    def test_no_default_canonical_collection_creates_no_collection(self):
        with override_settings(GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID=None):
            with self.assertRaises(ValueError):
                ImportGeneCoverageTask.process_items(self.file_upload)
        self.assertFalse(self._collections().exists())

    def test_failed_load_sets_collection_to_error(self):
        # Most genes have no canonical transcript match, which load_from_file rejects
        CanonicalTranscript.objects.filter(collection=self.canonical_collection).update(original_transcript="NM_000000.1")
        with self.assertRaises(ValueError):
            ImportGeneCoverageTask.process_items(self.file_upload)
        self.assertEqual(DataState.ERROR, self._collections().get().data_state)
