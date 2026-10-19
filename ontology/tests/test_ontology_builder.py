import os
import tempfile

from django.test import TestCase
from django.utils import timezone

from library.utils.hash_utils import file_sha256sum
from ontology.models import OntologyImport
from ontology.ontology_builder import OntologyBuilder, OntologyBuilderDataUpToDateException

ABC_MD5 = "900150983cd24fb0d6963f7d28e17f72"


class TestEnsureFileHashChanged(TestCase):
    def setUp(self):
        with tempfile.NamedTemporaryFile(delete=False) as f:
            f.write(b"abc")
        self.filename = f.name
        self.addCleanup(os.unlink, self.filename)

    def _builder(self) -> OntologyBuilder:
        return OntologyBuilder(filename=self.filename, context="test_file", import_source="test")

    def _previous_import(self, data_hash: str) -> OntologyImport:
        return OntologyImport.objects.create(import_source="test", context="test_file", filename="previous",
                                             hash=data_hash, processed_date=timezone.now(), completed=True)

    def test_legacy_md5_counts_as_unchanged_and_moves_to_sha256(self):
        previous = self._previous_import(ABC_MD5)
        with self.assertRaises(OntologyBuilderDataUpToDateException):
            self._builder().ensure_file_hash_changed(self.filename)
        previous.refresh_from_db()
        self.assertEqual(previous.hash, file_sha256sum(self.filename))

    def test_different_md5_is_changed(self):
        previous = self._previous_import("0" * 32)
        builder = self._builder()
        builder.ensure_file_hash_changed(self.filename)
        self.assertEqual(builder.data_hash, file_sha256sum(self.filename))
        previous.refresh_from_db()
        self.assertEqual(previous.hash, "0" * 32)
