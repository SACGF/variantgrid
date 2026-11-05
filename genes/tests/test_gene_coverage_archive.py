"""
Tests for GeneCoverageCollection archive helpers.

@see https://github.com/SACGF/variantgrid/issues/1536
"""

import tempfile

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from eventlog.models import Event
from genes.archive import archive_gene_coverage_collection
from genes.models import GeneCoverageCollection, GeneSymbol
from library.django_utils.data_archive_mixin import DataArchivedError
from snpdb.archive import ArchivePreconditionError
from snpdb.enums import DataState
from snpdb.models import GenomeBuild


class GCCArchiveTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="gcc_user")
        cls.gb = GenomeBuild.grch38()

    def _make_gcc(self, path="") -> GeneCoverageCollection:
        return GeneCoverageCollection.objects.create(
            path=path or "/tmp/nonexistent.tsv",
            data_state=DataState.COMPLETE,
            genome_build=self.gb,
        )

    def test_archive_raises_when_source_missing(self):
        gcc = self._make_gcc()
        with self.assertRaises(ArchivePreconditionError):
            archive_gene_coverage_collection(gcc, self.user, reason="test")

    def test_archive_idempotent(self):
        """ Archiving an archived collection keeps the first archive's record and logs nothing """
        with tempfile.NamedTemporaryFile(suffix=".tsv") as f:
            gcc = self._make_gcc(path=f.name)
            archived_date = timezone.now()
            gcc.data_archived_date = archived_date
            gcc.data_archived_by = self.user
            gcc.data_archive_reason = "first"
            gcc.save()
            other_user = User.objects.create_user(username="gcc_other_user")
            archive_gene_coverage_collection(gcc, other_user, reason="second")

        gcc.refresh_from_db()
        self.assertEqual((archived_date, self.user, "first"),
                         (gcc.data_archived_date, gcc.data_archived_by, gcc.data_archive_reason))
        self.assertFalse(Event.objects.filter(name="gene_coverage_collection_archived").exists())

    def test_get_uncovered_gene_symbols_raises_when_archived(self):
        """ Read entry point raises DataArchivedError when the GCC is archived. """
        gcc = self._make_gcc()
        gcc.data_archived_date = timezone.now()
        gcc.save()
        gene_symbols = GeneSymbol.objects.none()
        with self.assertRaises(DataArchivedError):
            gcc.get_uncovered_gene_symbols(gene_symbols, min_coverage=20)
