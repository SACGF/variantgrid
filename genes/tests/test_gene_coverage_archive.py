"""
Tests for GeneCoverageCollection archive helpers.

@see claude/issue_1536_data_archive_plan.md §4
"""

import os
import tempfile
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from genes.models import GeneCoverageCollection
from snpdb.archive import ArchivePreconditionError, DataArchivedError
from snpdb.models import DataState, GenomeBuild


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
        from genes.archive import archive_gene_coverage_collection
        gcc = self._make_gcc()
        with self.assertRaises(ArchivePreconditionError):
            archive_gene_coverage_collection(gcc, self.user, reason="test")

    def test_archive_idempotent(self):
        from genes.archive import archive_gene_coverage_collection
        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
            f.write(b"placeholder\n")
            tmp_path = f.name
        try:
            gcc = self._make_gcc(path=tmp_path)
            gcc.data_archived_date = timezone.now()
            gcc.save()
            with patch.object(GeneCoverageCollection, "delete_related_objects") as m:
                archive_gene_coverage_collection(gcc, self.user, reason="noop")
            m.assert_not_called()
        finally:
            os.unlink(tmp_path)

    def test_get_uncovered_gene_symbols_raises_when_archived(self):
        """ Read entry point raises DataArchivedError when the GCC is archived. """
        from genes.models import GeneSymbol
        gcc = self._make_gcc()
        gcc.data_archived_date = timezone.now()
        gcc.save()
        gene_symbols = GeneSymbol.objects.none()
        with self.assertRaises(DataArchivedError):
            gcc.get_uncovered_gene_symbols(gene_symbols, min_coverage=20)
