import importlib
import os
import tempfile

from django.apps import apps
from django.test import TestCase, override_settings

from snpdb.models import CachedGeneratedFile

migration = importlib.import_module("snpdb.migrations.0009_clear_unsalted_generated_files")


class ClearUnsaltedGeneratedFilesTest(TestCase):

    def _cgf(self, generated_dir, generator, params_hash, relative_filename):
        filename = os.path.join(generated_dir, generator, relative_filename)
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        with open(filename, "w") as f:
            f.write("x")
        return CachedGeneratedFile.objects.create(generator=generator, params_hash=params_hash, filename=filename)

    def test_clears_unsalted_graphs_and_pedigree_charts_only(self):
        with tempfile.TemporaryDirectory() as generated_dir, override_settings(GENERATED_DIR=generated_dir):
            unsalted = self._cgf(generated_dir, "HomozygosityPercentGraph", "aaa", "aaa.png")
            salted = self._cgf(generated_dir, "HomozygosityPercentGraph", "bbb", "c0ffee.png")
            pedigree = self._cgf(generated_dir, "PedigreeChart", "ddd", "ddd.svg")
            export = self._cgf(generated_dir, "export_sample_to_downloadable_file", "eee", "1234/eee.csv")
            orphan = os.path.join(generated_dir, "HomozygosityPercentGraph", "orphan.png")
            open(orphan, "w").close()

            migration._clear_unsalted_generated_files(apps, None)

            remaining = set(CachedGeneratedFile.objects.values_list("pk", flat=True))
            self.assertEqual(remaining, {salted.pk, export.pk})
            for cgf in (unsalted, pedigree):
                self.assertFalse(os.path.exists(cgf.filename))
            self.assertFalse(os.path.exists(orphan))
            for cgf in (salted, export):
                self.assertTrue(os.path.exists(cgf.filename))
