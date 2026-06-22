"""
Regression guard for SACGF/variantgrid#1581.

Chaining ``@classmethod`` on top of ``@property`` (in either order) was deprecated in
Python 3.11 and **removed in Python 3.13**. On 3.13+ accessing the attribute returns a
``method`` object instead of invoking the function, which fails at runtime (e.g.
``TypeError: 'method' object is not iterable``) and is only caught when the relevant
code path executes. This test statically scans the codebase so any reintroduction is
caught immediately rather than at runtime.
"""
import ast
import os

from django.test import SimpleTestCase, TestCase

from classification.models.classification_variant_info_models import ImportedAlleleInfo
from genes.hgvs.hgvs_matcher import ClinGenHGVSConverter
from snpdb.models import GenomeBuild
from variantgrid.settings.components.settings_paths import BASE_DIR

# Directories under BASE_DIR that are not our source (third-party / generated / tooling).
_SKIP_DIR_NAMES = {
    ".git", ".venv", "venv", "node_modules", "__pycache__",
    "static", "media_root", "private_data", ".idea", ".mypy_cache",
}


def _decorator_names(node: ast.AST) -> list[str]:
    names = []
    for dec in node.decorator_list:
        target = dec.func if isinstance(dec, ast.Call) else dec
        if isinstance(target, ast.Name):
            names.append(target.id)
        elif isinstance(target, ast.Attribute):
            names.append(target.attr)
    return names


def _iter_project_python_files():
    for root, dirs, files in os.walk(BASE_DIR):
        dirs[:] = [d for d in dirs if d not in _SKIP_DIR_NAMES]
        for fn in files:
            if fn.endswith(".py"):
                yield os.path.join(root, fn)


class ChainedClassmethodPropertyAuditTest(SimpleTestCase):
    """ Codebase-wide audit (no DB required). """

    def test_no_chained_classmethod_property(self):
        offenders = []
        for path in _iter_project_python_files():
            with open(path, encoding="utf-8") as f:
                source = f.read()
            try:
                tree = ast.parse(source, filename=path)
            except SyntaxError:
                continue  # not valid on this interpreter; not our concern here
            for node in ast.walk(tree):
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    names = set(_decorator_names(node))
                    if "classmethod" in names and {"property", "cached_property"} & names:
                        rel = os.path.relpath(path, BASE_DIR)
                        offenders.append(f"{rel}:{node.lineno} {node.name} -> {sorted(names)}")

        self.assertEqual(
            offenders, [],
            "Found @classmethod chained with @property/@cached_property (removed in "
            "Python 3.13, see issue #1581). Convert to a plain @classmethod and add () "
            "to call sites:\n" + "\n".join(offenders)
        )


class FixedClassMethodsCallableTest(TestCase):
    """ Exercises the methods that were converted from chained decorators in #1581.

        On Python 3.13+ a chained @classmethod/@property would not invoke the function,
        so calling these with () would return a descriptor/method object rather than the
        expected collection of genome builds. """

    def test_clingen_supported_builds_is_callable(self):
        builds = ClinGenHGVSConverter.supported_builds()
        self.assertEqual(set(builds), {GenomeBuild.grch37(), GenomeBuild.grch38()})
        self.assertTrue(ClinGenHGVSConverter.build_supported(GenomeBuild.grch37()))

    def test_imported_allele_info_supported_genome_builds_is_callable(self):
        builds = ImportedAlleleInfo.supported_genome_builds()
        self.assertIn(GenomeBuild.grch37(), builds)
        self.assertIn(GenomeBuild.grch38(), builds)
