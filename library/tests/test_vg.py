"""
Tests for the DB-free parts of `manage.py vg` (library/vg): import resolution, test selection rules,
the AST parsers behind the tasks / signals / settings maps, and Markdown rendering.
"""
import ast
import tempfile
from pathlib import Path

from django.test import SimpleTestCase

from library.vg.import_graph import (
    ImportGraph,
    iter_import_statements,
    relative_base,
    resolve_from_import,
)
from library.vg.maps import settings as settings_map
from library.vg.maps import signals as signals_map
from library.vg.maps import tasks as tasks_map
from library.vg.markdown import MapTable, render_markdown
from library.vg.test_selection import select_tests


def _graph(edges: dict[str, list[str]]) -> ImportGraph:
    graph = ImportGraph()
    for importer, imported in edges.items():
        for target in imported:
            graph.add_edge(importer, target)
    return graph


class ImportGraphTest(SimpleTestCase):

    def test_relative_import_base(self):
        self.assertEqual(relative_base("snpdb.models.models_variant", False, 1, "models_genome"), "snpdb.models.models_genome")
        self.assertEqual(relative_base("snpdb.models", True, 1, "models_genome"), "snpdb.models.models_genome")
        self.assertEqual(relative_base("snpdb.views.views", False, 2, "models"), "snpdb.models")

    def test_from_import_prefers_submodule_over_symbol(self):
        known = {"snpdb.models", "snpdb.models.models_variant"}
        self.assertEqual(resolve_from_import(known, "snpdb.models", "models_variant"), "snpdb.models.models_variant")
        self.assertEqual(resolve_from_import(known, "snpdb.models", "Variant"), "snpdb.models")
        self.assertIsNone(resolve_from_import(known, "django.db", "models"))

    def test_import_statements_found_inside_functions_and_try_blocks(self):
        tree = ast.parse("import a\ndef f():\n    from b import c\ntry:\n    import d\nexcept ImportError:\n    import e\n")
        names = sorted(getattr(n, "module", None) or n.names[0].name for n in iter_import_statements(tree))
        self.assertEqual(names, ["a", "b", "d", "e"])

    def test_dependents_are_transitive(self):
        graph = _graph({"app.tests.test_x": ["app.views"], "app.views": ["app.models"], "other": ["app.views"]})
        self.assertEqual(graph.dependents("app.models"), {"app.views", "app.tests.test_x", "other"})


class TestSelectionTest(SimpleTestCase):
    GRAPH = _graph({
        "snpdb.tests.test_variant": ["snpdb.models.models_variant"],
        "analysis.tests.test_merge_node": ["analysis.models.nodes.filters.merge_node"],
        "analysis.models.nodes.filters.merge_node": ["snpdb.models.models_variant"],
        "snpdb.views.views_data": ["snpdb.models.models_variant"],
    })

    def _labels(self, changed):
        return select_tests(changed=changed, graph=self.GRAPH).labels

    def test_changed_module_selects_transitive_test_importers(self):
        self.assertEqual(self._labels(["snpdb/models/models_variant.py"]),
                         ["analysis.tests.test_merge_node", "snpdb.tests.test_variant"])

    def test_changed_test_module_selects_itself(self):
        self.assertEqual(self._labels(["snpdb/tests/test_variant.py"]), ["snpdb.tests.test_variant"])

    def test_template_and_static_select_the_owning_apps_url_test(self):
        self.assertEqual(self._labels(["analysis/templates/analysis/analysis.html"]), ["analysis.tests.test_urls"])
        self.assertEqual(self._labels(["variantgrid/static_files/default_static/js/analysis/x.js"]), ["analysis.tests.test_urls"])

    def test_migration_selects_whole_app_and_absorbs_module_labels(self):
        selection = select_tests(changed=["snpdb/migrations/0001_initial.py", "snpdb/tests/test_variant.py"], graph=self.GRAPH)
        self.assertEqual(selection.labels, ["snpdb.tests"])
        self.assertEqual(selection.reasons["snpdb.tests"], ["snpdb/migrations/0001_initial.py", "snpdb/tests/test_variant.py"])

    def test_docs_and_data_select_nothing(self):
        self.assertEqual(self._labels(["claude/plans/x.md", "snpdb/CLAUDE.md", "requirements.txt"]), [])


class TasksMapTest(SimpleTestCase):

    def test_routes_resolve_queue_constants(self):
        routes = tasks_map.task_routes()
        self.assertEqual(routes["analysis.tasks.node_update_tasks.update_node_task"], "analysis_workers")

    def test_enqueue_targets_cover_delay_signature_and_task_classes(self):
        tree = ast.parse("a.delay(1)\nb.si(2)\nMyTask().apply_async()\nunrelated.filter()\nmod.c.s()\n")
        self.assertEqual(tasks_map.enqueue_targets(tree), {"a", "b", "MyTask", "c"})

    def test_decorator_queue(self):
        decorator = ast.parse("@celery.shared_task(queue='web_workers')\ndef f(): pass").body[0].decorator_list[0]
        self.assertEqual(tasks_map.decorator_queue(decorator), "web_workers")
        self.assertEqual(tasks_map.decorator_name(decorator), "shared_task")


class SignalsMapTest(SimpleTestCase):

    def test_receiver_decorator_signals_and_sender(self):
        decorator = ast.parse("@receiver([post_save, pre_delete], sender=Variant)\ndef f(): pass").body[0].decorator_list[0]
        self.assertEqual(signals_map.receiver_signal_names(decorator), ["post_save", "pre_delete"])
        self.assertEqual(signals_map.receiver_sender(decorator), "Variant")
        keyword_form = ast.parse("@receiver(signal=signals.my_signal)\ndef f(): pass").body[0].decorator_list[0]
        self.assertEqual(signals_map.receiver_signal_names(keyword_form), ["my_signal"])

    def test_signal_call_detection(self):
        self.assertTrue(signals_map.is_signal_call(ast.parse("x = django.dispatch.Signal()").body[0].value))
        self.assertFalse(signals_map.is_signal_call(ast.parse("x = Signal").body[0].value))


class SettingsMapTest(SimpleTestCase):

    def test_assigned_names_distinguish_reassignment_from_mutation(self):
        source = ("FOO = 1\nBAR['x'] = 2\nBAZ.update({'a': 1})\nQUX += [3]\nlower = 4\n_PRIVATE = 5\n"
                  "ANNOTATION[BUILD]['enabled'] = True\n")
        with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as f:
            f.write(source)
        names = settings_map.assigned_names(Path(f.name))
        self.assertEqual(names, [("FOO", False), ("BAR", True), ("BAZ", True), ("QUX", False), ("ANNOTATION", True)])


class MarkdownTest(SimpleTestCase):

    def test_pipes_and_newlines_are_escaped_in_cells(self):
        table = MapTable(title="T", columns=["A"], rows=[["a|b\nc"]])
        rendered = render_markdown("h", "intro", [table])
        self.assertIn("| a\\|b c |", rendered)
        self.assertIn("## T\n", rendered)

    def test_empty_table_renders_placeholder(self):
        rendered = render_markdown("h", "intro", [MapTable(title="Empty", columns=["A"])])
        self.assertIn("(none)", rendered)
