"""
The explicit search registry (#1690): every @search_receiver declared in an installed app is registered from that
app's AppConfig.ready(), and SearchReceiver applies the dispatch gates and error capture itself.
"""
import ast
import os
import re
import warnings

from django.apps import apps
from django.contrib.auth.models import User
from django.test import SimpleTestCase

from library.enums.log_level import LogLevel
from library.preview_request import PreviewData, PreviewProxyModel
from library.tests.test_signal_receiver_registration import SKIP_DIR_NAMES, _decorator_name
from snpdb.search import (
    SearchInput,
    SearchInputInstance,
    SearchRegistry,
    search_receiver,
    search_registry,
)


def _declared_search_receivers(source: str) -> list[str]:
    """ Names of the module-level functions decorated with @search_receiver """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SyntaxWarning)
            tree = ast.parse(source)
    except SyntaxError:
        return []
    return [node.name for node in tree.body
            if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef)
            and any(_decorator_name(d) == "search_receiver" for d in node.decorator_list)]


class SearchRegistryCompletenessTest(SimpleTestCase):
    @staticmethod
    def _declared() -> set[tuple[str, str]]:
        declared = set()
        for app_config in apps.get_app_configs():
            for dir_path, dir_names, filenames in os.walk(app_config.path):
                dir_names[:] = [d for d in dir_names if d not in SKIP_DIR_NAMES]
                for filename in filenames:
                    if not filename.endswith(".py"):
                        continue
                    path = os.path.join(dir_path, filename)
                    with open(path) as f:
                        names = _declared_search_receivers(f.read())
                    if not names:
                        continue
                    dotted = os.path.relpath(path, app_config.path)[:-len(".py")].replace(os.sep, ".")
                    module = f"{app_config.name}.{dotted}"
                    declared.update((module, name) for name in names)
        return declared

    def test_scan_finds_the_search_receivers(self):
        """ Guards the scan itself - a rename that breaks it would otherwise pass vacuously """
        declared = self._declared()
        self.assertGreater(len(declared), 40, f"Expected to find the search receivers, got {declared}")
        self.assertIn(("snpdb.signals.trio_search", "search_trio"), declared)

    def test_every_declared_receiver_is_registered(self):
        registered = {(r.func.__module__, r.func.__name__) for r in search_registry.receivers}
        missing = sorted(f"{module}.{name}" for module, name in self._declared() - registered)
        self.assertEqual([], missing,
                         "These @search_receiver functions are never registered - add each name to its app's "
                         "AppConfig.ready() search_registry.register(...) call")


def _stub_receiver(category: str = "Stub", **kwargs):
    stub_type = PreviewProxyModel(category, "fa-solid fa-flask")

    def stub_search(search_input: SearchInputInstance):
        yield PreviewData(category=category, identifier=search_input.search_string)

    return search_receiver(search_type=stub_type, **kwargs)(stub_search)


def _search_input(search_string: str = "stub", is_superuser: bool = False, classify: bool = False) -> SearchInput:
    return SearchInput(user=User(username="stub_user", is_superuser=is_superuser), search_string=search_string,
                       genome_build_preferred=None, classify=classify)


class SearchReceiverTest(SimpleTestCase):

    def test_disabled_is_not_visible(self):
        self.assertTrue(_stub_receiver().visible_to(_search_input()))
        self.assertFalse(_stub_receiver(enabled=False).visible_to(_search_input()))

    def test_admin_only_is_visible_only_to_superuser(self):
        receiver = _stub_receiver(admin_only=True)
        self.assertFalse(receiver.visible_to(_search_input()))
        self.assertTrue(receiver.visible_to(_search_input(is_superuser=True)))

    def test_classify_only_admits_variant_searches(self):
        self.assertFalse(_stub_receiver().visible_to(_search_input(classify=True)))
        self.assertTrue(_stub_receiver(category="Variant").visible_to(_search_input(classify=True)))

    def test_pattern_not_matched(self):
        response = _stub_receiver(pattern=re.compile(r"^\d+$")).search(_search_input("stub"))
        self.assertFalse(response.matched_pattern)
        self.assertEqual([], response.results)
        self.assertEqual(0.0, response.duration_seconds)

    def test_exception_becomes_error_message(self):
        def broken_search(search_input: SearchInputInstance):
            raise RuntimeError("stub failure")
            yield  # pylint: disable=unreachable

        receiver = search_receiver(search_type=PreviewProxyModel("Stub", "fa-solid fa-flask"))(broken_search)
        with self.assertLogs(level="ERROR"):
            response = receiver.search(_search_input())
        self.assertTrue(response.matched_pattern)
        self.assertEqual(["stub failure"], [m.message for m in response.messages_overall])
        self.assertEqual(LogLevel.ERROR, response.messages_overall[0].severity)

    def test_register_twice_raises(self):
        registry = SearchRegistry()
        receiver = _stub_receiver()
        registry.register(receiver)
        with self.assertRaises(ValueError):
            registry.register(receiver)
