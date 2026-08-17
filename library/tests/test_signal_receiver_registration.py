"""
Regression guard for the signal receivers that only register as a side effect of being imported.

Modules under ``<app>/signals/`` register their handlers by decorating module-level functions
with ``@receiver`` / ``@search_receiver`` / ``@worker_ready``. Nothing references the imported
name afterwards, so an unused-import autofix will happily delete the import from the app's
``ready()`` - which silently unregisters the handler. Nothing fails, nothing logs; searches just
stop returning results and health checks stop reporting.

This scans for modules that register on import and asserts each one actually got imported, so a
dropped ``ready()`` import is caught here instead of in production.
"""
import ast
import os
import sys
import warnings

from django.apps import apps
from django.test import SimpleTestCase

# Decorators that register a handler as a side effect of the module being imported
REGISTRATION_DECORATORS = {"receiver", "search_receiver", "worker_ready"}

# Not our runtime code, and migrations import lazily by design
SKIP_DIR_NAMES = {"migrations", "tests", "__pycache__", "static", "templates"}


def _decorator_name(decorator) -> str:
    """ The bare name of a decorator, whether it's called (@x()), dotted (@a.x) or plain (@x) """
    if isinstance(decorator, ast.Call):
        decorator = decorator.func
    if isinstance(decorator, ast.Attribute):
        return decorator.attr
    return getattr(decorator, "id", "")


def _registers_on_import(source: str) -> bool:
    try:
        with warnings.catch_warnings():
            # a few modules have regexes in non-raw strings - not this test's problem
            warnings.simplefilter("ignore", SyntaxWarning)
            tree = ast.parse(source)
    except SyntaxError:
        return False

    for node in tree.body:
        if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef):
            if any(_decorator_name(d) in REGISTRATION_DECORATORS for d in node.decorator_list):
                return True
        # module-level `some_signal.connect(handler)`
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            func = node.value.func
            if isinstance(func, ast.Attribute) and func.attr == "connect":
                return True
    return False


class SignalReceiverRegistrationTest(SimpleTestCase):
    @staticmethod
    def _registering_modules() -> list[str]:
        modules = []
        for app_config in apps.get_app_configs():
            for dir_path, dir_names, filenames in os.walk(app_config.path):
                dir_names[:] = [d for d in dir_names if d not in SKIP_DIR_NAMES]
                for filename in sorted(filenames):
                    if not filename.endswith(".py") or filename == "__init__.py":
                        continue
                    path = os.path.join(dir_path, filename)
                    with open(path) as f:
                        if not _registers_on_import(f.read()):
                            continue
                    dotted = os.path.relpath(path, app_config.path)[:-len(".py")].replace(os.sep, ".")
                    modules.append(f"{app_config.name}.{dotted}")
        return modules

    def test_scan_finds_the_signals_modules(self):
        """ Guards the scan itself - a rename that breaks it would otherwise pass vacuously """
        modules = self._registering_modules()
        self.assertGreater(len(modules), 20, f"Expected to find the signals modules, got {modules}")
        self.assertIn("snpdb.signals.trio_search", modules)

    def test_registering_modules_are_imported(self):
        """ Every module that registers on import must be imported by its app's ready() """
        not_imported = [m for m in self._registering_modules() if m not in sys.modules]
        self.assertEqual([], not_imported,
                         "These modules register signal receivers on import but were never imported - "
                         "add them to the app's AppConfig.ready() with a '# noqa: F401' so the "
                         "unused-import autofix leaves them alone")
