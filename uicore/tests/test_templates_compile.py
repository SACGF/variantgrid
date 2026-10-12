from pathlib import Path

from django.apps import apps
from django.conf import settings
from django.template.loader import get_template
from django.test import SimpleTestCase

# Tag libraries from optional apps (#1410): a template loading one only compiles when its app is installed
OPTIONAL_TAG_LIBRARY_APPS = {
    "leaflet_tags": "leaflet",
}


def _project_template_dirs() -> list[Path]:
    base_dir = Path(settings.BASE_DIR).resolve()
    dirs = [Path(d) for engine in settings.TEMPLATES for d in engine.get("DIRS", [])]
    dirs += [Path(app_config.path) / "templates" for app_config in apps.get_app_configs()]
    return [d for d in dirs
            if d.is_dir() and base_dir in d.resolve().parents and "site-packages" not in d.resolve().parts]


class TemplatesCompileTest(SimpleTestCase):
    """ Compiling catches an unknown tag/filter, a tag from a library the template doesn't {% load %},
        or an unclosed block - errors that otherwise only show when the page renders """

    def test_project_templates_compile(self):
        missing_libraries = [library for library, app in OPTIONAL_TAG_LIBRARY_APPS.items()
                             if not apps.is_installed(app)]
        for template_dir in _project_template_dirs():
            for path in template_dir.rglob("*.html"):
                if any(f"{{% load {library} %}}" in path.read_text() for library in missing_libraries):
                    continue
                with self.subTest(template=str(path)):
                    get_template(str(path.relative_to(template_dir)))
