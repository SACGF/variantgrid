from django.apps import AppConfig

from manual import checkpoint  # noqa: F401  # pylint: disable=unused-import  # registers its system check


class ManualConfig(AppConfig):
    name = 'manual'
