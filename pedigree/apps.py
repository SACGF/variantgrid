from django.apps import AppConfig


class PedigreeConfig(AppConfig):
    name = 'pedigree'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from pedigree.signals import pedigree_search  # noqa: F401
        # pylint: enable=import-outside-toplevel,unused-import
