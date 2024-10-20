from django.apps import AppConfig


class PedigreeConfig(AppConfig):
    name = 'pedigree'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel,unused-import
        from pedigree.signals import pedigree_search
        # pylint: enable=import-outside-toplevel,unused-import
