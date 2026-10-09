from django.apps import AppConfig


class PedigreeConfig(AppConfig):
    name = 'pedigree'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        from pedigree.signals.pedigree_search import search_pedigree
        from snpdb.search import search_registry
        # pylint: enable=import-outside-toplevel

        search_registry.register(search_pedigree)
