from django.apps import AppConfig


# noinspection PyUnresolvedReferences
class ClassificationConfig(AppConfig):
    name = 'classification'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        import classification.signals  # noqa: F401  # pylint: disable=unused-import
        import classification.user_awards  # noqa: F401  # pylint: disable=unused-import  # registers award definitions
        from classification.signals.classification_search import classification_search
        from classification.signals.discordance_report_search import discordance_report_search
        from snpdb.search import search_registry
        # pylint: enable=import-outside-toplevel

        search_registry.register(classification_search, discordance_report_search)
