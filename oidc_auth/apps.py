from django.apps import AppConfig

# TODO rename to OIDCAuth

# noinspection PyUnresolvedReferences
class OIDCAuthConfig(AppConfig):
    name = 'oidc_auth'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        # Registers receivers on import - noqa: F401 keeps the unused-import autofix from
        # silently unregistering them
        from oidc_auth.signals import (
            keycloak_uptime_check,  # noqa: F401  # pylint: disable=unused-import
        )
        # pylint: enable=import-outside-toplevel
