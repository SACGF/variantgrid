from django.apps import AppConfig

# TODO rename to OIDCAuth

# noinspection PyUnresolvedReferences
class OIDCAuthConfig(AppConfig):
    name = 'oidc_auth'

    # noinspection PyUnresolvedReferences
    def ready(self):
        # pylint: disable=import-outside-toplevel
        # imported to activate receivers
        from oidc_auth.signals import keycloak_uptime_check  # pylint: disable=unused-import
        # pylint: enable=import-outside-toplevel