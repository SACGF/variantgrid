import re
from urllib.parse import urlencode

import jwt
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import Group, User
from django.core.exceptions import PermissionDenied, SuspiciousOperation, ValidationError
from django.core.validators import EmailValidator
from django.utils.html import format_html, format_html_join
from django.utils.safestring import mark_safe
from mozilla_django_oidc.auth import OIDCAuthenticationBackend
from mozilla_django_oidc.utils import import_from_settings

from library.log_utils import report_message
from snpdb.models import UserSettingsOverride

_email_validator = EmailValidator()
_CONTROL_CHAR_RE = re.compile(r'[\x00-\x1f\x7f]')


def _sanitize_claim_str(value, max_length):
    """Strip control characters and enforce max length on a claim string."""
    return _CONTROL_CHAR_RE.sub('', value)[:max_length]


class VariantGridOIDCAuthenticationBackend(OIDCAuthenticationBackend):
    """
    Keycloak is the source of truth for users: every login (browser callback, or each DRF bearer-token request via
    mozilla's OIDCAuthentication) rewrites the local user, its groups and admin flags from the claims.
    """

    def create_user(self, claims):
        # Unsaved, so a refused first login (wrong environment, no labs, maintenance) leaves no row behind
        user = self.UserModel()
        user.set_unusable_password()
        return self.create_or_update(user, claims)

    def update_user(self, user, claims):
        return self.create_or_update(user, claims)

    def filter_users_by_claims(self, claims):
        """Check username first"""
        username = claims.get('preferred_username')
        if username:
            by_username = self.UserModel.objects.filter(username=username)
            if by_username.exists():
                return by_username

        # Return all users matching the specified email
        email = claims.get('email')
        if not email:
            return self.UserModel.objects.none()
        return self.UserModel.objects.filter(email__iexact=email)

    def get_userinfo(self, access_token, id_token, payload):
        user_info = super().get_userinfo(access_token, id_token, payload)
        if id_token is None:
            # DRF bearer token - the browser flow always has an ID token, verified against our client ID
            self._verify_bearer_token_client(access_token)
        return user_info

    def _verify_bearer_token_client(self, access_token):
        """
        Keycloak's userinfo endpoint accepts an access token issued to any client in the realm, and one realm serves
        every Shariant environment - so also require the token to have been issued to (azp) or for (aud) this deployment.
        Userinfo has already checked this exact token's signature and expiry, so its claims are read without a JWKS fetch.
        """
        try:
            token_claims = jwt.decode(access_token, options={"verify_signature": False})
        except jwt.DecodeError as e:
            raise SuspiciousOperation("Bearer token is not a JWT") from e

        allowed_clients = {self.OIDC_RP_CLIENT_ID, *settings.OIDC_API_EXTRA_CLIENT_IDS}
        azp = token_claims.get("azp")
        audience = token_claims.get("aud") or []
        if isinstance(audience, str):
            audience = [audience]
        if azp not in allowed_clients and self.OIDC_RP_CLIENT_ID not in audience:
            report_message(f"Rejected API bearer token issued to client {azp!r} (aud={audience}), "
                           f"allowed clients: {sorted(allowed_clients)}", level='warning')
            raise SuspiciousOperation(f"Bearer token was issued to client {azp!r}")

    def _refuse_login(self, user: User, message, deactivate: bool, extra_tags: str = "") -> User:
        """
        Browser login: returns the user inactive, which mozilla's callback turns into a login failure showing `message`.
        Returning None instead would count as a failed login to django-axes, keyed by IP with no username.
        DRF bearer auth has no request (it calls get_or_create_user without authenticate()) and would accept an inactive
        user as authenticated, so that raises PermissionDenied - a 403.
        """
        user.is_active = False
        if deactivate and user.pk:
            user.save()

        request = getattr(self, "request", None)
        if request is None:
            raise PermissionDenied()
        messages.add_message(request, messages.ERROR, message, extra_tags=extra_tags)
        return user

    def create_or_update(self, user: User, claims):

        missing_claims = set()
        for claim in ('preferred_username', 'email', 'sub', 'groups'):
            if not claims.get(claim):
                missing_claims.add(claim)
        if missing_claims:
            report_message(f"Missing claims {missing_claims}", level='error')

        username = claims.get('preferred_username')
        if not username:
            raise SuspiciousOperation("No preferred_username in OIDC claims")

        email = claims.get('email', '')
        try:
            _email_validator(email)
        except ValidationError:
            raise SuspiciousOperation(f"Invalid email in OIDC claims: {email!r}") from None
        email = _sanitize_claim_str(email, 254)

        # Copy over basic details from open ID connect
        # Assume there will be no user-name clashes
        user.username = _sanitize_claim_str(username, 254)
        user.email = email
        user.first_name = _sanitize_claim_str(claims.get('given_name', ''), 150)
        user.last_name = _sanitize_claim_str(claims.get('family_name', ''), 150)
        sub = claims.get('sub', None)

        # Work out what groups the user has joined/left since their last login
        django_groups = set(user.groups.values_list("name", flat=True)) if user.pk else set()
        # groups with be in the form of '/variantgrid/some_group_1', '/variantgrid/some_group_2', '/unrelated'
        # convert it so we get 'some_group_1', 'some_group_2'

        user.is_active = True
        all_claim_groups = claims.get("groups") or []
        if settings.OIDC_REQUIRED_GROUP and settings.OIDC_REQUIRED_GROUP not in all_claim_groups:
            report_message(f"User {user.username} attempted to login but lacked the group permission {settings.OIDC_REQUIRED_GROUP} - belongs to groups {all_claim_groups}", level='error')

            # Please try our test environment <a href="https://test.shariant.org.au">https://test.shariant.org.au</a>

            allowed_environments_map = {
                "/variantgrid/shariant_demo": """Please try out demo environment <a href="https://demo.shariant.org.au">https://demo.shariant.org.au</a>""",
                "/variantgrid/shariant_test": """Please try out test environment <a href="https://test.shariant.org.au">https://test.shariant.org.au</a>""",
                "/variantgrid/shariant_production": """Please try out production environment <a href="https://shariant.org.au">https://shariant.org.au</a>""",
                "/variantgrid/shariant_security": """Please try out security testing environment <a href="https://test2.shariant.org.au">https://test2.shariant.org.au</a>""",
            }
            allowed_environment_list = []
            for group, message in allowed_environments_map.items():
                if group in all_claim_groups:
                    allowed_environment_list.append(message)

            # Note that the user has provided a correct username and password from our system, but tried to log into the wrong account
            # No security issue reflecting their email back to them
            message = format_html("This account <i>{}</i> is not authorised for this environment.", user.email)
            message += format_html_join("", "<br/>{}", ((mark_safe(env),) for env in allowed_environment_list))

            return self._refuse_login(user, message, deactivate=True, extra_tags="html")

        oauth_groups = [g.split('/')[1:] for g in all_claim_groups]

        associations = [g[1:] for g in oauth_groups if len(g) > 1 and g[0] == 'associations']
        # Should make 'variantgrid' setting configurable at some point in case there are 2 variantgrid installations on the same OAuth instance
        variant_grid_groups = [g[1:] for g in oauth_groups if len(g) > 1 and g[0] == 'variantgrid']

        if not associations and not variant_grid_groups:
            return self._refuse_login(user, "This account doesn't belong to any labs.", deactivate=True)

        # everyone with a login is considered part of the public group
        groups: set[str] = set()
        groups.add(settings.PUBLIC_GROUP_NAME)
        groups.add(settings.LOGGED_IN_USERS_GROUP_NAME)

        # currently only variantgrid permission
        is_super_user = False
        is_bot = False
        is_tester = False
        for vg in variant_grid_groups:
            permission = '/'.join(vg)
            if permission == 'admin':
                is_super_user = True
            elif permission == 'bot':
                groups.add('variantgrid/bot')
                is_bot = True
            elif permission == 'tester':
                groups.add('variantgrid/tester')
                is_tester = True

        if settings.MAINTENANCE_MODE:
            if is_tester:
                # testers are allowed to login during maintenance mode
                pass
            elif (not is_super_user) or is_bot:
                # don't want bots logging in during maintenance mode
                return self._refuse_login(user, "Non-administrator logins have temporary been disabled.",
                                          deactivate=False)

        user.is_superuser = is_super_user
        user.is_staff = is_super_user
        user.save()

        user_settings_override, _ = UserSettingsOverride.objects.get_or_create(user=user)
        user_settings_override.oauth_sub = sub

        # for nested groups, adds each level of nesting as its own group
        # e.g. association/fake_pathology/lab_1 will be added as
        # "fake_pathology/lab_1"
        # "fake_pathology"
        # (remember that the 'association' prefix is removed)
        for assoc in associations:
            assoc_groups = set()
            for i in range(len(assoc) + 1):
                parts = assoc[0:i]
                assoc_groups.add('/'.join(parts))
                groups.update(assoc_groups)

        removed_groups = django_groups.difference(groups)
        added_groups = groups.difference(django_groups)

        for removed_group in removed_groups:
            group = Group.objects.get(name=removed_group)
            user.groups.remove(group)

        for added_group in added_groups:
            # note that we trust the OIDC connector as it can already make admins
            # and sometimes group permissions are setup in KeyCloak before they are in the app
            # so happy for this to make users
            # Groups have already been verified to be inside variantgrid/ or associations/
            group, _ = Group.objects.get_or_create(name=added_group)
            user.groups.add(group)

        # ensures default lab is valid and sets it if there are labs
        # and the default lab is blank
        user_settings_override.auto_set_default_lab()
        user_settings_override.save()
        return user


def provider_logout(request) -> str:
    redirect = import_from_settings("LOGOUT_REDIRECT_URL", "")
    oidc_id_token = request.session.get("oidc_id_token")
    if not oidc_id_token:
        # Password login through ModelBackend (/accounts/, /admin/) - there's no Keycloak session to end
        return redirect or "/"

    oidc_logout = import_from_settings("KEY_CLOAK_PROTOCOL_BASE", "") + "/logout"
    if redirect:
        oidc_logout += "?" + urlencode({
            "id_token_hint": oidc_id_token,
            "post_logout_redirect_uri": redirect
        })
    return oidc_logout
