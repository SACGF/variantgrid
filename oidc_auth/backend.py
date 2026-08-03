from urllib.parse import urlencode
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import Group, User
from django.core.validators import EmailValidator
from mozilla_django_oidc.auth import OIDCAuthenticationBackend
from mozilla_django_oidc.utils import import_from_settings
from library.log_utils import report_message
from snpdb.models import UserSettingsOverride
from django.core.exceptions import SuspiciousOperation, ValidationError
import re


_email_validator = EmailValidator()
_CONTROL_CHAR_RE = re.compile(r'[\x00-\x1f\x7f]')


def _sanitize_claim_str(value, max_length):
    """Strip control characters and enforce max length on a claim string."""
    return _CONTROL_CHAR_RE.sub('', value)[:max_length]


class VariantGridOIDCAuthenticationBackend(OIDCAuthenticationBackend):

    def create_user(self, claims):
        user = super().create_user(claims)
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

    def create_or_update(self, user: User, claims):

        missing_claims = set()
        for claim in ('preferred_username', 'email', 'sub', 'groups'):
            if not claims.get(claim):
                missing_claims.add(claim)
        if missing_claims:
            report_message(f"Missing claims {missing_claims}", level='error')

        email = claims.get('email', '')
        try:
            _email_validator(email)
        except ValidationError:
            raise SuspiciousOperation(f"Invalid email in OIDC claims: {email!r}")
        email = _sanitize_claim_str(email, 254)

        # Copy over basic details from open ID connect
        # Assume there will be no user-name clashes
        user.username = _sanitize_claim_str(claims.get('preferred_username'), 254)
        user.email = email
        user.first_name = _sanitize_claim_str(claims.get('given_name', ''), 150)
        user.last_name = _sanitize_claim_str(claims.get('family_name', ''), 150)
        sub = claims.get('sub', None)

        # Work out what groups the user has joined/left since their last login
        django_groups = {g.name for g in user.groups.all()}
        # groups with be in the form of '/variantgrid/some_group_1', '/variantgrid/some_group_2', '/unrelated'
        # convert it so we get 'some_group_1', 'some_group_2'

        user.is_active = True
        all_claim_groups = claims.get("groups")
        if settings.OIDC_REQUIRED_GROUP and settings.OIDC_REQUIRED_GROUP not in all_claim_groups:
            user.is_active = False
            user.save()

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
            message = f"This account <i>{user.email}</i> is not authorised for this environment."
            for allowed_environment in allowed_environment_list:
                message += "<br/>" + allowed_environment

            messages.add_message(self.request, messages.ERROR, message, extra_tags="html")
            return user

        oauth_groups = [g.split('/')[1:] for g in claims['groups']]

        associations = [g[1:] for g in oauth_groups if len(g) > 1 and g[0] == 'associations']
        # Should make 'variantgrid' setting configurable at some point in case there are 2 variantgrid installations on the same OAuth instance
        variant_grid_groups = [g[1:] for g in oauth_groups if len(g) > 1 and g[0] == 'variantgrid']

        if not associations and not variant_grid_groups:
            user.is_active = False
            user.save()
            messages.add_message(self.request, messages.ERROR, "This account doesn't belong to any labs.")
            return user
            # raise ValueError(f'Could not find any valid groups in {str_groups}')

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
                messages.add_message(self.request, messages.ERROR, "Non-administrator logins have temporary been disabled.")
                return None

        user.is_superuser = is_super_user
        user.is_staff = is_super_user

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
            # Group may have been deleted out-of-band between logins - skip if gone
            if group := Group.objects.filter(name=removed_group).first():
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
        user.save()
        return user


def provider_logout(request) -> str:
    oidc_logout = import_from_settings("KEY_CLOAK_PROTOCOL_BASE", "") + "/logout"
    if redirect := import_from_settings("LOGOUT_REDIRECT_URL", ""):
        if oidc_id_token := request.session["oidc_id_token"]:
            oidc_logout += "?" + urlencode({
                "id_token_hint": oidc_id_token,
                "post_logout_redirect_uri": redirect
            })
    return oidc_logout
