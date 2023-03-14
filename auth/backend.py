from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import Group, User
from mozilla_django_oidc.auth import OIDCAuthenticationBackend

from snpdb.models import UserSettingsOverride


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

        # Copy over basic details from open ID connect
        # Assume there will be no user-name clashes
        user.username = claims.get('preferred_username')
        user.email = claims.get('email', '')
        user.first_name = claims.get('given_name', '')
        user.last_name = claims.get('family_name', '')
        sub = claims.get('sub', None)

        # Work out what groups the user has joined/left since their last login
        django_groups = {g.name for g in user.groups.all()}
        # groups with be in the form of '/variantgrid/some_group_1', '/variantgrid/some_group_2', '/unrelated'
        # convert it so we get 'some_group_1', 'some_group_2'

        user.is_active = True
        if settings.OIDC_REQUIRED_GROUP and settings.OIDC_REQUIRED_GROUP not in claims['groups']:
            user.is_active = False
            user.save()
            messages.add_message(self.request, messages.ERROR, "This account is not authorised for this environment.")
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
        groups = set()
        groups.add(settings.PUBLIC_GROUP_NAME)
        groups.add(settings.LOGGED_IN_USERS_GROUP_NAME)

        # currently only variantgrid permission
        is_super_user = False
        is_bot = False
        for vg in variant_grid_groups:
            permission = '/'.join(vg)
            if permission == 'admin':
                is_super_user = True
            elif permission == 'bot':
                groups.add('variantgrid/bot')
                is_bot = True

        if settings.MAINTENANCE_MODE:
            if 'variantgrid/tester' in groups:
                # testers are allowed to login during maintenance mode
                pass
            elif not is_super_user or is_bot:
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
                #might want to check if valid lab before adding?
                groups.update(assoc_groups)

        removed_groups = django_groups.difference(groups)
        added_groups = groups.difference(django_groups)

        for removed_group in removed_groups:
            # print("removing from group %s" % removed_group)
            group = Group.objects.get(name=removed_group)
            user.groups.remove(group)

        for added_group in added_groups:
            # print("adding to group %s" % added_group)
            group, _ = Group.objects.get_or_create(name=added_group)
            user.groups.add(group)

        # ensures default lab is valid and sets it if there are labs
        # and the default lab is blank
        user_settings_override.auto_set_default_lab()
        user_settings_override.save()
        user.save()
        return user
