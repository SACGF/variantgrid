import json
import urllib
from typing import List, Optional, Dict

import requests
from django.conf import settings
from django.contrib.auth.models import User

from library.constants import MINUTE_SECS
from library.oauth import ServerAuth
from snpdb.models.models import Lab
from snpdb.models.models_user_settings import UserSettings


class KeycloakError(Exception):
    pass


class KeycloakNewUser:

    def __init__(self,
                 first_name: str,
                 last_name: str,
                 email: str,
                 lab: Lab):
        self.first_name = first_name
        self.last_name = last_name
        self.email = email
        self.lab = lab
        self.welcome_email = None


class Keycloak:

    def __init__(self, connector: ServerAuth = None):
        if not connector:
            connector = ServerAuth.keycloak_connector()
        self.connector = connector
        self.realm = settings.KEY_CLOAK_REALM

    def change_password(self, user: User):
        user_settings = UserSettings.get_for_user(user)
        if not user_settings.oauth_sub:
            raise ValueError('No sub recorded against this user, cannot change password')

        user_id = user_settings.oauth_sub

        response = requests.put(
            auth=self.connector.auth(),
            url=self.connector.url(f'/auth/admin/realms/{self.realm}/users/{user_id}/execute-actions-email'),
            json=['UPDATE_PASSWORD'],
            timeout=MINUTE_SECS,
        )

        response.raise_for_status()

    def check_groups(self, groups: List[str]) -> List[str]:
        groups.sort()
        response = requests.get(
            auth=self.connector.auth(),
            url=self.connector.url(f'/auth/admin/realms/{self.realm}/groups'),
            timeout=MINUTE_SECS,
        )
        group_array = json.loads(response.text)
        group_dict = {}

        def recurse_subgroups(groups: List[dict]):
            if groups:
                for g in groups:
                    group_dict[g.get('path')] = g.get('id')
                    recurse_subgroups(g.get('subGroups'))

        recurse_subgroups(group_array)

        missing_groups = []
        for group_path in groups:
            if group_path not in group_dict:
                missing_groups.append(group_path)

        for missing_group in groups:
            parts = missing_group.split('/')
            parent = None
            for i in range(1, len(parts)):
                combined = '/'.join(parts[0:i + 1])
                existing_id = group_dict.get(combined)
                if existing_id:
                    parent = existing_id
                else:
                    if not parent:
                        print(group_dict)
                        raise ValueError(f'No parent group found for {missing_group} ({combined})')
                    response = requests.post(
                        auth=self.connector.auth(),
                        url=self.connector.url(f'/auth/admin/realms/{self.realm}/groups/{parent}/children'),
                        json={
                            "name": combined.rsplit('/', maxsplit=1)[-1],
                        },
                        timeout=MINUTE_SECS,
                    )
                    response.raise_for_status()
                    new_group_json = json.loads(response.text)
                    parent = new_group_json.get('id')
                    group_dict[new_group_json.get('path')] = parent

        return [group_dict[gp] for gp in groups]

    def existing_user(self, email: str) -> Optional[Dict]:
        params = {'email': email}
        params_str = urllib.parse.urlencode(params)
        response = requests.get(
            auth=self.connector.auth(),
            url=self.connector.url(f'/auth/admin/realms/{self.realm}/users?{params_str}'),
            timeout=MINUTE_SECS,
        )
        response.raise_for_status()

        user_json = json.loads(response.text)
        if user_json:
            return user_json[0]
        return None

    def welcome_user(self, user: KeycloakNewUser):
        pass

    def add_user(self, user: KeycloakNewUser) -> str:
        # TODO check to see if the user already exists
        if self.existing_user(user.email):
            raise KeycloakError(f'User with email "{user.email}" already exists. Login to Keycloak to see what groups they are associated with')

        lab = user.lab
        groups = [f'/associations/{lab.organization.group_name}', f'/associations/{lab.group_name}']
        if settings.OIDC_REQUIRED_GROUP:
            groups.append(settings.OIDC_REQUIRED_GROUP)

        group_ids = self.check_groups(groups)

        user_rep = {
            'email': user.email,
            'username': user.email,
            'firstName': user.first_name,
            'lastName': user.last_name,
            'enabled': True
        }

        response = requests.post(
            auth=self.connector.auth(),
            url=self.connector.url(f'/auth/admin/realms/{self.realm}/users'),
            json=user_rep,
            timeout=MINUTE_SECS,
        )
        response.raise_for_status()

        # create user just returns 200 success, not details about the user
        # have to search for the user manually
        user_id = self.existing_user(user.email).get('id')

        # add groups
        for group_id in group_ids:
            response = requests.put(
                auth=self.connector.auth(),
                url=self.connector.url(f'/auth/admin/realms/{self.realm}/users/{user_id}/groups/{group_id}'),
                timeout=MINUTE_SECS,
            )
            response.raise_for_status()
            print(response.text)

        if user.welcome_email:
            user.welcome_email.send()

        # send password reset email (might need to double check how long this stays open for)
        response = requests.put(
            auth=self.connector.auth(),
            url=self.connector.url(f'/auth/admin/realms/{self.realm}/users/{user_id}/execute-actions-email'),
            json=['UPDATE_PASSWORD'],
            timeout=MINUTE_SECS,
        )

        response.raise_for_status()

        return response.text
