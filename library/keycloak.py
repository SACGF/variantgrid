import json
import logging
import urllib
from typing import Optional, Callable

import requests
from django.conf import settings
from django.contrib.auth.models import User
from django.template.loader import render_to_string
from requests import Response

from email_manager.models import EmailLog
from library.constants import MINUTE_SECS
from library.email import Email
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
        self.welcome_email: Email = None


class Keycloak:

    def __init__(self, connector: Optional[ServerAuth] = None):
        if not connector:
            connector = ServerAuth.keycloak_connector()
        self.connector = connector
        self.realm = settings.KEY_CLOAK_REALM

    def ping(self):
        # Method does a basic request to KeyCloak, should raise an exception if anything goes wrong
        # returns nothing otherwise
        self.request(
            'GET',
            url=f'/admin/realms/{self.realm}/clients',
        )

    def change_password(self, user: User):
        user_settings = UserSettings.get_for_user(user)
        if not user_settings.oauth_sub:
            raise ValueError('No sub recorded against this user, cannot change password')

        user_id = user_settings.oauth_sub

        self.request(
            'PUT',
            url=f'/admin/realms/{self.realm}/users/{user_id}/execute-actions-email',
            json_data=['UPDATE_PASSWORD']
        )

    def check_groups(self, groups: list[str]) -> list[str]:
        """
        Makes sure the groups exist within KeyCloak and returns their IDs
        :param groups: list of group names
        """
        groups.sort()
        # Get ALL groups within KeyCloak, so we can see if the groups we want exist in there yet
        response = self.request(
            'GET',
            url=f'/admin/realms/{self.realm}/groups?q=*',
        )
        group_array = json.loads(response.text)
        group_dict = {}

        def recurse_subgroups(sub_groups: list[dict]):
            if sub_groups:
                for g in sub_groups:
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
                        logging.error("No parent group found for %s (%s), group_dict: %s", missing_group, combined, group_dict)
                        raise ValueError(f'No parent group found for {missing_group} ({combined})')
                    response = self.request(
                        'POST',
                        url=f'/admin/realms/{self.realm}/groups/{parent}/children',
                        json_data={
                            "name": combined.rsplit('/', maxsplit=1)[-1],
                        }
                    )
                    new_group_json = json.loads(response.text)
                    parent = new_group_json.get('id')
                    group_dict[new_group_json.get('path')] = parent

        return [group_dict[gp] for gp in groups]

    def request(self, method: str, url: str, json_data=None) -> Response:
        # :param method: method for the new :class:`Request` object: ``GET``, ``OPTIONS``, ``HEAD``, ``POST``, ``PUT``, ``PATCH``, or ``DELETE``.
        response = requests.request(
            method,
            auth=self.connector.auth,
            url=self.connector.url(url),
            json=json_data,
            timeout=MINUTE_SECS
        )
        response.raise_for_status()
        return response

    def existing_user(self, email: str) -> Optional[dict]:
        params = {'email': email}
        params_str = urllib.parse.urlencode(params)
        response = self.request(
            'GET',
            url=f'/admin/realms/{self.realm}/users?{params_str}',
        )
        if json_response := response.json():
            return json_response[0]
        else:
            return None

    def welcome_user(self, user: KeycloakNewUser):
        pass

    def add_user(self, user: KeycloakNewUser, pre_password_reset: Optional[Callable] = None) -> str:
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

        self.request(
            'POST',
            url=f'/admin/realms/{self.realm}/users',
            json_data=user_rep
        )

        # create user just returns 200 success, not details about the user
        # have to search for the user manually
        user_id = self.existing_user(user.email).get('id')

        # add groups
        for group_id in group_ids:
            response = self.request(
                'PUT',
                url=f'/admin/realms/{self.realm}/users/{user_id}/groups/{group_id}',
            )
            logging.info("Keycloak add user to group response: %s", response.text)

        if pre_password_reset:
            pre_password_reset()

        # send password reset email (might need to double check how long this stays open for)
        response = self.request(
            'PUT',
            url=f'/admin/realms/{self.realm}/users/{user_id}/execute-actions-email',
            json_data=['UPDATE_PASSWORD']
        )

        return response.text
