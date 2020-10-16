import logging
from typing import Dict, Optional, List, Any, Tuple
import json
import os


def _get_env_variable(key: str) -> Tuple[Any, bool]:
    try:
        e_value = os.environ[key]
        if e_value:
            if e_value == 'false':
                e_value = False
            elif e_value == 'true':
                e_value = True
            elif e_value.startswith('"') and e_value.endswith('"'):
                e_value = e_value[1:-1]
            return e_value, True
    except KeyError:
        return None, False


def _load_settings():
    filename = '/etc/variantgrid/settings_config.json'
    try:
        value, found = _get_env_variable('SETTINGS_CONFIG')
        if found:
            filename = value

        logging.info(f'Loading config from {filename}')
        with open(filename) as f:
            return json.load(f)
    except Exception as e:
        logging.info(f"Could not load settings_config from {filename} : {str(e)}\n")
        return {}


_settings_json = _load_settings()

# CELERY.broker_url: guest:guest default account see: https://www.rabbitmq.com/access-control.html
# TODO: migrate CLINGEN_ALLELE_REGISTRY to each environment's settings_config.json
_default_settings = {
    "DB": {
        "host": "localhost",
        "name": "snpdb",
        "user": "snpdb",
        "password": "snpdb",
        "port": ""
    },
    "CELERY": {
        "broker_url": "amqp://guest:guest@localhost"
    },
    "CLINGEN_ALLELE_REGISTRY": {
        "login": "",
        "password": ""
    },
    "ROLLBAR": {
        "access_token": "",
        "client_access_token": ""
    },
    "ENTREZ": {
        "api_key": None,
        "email": None
    }
}
"""
    "AWS": {
        "SES": {
            "region": e.g. "eu-west-1",
            "access_key": None,
            "key_id": None
        }
    },
    "OIDC": {
        "client_secret": as assigned by keycloak or similar,
    },
    # Used for SAPath's Helix
    "SAPATH": {
        "HELIX": {
            "user": None,
            "password": None
        }
    },
    # Used to connect to Keycloak for admin powers (not for authentication)
    "KEYCLOAK": {
        "username": None,
        "password": None,
        "host": e.g. "https://shariant.org.au",
        "oauth_url": e.g. "https://shariant.org.au/auth/realms/master/protocol/openid-connect/token",
        "client_id": e.g. "admin-cli"
    }
"""


def _get_nested(parts: List[str], value) -> Tuple[Any, bool]:
    if not len(parts):
        return value, True
    elif isinstance(value, dict) and parts[0] in value:
        return _get_nested(parts[1:], value[parts[0]])
    else:
        return None, False


def get_secret(key: str) -> Optional[Any]:
    value, found = _get_env_variable(key)
    if found:
        return value

    parts = key.split('.')
    value, found = _get_nested(parts, _settings_json)
    if found:
        #print(f'{key} = {value}')
        return value

    value, found = _get_nested(parts, _default_settings)
    if found:
        if key != 'DB.port':
            logging.warning(f'Warning {key} was loaded from default_settings, migrate to settings config or env variable')
        return value

    raise ValueError(f'No config value found for {key}')


def get_secrets(prefix: str, leafs: List[str]) -> Dict[str, Any]:
    secret_dict = {}
    for leaf in leafs:
        secret_dict[leaf] = get_secret(prefix + '.' + leaf)
    return secret_dict
