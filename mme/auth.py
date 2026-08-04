""" Authenticate inbound MatchMaker Exchange requests.

Peers present the static `X-Auth-Token` we issued them (mme-apis join-protocol). One token
per peer, so every inbound query can be attributed to a named node.
"""
import secrets

from django.conf import settings
from rest_framework.exceptions import AuthenticationFailed


def peer_node_id_for_token(token: str | None) -> str | None:
    """ Return the peer node id whose issued token matches, else None.
        Constant-time compare - `token` is an attacker-supplied bearer secret. """
    if not token:
        return None
    for node_id, expected in (settings.MME_INBOUND_TOKENS or {}).items():
        if expected and secrets.compare_digest(token, expected):
            return node_id
    return None


def authenticate_mme_peer(request) -> str:
    """ Shared guard for every inbound MME endpoint. Returns the peer node id.
        Raises AuthenticationFailed (401) when MME is off or the token is unknown. """
    if not settings.MME_ENABLED:
        raise AuthenticationFailed("MatchMaker Exchange is not enabled")
    peer_node_id = peer_node_id_for_token(request.headers.get("X-Auth-Token"))
    if peer_node_id is None:
        raise AuthenticationFailed("Invalid MME auth token")
    return peer_node_id
