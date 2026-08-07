from types import SimpleNamespace

from django.test import TestCase, override_settings
from rest_framework.exceptions import AuthenticationFailed

from mme.auth import peer_node_id_for_token, authenticate_mme_peer

TOKENS = {"genematcher": "gm-secret", "decipher": "dc-secret"}


def _request(token=None) -> SimpleNamespace:
    return SimpleNamespace(headers={"X-Auth-Token": token} if token is not None else {})


@override_settings(MME_ENABLED=True, MME_INBOUND_TOKENS=TOKENS)
class PeerTokenTestCase(TestCase):

    def test_each_token_resolves_to_its_own_node(self):
        self.assertEqual(peer_node_id_for_token("gm-secret"), "genematcher")
        self.assertEqual(peer_node_id_for_token("dc-secret"), "decipher")

    def test_unknown_token_resolves_to_nothing(self):
        self.assertIsNone(peer_node_id_for_token("some-other-secret"))

    def test_blank_and_missing_tokens_resolve_to_nothing(self):
        self.assertIsNone(peer_node_id_for_token(""))
        self.assertIsNone(peer_node_id_for_token(None))

    @override_settings(MME_INBOUND_TOKENS={"unissued": None, "alsounissued": ""})
    def test_node_configured_without_a_token_never_authenticates(self):
        # A peer we have not yet issued a token to must not be matchable by a blank token.
        self.assertIsNone(peer_node_id_for_token(""))
        self.assertIsNone(peer_node_id_for_token("anything"))


@override_settings(MME_ENABLED=True, MME_INBOUND_TOKENS=TOKENS)
class AuthenticateMMEPeerTestCase(TestCase):

    def test_returns_peer_node_id(self):
        self.assertEqual(authenticate_mme_peer(_request("dc-secret")), "decipher")

    def test_unknown_token_raises(self):
        with self.assertRaises(AuthenticationFailed):
            authenticate_mme_peer(_request("nope"))

    @override_settings(MME_ENABLED=False)
    def test_disabled_deployment_raises_even_with_a_valid_token(self):
        with self.assertRaises(AuthenticationFailed):
            authenticate_mme_peer(_request("gm-secret"))
