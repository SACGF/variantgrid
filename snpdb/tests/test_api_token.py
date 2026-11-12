from datetime import timedelta
from importlib import import_module

from django.apps import apps
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.urls import reverse
from django.utils import timezone
from rest_framework.authtoken.models import Token
from rest_framework.response import Response
from rest_framework.test import APIRequestFactory
from rest_framework.views import APIView

from eventlog.models import Event
from snpdb.forms import ApiTokenCreateForm
from snpdb.models.models_api_token import ApiToken


class WhoAmIView(APIView):
    """ Authenticates through REST_FRAMEWORK's DEFAULT_AUTHENTICATION_CLASSES, as every API view does """

    def get(self, request):
        return Response({"username": request.user.username})


class ApiTokenAuthenticationTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user("api_token_user")

    def _call(self, key: str):
        request = APIRequestFactory().get("/api/whoami", HTTP_AUTHORIZATION=f"Token {key}")
        return WhoAmIView.as_view()(request)

    def _create_token(self, **kwargs) -> tuple[ApiToken, str]:
        return ApiToken.create_for_user(self.user, "pipeline", kwargs.get("expiry", timedelta(days=30)),
                                        created_by=self.user)

    def test_token_authenticates_and_records_one_use_per_hour(self):
        token, key = self._create_token()
        self.assertTrue(key.startswith("vg_"))
        for _ in range(3):
            response = self._call(key)
            self.assertEqual(200, response.status_code)
            self.assertEqual("api_token_user", response.data["username"])
        self.assertEqual(1, token.uses.count())

        ApiToken.objects.filter(pk=token.pk).update(last_used=timezone.now() - timedelta(hours=2))
        self._call(key)
        self.assertEqual(2, token.uses.count())
        self.assertEqual("/api/whoami", token.uses.first().path)

    def test_revoked_and_expired_tokens_are_refused_and_kept(self):
        revoked, revoked_key = self._create_token()
        revoked.revoke(self.user)
        expired, expired_key = self._create_token(expiry=timedelta(days=-1))
        for key in (revoked_key, expired_key):
            self.assertEqual(401, self._call(key).status_code)
        self.assertEqual(2, ApiToken.objects.filter(pk__in=[revoked.pk, expired.pk]).count())

    def test_imported_legacy_token_authenticates(self):
        legacy = Token.objects.create(user=self.user)
        import_module("snpdb.migrations.0016_api_token")._import_legacy_tokens(apps, None)
        self.assertFalse(Token.objects.exists())
        token = ApiToken.objects.get(user=self.user)
        self.assertEqual(("Legacy token", None, legacy.created), (token.name, token.expiry, token.created))
        self.assertEqual(200, self._call(legacy.key).status_code)



class ApiTokenViewsTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user("api_token_owner")
        cls.other_user = User.objects.create_user("api_token_other")
        cls.admin = User.objects.create_superuser("api_token_admin")

    def test_create_shows_key_once_and_audits(self):
        self.client.force_login(self.user)
        response = self.client.post(reverse("api_tokens"), {"name": "my script", "expiry_days": 30})
        token = ApiToken.objects.get(user=self.user)
        key = response.context["new_token_key"]
        self.assertTrue(token.digest and key not in token.digest)
        self.assertContains(response, key)
        self.assertAlmostEqual(30, (token.expiry - timezone.now()).days, delta=1)
        self.assertTrue(Event.objects.filter(name="api_token_created", user=self.user).exists())

        response = self.client.get(reverse("api_tokens"))
        self.assertNotContains(response, key)

    @override_settings(API_TOKEN_MAX_LIFETIME_DAYS=60)
    def test_expiry_capped_at_max_lifetime(self):
        self.assertFalse(ApiTokenCreateForm({"name": "x", "expiry_days": 90}).is_valid())
        self.assertFalse(ApiTokenCreateForm({"name": "x", "expiry_days": ApiTokenCreateForm.NEVER}).is_valid())
        self.assertTrue(ApiTokenCreateForm({"name": "x", "expiry_days": 60}).is_valid())

    @override_settings(API_TOKEN_MAX_LIFETIME_DAYS=None)
    def test_never_expires_only_without_max_lifetime(self):
        form = ApiTokenCreateForm({"name": "x", "expiry_days": ApiTokenCreateForm.NEVER})
        self.assertTrue(form.is_valid())
        self.assertIsNone(form.expiry)

    def test_revoke_own_not_others(self):
        token, _ = ApiToken.create_for_user(self.user, "mine", timedelta(days=7), created_by=self.user)
        self.client.force_login(self.other_user)
        response = self.client.post(reverse("api_token_revoke"), {"digest": token.digest})
        self.assertEqual(403, response.status_code)
        self.assertEqual(403, self.client.get(reverse("api_token_detail", kwargs={"digest": token.digest})).status_code)

        self.client.force_login(self.user)
        self.client.post(reverse("api_token_revoke"), {"digest": token.digest})
        token.refresh_from_db()
        self.assertEqual(self.user, token.revoked_by)

    def test_admin_revokes_any_token(self):
        token, _ = ApiToken.create_for_user(self.user, "theirs", timedelta(days=7), created_by=self.user)
        self.client.force_login(self.admin)
        response = self.client.post(reverse("api_token_revoke"), {"digest": token.digest, "admin": "1"})
        self.assertRedirects(response, reverse("api_tokens_admin"), fetch_redirect_response=False)
        token.refresh_from_db()
        self.assertEqual(self.admin, token.revoked_by)
        self.assertTrue(Event.objects.filter(name="api_token_revoked", user=self.admin).exists())
