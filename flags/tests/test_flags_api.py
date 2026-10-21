"""FlagsView / FlagView: bad parameters answer 400, other users' labs stay hidden, and collection permissions apply."""
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from flags.models import FlagCollection
from flags.models.enums import FlagStatus
from flags.models.models import Flag, FlagResolution, FlagType, FlagTypeContext, FlagTypeResolution
from snpdb.models.models_variant import Allele


class FlagsApiTest(TestCase):

    def setUp(self):
        self.user = User.objects.create_superuser("flag_api_admin", "admin@example.com", "pw")
        self.other = User.objects.create_user("flag_api_other", "other@example.com", "pw")
        self.client.force_login(self.user)
        self.fc = FlagCollection.objects.create(context=FlagTypeContext.objects.get(pk="allele"))
        self.url = reverse("flags_api", kwargs={"flag_collection_id": self.fc.pk})
        self.flag_type = FlagType.objects.create(id="flag_api_test", context_id="allele", label="Test",
                                                 description="Test", raise_permission="U", permission="U")
        open_resolution = FlagResolution.objects.filter(status=FlagStatus.OPEN).first()
        FlagTypeResolution.objects.create(flag_type=self.flag_type, resolution=open_resolution)

    def _post(self, url, data):
        return self.client.post(url, data, content_type="application/json")

    def test_invalid_history_and_since(self):
        for params in [{"history": "abc"}, {"since": "xyz"}, {"since": "1e20"}]:
            with self.subTest(params=params):
                self.assertEqual(self.client.get(self.url, params).status_code, 400)

    def test_post_without_flag_type(self):
        """ The retired 'watch' POST has no flag_type """
        response = self._post(self.url, {"watch": True})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json(), {"error": "Invalid flag data"})

    def test_post_invalid_flag_type_or_resolution(self):
        for data in [{"flag_type": "no_such_type"}, {"flag_type": self.flag_type.pk, "resolution": "no_such"}]:
            with self.subTest(data=data):
                self.assertEqual(self._post(self.url, data).status_code, 400)
        self.assertFalse(Flag.objects.filter(collection=self.fc).exists())

    def test_flag_invalid_resolution(self):
        self.assertEqual(self._post(self.url, {"flag_type": self.flag_type.pk}).status_code, 200)
        flag = Flag.objects.get(collection=self.fc)
        flag_url = reverse("flag_api", kwargs={"flag_id": flag.pk})
        self.assertEqual(self._post(flag_url, {"resolution": "no_such"}).status_code, 400)

    def test_lab_only_for_requesting_user(self):
        self._post(self.url, {"flag_type": self.flag_type.pk})
        Flag.objects.filter(collection=self.fc).update(user=self.other)
        users = {u["id"]: u for u in self.client.get(self.url).json()["users"]}
        self.assertIn("lab", users[self.user.pk])
        self.assertNotIn("lab", users[self.other.pk])

    def test_no_perm_user(self):
        """ The collection has no source object, so a non-superuser has NO_PERM on it """
        self._post(self.url, {"flag_type": self.flag_type.pk})
        flag_url = reverse("flag_api", kwargs={"flag_id": Flag.objects.get(collection=self.fc).pk})
        self.client.force_login(self.other)

        data = self.client.get(self.url).json()
        self.assertEqual(data["collections"], [])
        self.assertEqual(data["flags"], [])
        self.assertEqual(self.client.get(self.url, {"history": self.fc.pk}).status_code, 403)
        self.assertEqual(self._post(self.url, {"flag_type": self.flag_type.pk}).status_code, 403)
        self.assertEqual(self.client.get(flag_url).status_code, 403)
        self.assertEqual(self._post(flag_url, {"comment": "hello"}).status_code, 403)
        self.assertEqual(Flag.objects.filter(collection=self.fc).count(), 1)

    def test_raise_permission(self):
        """ Allele collections give every user USERS, which can't raise an ADMIN-only flag type """
        allele_fc = Allele.objects.create().flag_collection_safe
        url = reverse("flags_api", kwargs={"flag_collection_id": allele_fc.pk})
        admin_type = FlagType.objects.create(id="flag_api_test_admin", context_id="allele", label="Admin",
                                             description="Admin", raise_permission="X", permission="X")
        self.client.force_login(self.other)
        self.assertEqual(self._post(url, {"flag_type": self.flag_type.pk}).status_code, 200)
        self.assertEqual(self._post(url, {"flag_type": admin_type.pk}).status_code, 403)
        self.assertEqual(self.client.get(url, {"history": self.fc.pk}).status_code, 403)
