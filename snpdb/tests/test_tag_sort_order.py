from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from snpdb.models import Tag, TagColor, TagColorsCollection
from snpdb.utils import get_all_tags_and_user_colors

TEST_TAG_IDS = {"aaaSortTest", "bbbSortTest", "cccSortTest"}


class TagSortOrderTest(TestCase):
    """ Custom tag ordering from the tag colors page - issue #343 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("tag_sort_order_user")
        cls.tag_a = Tag.objects.create(pk="aaaSortTest")
        cls.tag_b = Tag.objects.create(pk="bbbSortTest")
        cls.tag_c = Tag.objects.create(pk="cccSortTest")
        cls.collection = TagColorsCollection.objects.create(name="test sort colors", user=cls.user)

    def _test_tag_ids_in_order(self):
        user_tag_colors = get_all_tags_and_user_colors(self.user, self.collection)
        return [tag.pk for tag in user_tag_colors if tag.pk in TEST_TAG_IDS]

    def test_default_order_is_by_tag_name(self):
        self.assertEqual(self._test_tag_ids_in_order(), ["aaaSortTest", "bbbSortTest", "cccSortTest"])

    def test_custom_sort_order_with_unset_tags_first(self):
        # Tags without sort_order sort as 0, so a positive value pushes a tag past them
        TagColor.objects.create(collection=self.collection, tag=self.tag_a, rgb="", sort_order=10)
        self.assertEqual(self._test_tag_ids_in_order(), ["bbbSortTest", "cccSortTest", "aaaSortTest"])

    def test_sort_order_only_rows_emit_no_color(self):
        TagColor.objects.create(collection=self.collection, tag=self.tag_a, rgb="", sort_order=10)
        self.assertEqual(self.collection.get_user_colors_by_tag(), {})

    def test_save_tag_order_via_view(self):
        self.client.force_login(self.user)
        url = reverse("view_tag_colors_collection", kwargs={"tag_colors_collection_id": self.collection.pk})
        response = self.client.post(url, {"tag_order": "cccSortTest,aaaSortTest,bbbSortTest,notARealTag"})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(self.collection.get_sort_order_by_tag(),
                         {"cccSortTest": 0, "aaaSortTest": 1, "bbbSortTest": 2})
        self.assertEqual(self._test_tag_ids_in_order(), ["cccSortTest", "aaaSortTest", "bbbSortTest"])

    def test_save_name_via_view(self):
        self.client.force_login(self.user)
        url = reverse("view_tag_colors_collection", kwargs={"tag_colors_collection_id": self.collection.pk})
        self.client.post(url, {"name": "renamed colors"})
        self.collection.refresh_from_db()
        self.assertEqual(self.collection.name, "renamed colors")
