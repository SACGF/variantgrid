import unittest

from django.contrib.auth.models import User

from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.test_data import create_fake_pedigree


class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        cls.user_owner = User.objects.get_or_create(username='testuser')[0]
        cls.user_non_owner = User.objects.get_or_create(username='different_user')[0]
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.pedigree = create_fake_pedigree(cls.user_owner, grch37)
        cls.ped_file = cls.pedigree.ped_file_family.ped_file

        cls.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS = [
            # ('view_pedigree', {"pedigree_id": cls.pedigree.pk}, 200),
            ('view_ped_file', {"ped_file_id": cls.ped_file.pk}, 200),
        ]

        # (url_name, url_kwargs, object to check appears in grid pk column or (grid column, object)
        cls.PRIVATE_GRID_LIST_URLS = [
            ("ped_files_grid", {}, cls.ped_file),
            ("pedigree_grid", {}, cls.pedigree),
        ]

        cls.PRIVATE_AUTOCOMPLETE_URLS = [
            ('pedigree_autocomplete', cls.pedigree, {"q": cls.pedigree.name}),
        ]

    def testUrls(self):
        URL_NAMES_AND_KWARGS = [
            ("pedigrees", {}, 200),
            ("ped_files", {}, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user_non_owner)

    def testPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_owner)

    @prevent_request_warnings
    def testNoPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_non_owner, expected_code_override=403)

    def testAutocompletePermission(self):
        self._test_autocomplete_urls(self.PRIVATE_AUTOCOMPLETE_URLS, self.user_owner, True)

    @prevent_request_warnings
    def testAutocompleteNoPermission(self):
        self._test_autocomplete_urls(self.PRIVATE_AUTOCOMPLETE_URLS, self.user_non_owner, False)

    def testGridListPermission(self):
        self._test_grid_list_urls(self.PRIVATE_GRID_LIST_URLS, self.user_owner, True)

    @prevent_request_warnings
    def testGridListNoPermission(self):
        self._test_grid_list_urls(self.PRIVATE_GRID_LIST_URLS, self.user_non_owner, False)

if __name__ == "__main__":
    unittest.main()
