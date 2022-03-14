import unittest

from django.contrib.auth.models import User

from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings
from library.guardian_utils import assign_permission_to_user_and_groups
from pedigree.models import PedFile, Pedigree, PedFileFamily
from snpdb.models import ImportStatus, Cohort
from snpdb.models.models_genome import GenomeBuild


class Test(URLTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls.user_owner = User.objects.get_or_create(username='testuser')[0]
        cls.user_non_owner = User.objects.get_or_create(username='different_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = Cohort.objects.get_or_create(name="fake cohort",
                                                  user=cls.user_owner,
                                                  version=1,
                                                  import_status=ImportStatus.SUCCESS,
                                                  genome_build=cls.grch37)[0]
        cls.ped_file = PedFile.objects.get_or_create(name="fakepf", user=cls.user_owner,
                                                     import_status=ImportStatus.SUCCESS)[0]
        assign_permission_to_user_and_groups(cls.user_owner, cls.ped_file)
        cls.ped_file_family = PedFileFamily.objects.get_or_create(name="fake family", ped_file=cls.ped_file)[0]
        cls.pedigree = Pedigree.objects.get_or_create(user=cls.user_owner, name="fake pedigree",
                                                      cohort=cls.cohort, ped_file_family=cls.ped_file_family)[0]

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
