import unittest

from django.contrib.auth.models import User

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from library.django_utils.unittest_utils import prevent_request_warnings, URLTestCase
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.models.models_cohort import Cohort
from snpdb.models.models_columns import CustomColumnsCollection
from snpdb.models.models_enums import ImportStatus
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_genomic_interval import GenomicIntervalsCollection, GenomicIntervalsCategory
from snpdb.tests.test_data import create_fake_trio


class Test(URLTestCase):
    """ Need to override settings as ManifestStaticFilesStorage expects staticfiles.json to exist
        and contain the file asked. @see https://stackoverflow.com/a/51580328/295724
    """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)

        owner_username = f"test_user_{__file__}_owner"
        non_owner_username = f"test_user_{__file__}_non_owner"
        cls.user_owner = User.objects.get_or_create(username=owner_username)[0]
        cls.user_non_owner = User.objects.get_or_create(username=non_owner_username)[0]

        cls.trio = create_fake_trio(cls.user_owner, grch37)
        cls.cohort = cls.trio.cohort
        cls.vcf = cls.cohort.vcf
        cls.sample = cls.vcf.sample_set.first()

        # Auto cohorts don't show on list
        cls.cohort2 = Cohort.objects.create(name="blah cohort", user=cls.user_owner, vcf=None, genome_build=grch37,
                                            import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.user_owner, cls.cohort2)

        gic = GenomicIntervalsCategory.objects.get_or_create(name="blah", description="blah")[0]
        cls.genomic_intervals_collection = GenomicIntervalsCollection.objects.create(name="",
                                                                                     category=gic,
                                                                                     user=cls.user_owner)
        transcript_version = create_fake_transcript_version(grch37)
        gene_symbol = transcript_version.gene_version.gene_symbol

        cls.custom_columns_collection = CustomColumnsCollection.objects.create(name="blah", user=cls.user_owner)

        cls.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS = [
            ('view_vcf', {"vcf_id": cls.vcf.pk}, 200),
            ('get_patient_upload_csv_for_vcf', {"pk": cls.vcf.pk}, 200),

            # Sample related
            ('view_sample', {"sample_id": cls.sample.pk}, 200),
            ('sample_variants_tab', {"sample_id": cls.sample.pk}, 200),
            ('sample_variants_gene_detail', {"sample_id": cls.sample.pk, "gene_symbol": gene_symbol}, 200),
            ('sample_graphs_tab', {"sample_id": cls.sample.pk}, 200),
            ('sample_permissions_tab', {"sample_id": cls.sample.pk}, 200),
            # Cohort
            ('view_cohort_details_tab', {"cohort_id": cls.cohort.pk}, 200),
            ('view_cohort', {"cohort_id": cls.cohort.pk}, 302),
            ('cohort_hotspot', {"cohort_id": cls.cohort.pk}, 200),
            ('cohort_gene_counts', {"cohort_id": cls.cohort.pk}, 200),
            ('cohort_sort', {"cohort_id": cls.cohort.pk}, 200),
            ('cohort_sample_count', {"cohort_id": cls.cohort.pk}, 200),
            ('cohort_sample_edit', {"cohort_id": cls.cohort.pk}, 200),

            # Trio
            ('view_trio', {"pk": cls.trio.pk}, 200),

            # Data objects
            ('view_genomic_intervals', {"genomic_intervals_collection_id": cls.genomic_intervals_collection.pk}, 200),

            # Grids for objects
            ("cohort_sample_grid", {"cohort_id": cls.cohort.pk, "op": "config"}, 200),
            ("cohort_sample_grid", {"cohort_id": cls.cohort.pk, "op": "handler"}, 200),
        ]

        cls.PRIVATE_AUTOCOMPLETE_URLS = [
            ('vcf_autocomplete', cls.vcf, {"q": cls.vcf.name}),
            ('sample_autocomplete', cls.sample, {"q": cls.sample.name}),
            ('cohort_autocomplete', cls.cohort, {"q": cls.cohort.name}),
            ('trio_autocomplete', cls.trio, {"q": cls.trio.name}),
        ]

        cls.PRIVATE_GRID_LIST_URLS = [
            ("vcfs_grid", {}, cls.vcf),
            ("samples_grid", {}, cls.sample),
            ("cohort_grid", {}, cls.cohort2),
            ("trio_grid", {}, cls.trio),
            ("genomic_intervals_grid", {}, cls.genomic_intervals_collection),
            ("custom_columns_grid", {}, cls.custom_columns_collection),
        ]

    @classmethod
    def tearDownClass(cls):
        cls.user_owner.delete()
        cls.user_non_owner.delete()
        cls.vcf.delete()  # Will cascade sample/cohort/trio
        cls.genomic_intervals_collection.delete()

        super().tearDownClass()

    def testUrls(self):
        URL_NAMES_AND_KWARGS = [
            ("data", {}, 200),
            ("vcfs", {}, 200),
            ("samples", {}, 200),
            ("bed_files", {}, 200),
            ("cohorts", {}, 200),
            ("trios", {}, 200),
            ("manual_variant_entry", {}, 200),
            ("custom_columns", {}, 200),
            ("tag_settings", {}, 200),
            ("view_user_settings", {}, 200),
            ("labs", {}, 200),
            ("index", {}, 200),
            ("variant_tags", {}, 200),
        ]

        self._test_urls(URL_NAMES_AND_KWARGS, self.user_non_owner)

        # Make sure that GlobalLoginRequiredMiddleware bounces unauth users
        self._test_urls(URL_NAMES_AND_KWARGS, expected_code_override=302)

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
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
