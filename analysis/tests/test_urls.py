import os

from django.conf import settings
from django.contrib.auth.models import User
from django.test.client import Client
from django.urls.base import reverse
from django.utils import timezone
import unittest

from analysis.models import Analysis, SampleNode, KaryomappingAnalysis
from analysis.models.enums import SNPMatrix
from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils.unittest_utils import prevent_request_warnings, URLTestCase
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.models import Variant
from snpdb.models.models_cohort import Cohort, Trio, CohortSample, CohortGenotypeCollection, CohortGenotype
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_vcf import VCF, Sample
from snpdb.models.models_enums import ImportStatus
from snpdb.tests.utils.vcf_testing_utils import slowly_create_loci_and_variants_for_vcf


class Test(URLTestCase):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        annotation_version = get_fake_annotation_version(grch37)  # Needed in cohort_hotspot_graph

        cls.user_owner = User.objects.get_or_create(username='testuser')[0]
        cls.user_non_owner = User.objects.get_or_create(username='different_user')[0]
        cls.vcf = VCF.objects.create(name="test_urls_vcf", genotype_samples=1, genome_build=grch37,
                                     import_status=ImportStatus.SUCCESS,
                                     user=cls.user_owner, date=timezone.now())
        cls.sample = Sample.objects.create(name="sample1", vcf=cls.vcf, import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.user_owner, cls.vcf)
        assign_permission_to_user_and_groups(cls.user_owner, cls.sample)

        mother_sample = Sample.objects.create(name="mother", vcf=cls.vcf)
        father_sample = Sample.objects.create(name="father", vcf=cls.vcf)
        cls.cohort = Cohort.objects.create(name="test_urls_cohort", vcf=cls.vcf, genome_build=grch37,
                                           import_status=ImportStatus.SUCCESS)

        proband_cs = CohortSample.objects.create(cohort=cls.cohort, sample=cls.sample,
                                                 cohort_genotype_packed_field_index=0, sort_order=1)
        mother_cs = CohortSample.objects.create(cohort=cls.cohort, sample=mother_sample,
                                                cohort_genotype_packed_field_index=1, sort_order=2)
        father_cs = CohortSample.objects.create(cohort=cls.cohort, sample=father_sample,
                                                cohort_genotype_packed_field_index=2, sort_order=3)

        assign_permission_to_user_and_groups(cls.user_owner, cls.cohort)

        # Cohort version has been bumped every time a cohort sample has been added
        collection = CohortGenotypeCollection.objects.create(cohort=cls.cohort,
                                                             cohort_version=cls.cohort.version,
                                                             num_samples=cls.cohort.cohortsample_set.count())

        cls.trio = Trio.objects.create(name="test_urls_trio",
                                       cohort=cls.cohort,
                                       mother=mother_cs,
                                       mother_affected=True,
                                       father=father_cs,
                                       father_affected=False,
                                       proband=proband_cs)
        vcf_filename = os.path.join(settings.BASE_DIR, "annotation/tests/test_data/test_grch37.vep_annotated.vcf")
        slowly_create_loci_and_variants_for_vcf(grch37, vcf_filename, get_variant_id_from_info=True)
        variant = Variant.objects.filter(Variant.get_no_reference_q()).first()
        CohortGenotype.objects.create(collection=collection,
                                      variant=variant,
                                      ref_count=1,
                                      het_count=1,
                                      hom_count=1,
                                      filters="&",
                                      samples_zygosity="ERO",
                                      samples_allele_depth=[42, 22, 32],
                                      samples_allele_frequency=[100, 100, 100],
                                      samples_read_depth=[42, 22, 32],
                                      samples_genotype_quality=[20, 20, 20],
                                      samples_phred_likelihood=[0, 0, 0])

        # Auto cohorts don't show on list
        cls.cohort2 = Cohort.objects.create(name="blah cohort", vcf=None, genome_build=grch37,
                                            import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.user_owner, cls.cohort2)

        cls.analysis = Analysis.objects.create(genome_build=grch37,
                                               annotation_version=annotation_version,
                                               user=cls.user_owner)
        cls.node = SampleNode.objects.create(analysis=cls.analysis,
                                             sample=cls.sample,
                                             count=1)

        cls.karyomapping_analysis = KaryomappingAnalysis.objects.create(user=cls.user_owner,
                                                                        name="test karyomapping",
                                                                        trio=cls.trio)

        grid_column_name = "variantannotation__transcript_version__gene_version__gene_symbol"
        analysis_params = {"analysis_id": cls.analysis.pk}
        node_version_params = {"node_id": cls.node.pk, "node_version": cls.node.version}
        analysis_version_and_node_version_params = {**node_version_params,
                                                    "analysis_version": cls.analysis.version,
                                                    "extra_filters": "default"}

        cls.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS = [
            ('analysis', {"analysis_id": cls.analysis.pk}, 200),

            # analysis templates - TODO: Need to make Template and verify it is kept private
            # ('analysis_templates_list', {"analysis_template_id": cls.analysis.pk}, 200),

            # Node editor
            ('node_view', analysis_version_and_node_version_params, 200),
            ('node_debug', analysis_version_and_node_version_params, 200),

            ('node_doc', {"node_id": cls.node.pk}, 200),
            ('node_load', {"node_id": cls.node.pk}, 302),

            ('node_column_summary', {**node_version_params,
                                     "analysis_version": cls.analysis.version,
                                     "extra_filters": "default",
                                     "grid_column_name": grid_column_name,
                                     "significant_figures": 2}, 200),
            ('node_snp_matrix', {**node_version_params,
                                 "conversion": SNPMatrix.TOTAL_PERCENT,
                                 "significant_figures": 2}, 200),

            ('node_data', {"node_id": cls.node.pk}, 200),
            ('analysis_node_versions', analysis_params, 200),
            ('analysis_editor_and_grid', analysis_params, 200),
            ('analysis_settings', analysis_params, 200),
            ('analysis_settings_details_tab', analysis_params, 200),
            ('analysis_settings_node_counts_tab', analysis_params, 200),
            ('analysis_input_samples', analysis_params, 200),

            # Node data
            ('node_data_grid', analysis_version_and_node_version_params, 200),
            ('node_async_wait', analysis_version_and_node_version_params, 200),
            ('node_errors', analysis_version_and_node_version_params, 200),
            ('node_method_description', node_version_params, 200),

            ('view_karyomapping_analysis', {"pk": cls.karyomapping_analysis.pk}, 200),
        ]

        cls.PRIVATE_GRID_LIST_URLS = [
            #("vcfs_grid", {}, cls.vcf),
        ]

    @classmethod
    def tearDownClass(cls):
        cls.user_owner.delete()
        cls.user_non_owner.delete()
        cls.vcf.delete()  # Will cascade sample/cohort/trio

        super().tearDownClass()

    def testUrls(self):
        URL_NAMES_AND_KWARGS = [
            ("analyses", {}, 200),
            ("analysis_templates", {}, 200),
            ("analyses_variant_tags", {}, 200),
            ("karyomapping_analyses", {}, 200),
        ]

        self._test_urls(URL_NAMES_AND_KWARGS, self.user_non_owner)

    def testPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_owner)

    @prevent_request_warnings
    def testNoPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_non_owner, expected_code_override=403)

    def testGridListPermission(self):
        self._test_grid_list_urls(self.PRIVATE_GRID_LIST_URLS, self.user_owner, True)

    @prevent_request_warnings
    def testGridListNoPermission(self):
        self._test_grid_list_urls(self.PRIVATE_GRID_LIST_URLS, self.user_non_owner, False)

    def _testVariantGridExport(self, export_type: str):
        client = Client()
        client.force_login(self.user_owner)
        url = reverse("node_grid_export")
        params = {"node_id": self.node.pk,
                  "version_id": self.node.version,
                  "extra_filters": "",
                  "export_type": export_type}
        url = url + "?" + "&".join([f"{k}={v}" for k, v in params.items()])
        response = client.get(url)
        response.getvalue()  # Read streaming content
        self.assertEqual(response.status_code, 200)

    def testVariantGridExport(self):
        for export_type in ["vcf", "csv"]:
            self._testVariantGridExport(export_type)


if __name__ == "__main__":
    unittest.main()
