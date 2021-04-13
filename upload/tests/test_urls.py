from django.contrib.auth.models import User
import unittest

from django.utils.timezone import localdate

from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings
from snpdb.models import ImportSource, ProcessingStatus, VCF
from snpdb.models.models_genome import GenomeBuild
from upload.models import UploadedFile, UploadPipeline, UploadedVCF, UploadedFileTypes


class Test(URLTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls.user_owner = User.objects.get_or_create(username='testuser')[0]
        cls.user_non_owner = User.objects.get_or_create(username='different_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        uploaded_file = UploadedFile.objects.create(user=cls.user_owner,
                                                    name="fake uploaded file",
                                                    path="/tmp/foo.vcf",
                                                    file_type=UploadedFileTypes.VCF,
                                                    import_source=ImportSource.COMMAND_LINE)
        upload_pipeline = UploadPipeline.objects.create(status=ProcessingStatus.PROCESSING,
                                                        uploaded_file=uploaded_file)
        upload_step = upload_pipeline.uploadstep_set.create(name="fake step", sort_order=0)
        vcf = VCF.objects.create(name="fake_vcf.vcf", date=localdate(), user=cls.user_owner,
                                 genotype_samples=0, genome_build=cls.grch37)
        UploadedVCF.objects.create(uploaded_file=uploaded_file,
                                   vcf=vcf,
                                   upload_pipeline=upload_pipeline)

        upload_pipeline_kwargs = {"upload_pipeline_id": upload_pipeline.pk}
        cls.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS = [
            ('view_uploaded_file', {"uploaded_file_id": uploaded_file.pk}, 200),
            ('view_upload_pipeline', upload_pipeline_kwargs, 200),
            ('view_upload_pipeline_warnings_and_errors', upload_pipeline_kwargs, 200),
        ]

        # (url_name, url_kwargs, object to check appears in grid pk column or (grid column, object)
        cls.PRIVATE_GRID_LIST_URLS = [
            ("upload_step_grid", upload_pipeline_kwargs, upload_step),
            ("upload_pipeline_modified_variants_grid", upload_pipeline_kwargs, None),
        ]

    def testUrls(self):
        URL_NAMES_AND_KWARGS = [
            ("upload", {}, 200),
            ("upload_poll", {}, 200),
            ("view_upload_stats", {}, 200),
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


if __name__ == "__main__":
    unittest.main()
