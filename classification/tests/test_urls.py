import unittest

from django.contrib.auth.models import User

from annotation.fake_annotation import get_fake_annotation_version, create_fake_variants
from classification.autopopulate_evidence_keys.autopopulate_evidence_keys import \
    create_classification_for_sample_and_variant_objects
from classification.models import EvidenceKey
from library.django_utils.unittest_utils import prevent_request_warnings, URLTestCase
from snpdb.models import GenomeBuild, Variant, ClinGenAllele, Allele, VariantAllele, AlleleOrigin, Lab, Organization, \
    Country


class Test(URLTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls.user_admin = User.objects.create_superuser("admin_user", "foo@bar", "password")
        cls.user_owner = User.objects.get_or_create(username='testuser')[0]
        cls.user_non_owner = User.objects.get_or_create(username='different_user')[0]

        organization = Organization.objects.get_or_create(name="Fake Org", group_name="fake_org")[0]
        australia = Country.objects.get_or_create(name="Australia")[0]
        cls.lab = Lab.objects.get_or_create(name="Fake Lab", city="Adelaide", country=australia,
                                            organization=organization, group_name="fake_org/fake_lab")[0]
        cls.lab.group.user_set.add(cls.user_owner)

        cls.other_lab = Lab.objects.get_or_create(name="Fake Lab 2", city="Adelaide", country=australia,
                                            organization=organization, group_name="fake_org/fake_lab2")[0]
        cls.other_lab.group.user_set.add(cls.user_non_owner)

        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        annotation_version = get_fake_annotation_version(grch37)
        create_fake_variants(grch37)
        variant = Variant.objects.get(locus__contig__name='13', locus__position=95839002, locus__ref__seq='C',
                                      alt__seq='T')
        clingen_allele = ClinGenAllele.objects.get_or_create(id=7019515, api_response={})[0]
        allele = Allele.objects.get_or_create(clingen_allele=clingen_allele)[0]
        VariantAllele.objects.get_or_create(variant=variant,
                                            genome_build=grch37,
                                            allele=allele,
                                            origin=AlleleOrigin.IMPORTED_TO_DATABASE)

        sample = None
        # ensembl_transcript_accession = "ENST00000376887.4:c.1498G>A"
        classification = create_classification_for_sample_and_variant_objects(cls.user_owner,
                                                                              sample,
                                                                              variant,
                                                                              grch37,
                                                                              annotation_version=annotation_version)
        # Need a modification to be able to export
        classification.patch_value({"clinical_significance": "VUS"}, user=cls.user_owner, save=True)
        classification.publish_latest(cls.user_owner)

        classification_kwargs = {"classification_id": classification.pk}
        #  = {"report_id": None}
        cls.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS = [
            # ('classification_diff', {"GET_PARAMS": {"history": classification.pk}}, 200),
            # ('classification_diff', {"GET_PARAMS": {"record_ids": f"[{classification.pk}]"}}, 200),
            # ('classification_csv', record_kwargs, 200),
            # ('view_template_report', record_kwargs, 200),
            ('classification_history', classification_kwargs, 200),
            ('view_classification', classification_kwargs, 200),
            # ('discordance_report', report_kwargs, 200),
            # ('discordance_export', report_kwargs, 200),
        ]

    def testAdminUrls(self):
        ADMIN_URL_NAMES_AND_KWARGS = [
            ("activity", {}, 200),
            ("classification_dashboard", {"lab_id": self.lab.pk}, 200),
            ('classification_import_tool', {}, 200),
            ('hgvs_issues', {}, 200),
            # ("clinical_context", {"pk": self.allele.pk}, 200), # Needs clinical context on allele
        ]
        self._test_urls(ADMIN_URL_NAMES_AND_KWARGS, self.user_admin)

    def testUrls(self):
        URL_NAMES_AND_KWARGS = [
            ("classifications", {}, 200),
            ("export_classifications_grid", {}, 200),
            ("export_classifications_grid_redcap", {}, 200),
            ("redcap_data_dictionary", {}, 200),
            ("evidence_keys", {}, 200),
            ('classification_export', {}, 200),
            ("classification_graphs", {}, 200),
            ("summary_email_html", {}, 200),
            ("summary_email_text", {}, 200),
            ("overlaps", {}, 200),
            ("evidence_keys_api", {}, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user_non_owner)

    def testAutocompleteUrls(self):
        """ Autocompletes w/o permissions """
        evidence_key = EvidenceKey.objects.first()

        AUTOCOMPLETE_URLS = [
            ('evidence_key_autocomplete', evidence_key, {"q": evidence_key.key}),
        ]
        self._test_autocomplete_urls(AUTOCOMPLETE_URLS, self.user_non_owner, True)

    def testPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_owner)

    @prevent_request_warnings
    def testNoPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_non_owner, expected_code_override=403)


# TODO: Still need to test below
# ('api/classifications/auto_populate', AutopopulateView.as_view(), name='classification_auto_populate_api'),
# ('api/classifications/record/', ClassificationView.as_view(), name='classification_api'),
# ('api/classifications/record/<record_id>', ClassificationView.as_view(), name='classification_with_record_api'),
# ('api/classifications/datatables/', ClassificationModificationDatatableView.as_view(), name='classification_datatables')


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
