from django.test import TestCase, override_settings
from django.urls import reverse

from classification.enums import SpecialEKeys
from classification.models import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from snpdb.models import GenomeBuild


@override_settings(CLASSIFICATION_WEB_FORM_CREATE=True,
                   CLASSIFICATION_WEB_FORM_CREATE_BY_NON_ADMIN=True,
                   CLASSIFICATION_WEB_FORM_CREATE_ALLOW_NO_VARIANT=True)
class CreateClassificationFromHGVSTest(TestCase):
    """ HGVS that can't be resolved is classified without a variant, via its own create page """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()
        self.client.force_login(self.user)

    def test_create_without_variant(self):
        genome_build = GenomeBuild.grch38()
        hgvs_string = "NM_999999.1:c.100A>G"
        url = reverse("create_classification_from_hgvs",
                      kwargs={"genome_build_name": genome_build.name, "hgvs_string": hgvs_string})
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["lab_form"]["lab"].value(), self.lab.pk)

        response = self.client.post(reverse("create_classification"), {
            "genome_build_name": genome_build.name,
            "evidence_json": response.context["evidence_json"],
            "refseq_transcript_accession": response.context["refseq_transcript_accession"],
            "lab": self.lab.pk,
        })
        self.assertEqual(response.status_code, 302)
        classification = Classification.objects.get(lab=self.lab, user=self.user)
        self.assertIsNone(classification.variant)
        self.assertEqual(classification.get(SpecialEKeys.C_HGVS), hgvs_string)
