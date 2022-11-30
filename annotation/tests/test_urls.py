import unittest

from django.contrib.auth.models import User

from annotation.fake_annotation import get_fake_annotation_version, create_fake_clinvar_data, \
    create_fake_variant_annotation
from annotation.models import HumanProteinAtlasTissueSample, ClinVar, Citation, CitationSource
from library.django_utils.unittest_utils import URLTestCase
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild


class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        owner_username = f"test_user_{__file__}_owner"
        admin_username = f"test_user_{__file__}_admin"
        cls.user = User.objects.get_or_create(username=owner_username)[0]
        cls.admin_user = User.objects.create_superuser(admin_username)
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")

        cls.annotation_version_grch37 = get_fake_annotation_version(cls.grch37)
        cls.annotation_version_grch38 = get_fake_annotation_version(cls.grch38)

        cls.human_protein_atlas_tissue_sample = HumanProteinAtlasTissueSample.objects.get_or_create(name="foo")[0]

        create_fake_clinvar_data(cls.annotation_version_grch37.clinvar_version)
        clinvar = ClinVar.objects.filter(version=cls.annotation_version_grch37.clinvar_version).first()
        q = Variant.get_contigs_q(cls.grch37) & Variant.get_no_reference_q()
        variant = Variant.objects.filter(q).first()
        create_fake_variant_annotation(variant, cls.annotation_version_grch37.variant_annotation_version)
        citation = Citation.objects.filter(citation_source=CitationSource.PUBMED).first()
        pubmed_citation = f"PubMed:{citation.citation_id}"

        cls.clinvar_id = clinvar.pk
        cls.variant_string = str(variant)
        cls.pubmed_citations = "&".join((str(c) for c in Citation.objects.all().values_list("citation_id", flat=True)[:2]))
        cls.citations_ids_list = "/".join((str(c) for c in Citation.objects.all().values_list("pk", flat=True)[:2]))
        cls.citations_ids_list_pubmed = pubmed_citation

    def testUrls(self):
        """ No permissions to test """
        URL_NAMES_AND_KWARGS = [
            # Disable annotation as it kicks off cached web tasks
            ("annotation", {}, 200),
            # Don't run annotation versions as it kicks off new versions after loading VEP
            # ("annotation_versions", {}, 200),
            # ("version_diffs", {}, 200),
            ("variant_annotation_runs", {}, 403),
            ("view_annotation_descriptions", {}, 200),
            ("about_new_vep_columns", {}, 200),
            ("view_annotation_version_details", {"annotation_version_id": self.annotation_version_grch37.pk}, 200),
            ("clinvar_citations_tab", {"clinvar_id": self.clinvar_id}, 200),
            ("pubmed_citations_tab", {"pubmed_citations": self.pubmed_citations}, 200),
            ("citations_tab", {"citations_ids_list": self.citations_ids_list}, 200),
            ("citations_json", {"citations_ids_list": self.citations_ids_list_pubmed}, 200),
            # API
            ("api_variant_annotation", {"genome_build_name": self.grch37.name, "variant_string": self.variant_string}, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user)

    def testAdminUrls(self):
        URL_NAMES_AND_KWARGS = [
            ("variant_annotation_runs", {}, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.admin_user)

    def testGridUrls(self):
        build_kwargs = {"genome_build_name": self.grch37.name, "op": "config"}
        GRID_LIST_URLS = [
            ("variant_annotation_version_grid", build_kwargs, 200),
        ]
        self._test_urls(GRID_LIST_URLS, self.user)


if __name__ == "__main__":
    unittest.main()
