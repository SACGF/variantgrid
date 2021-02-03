from django.contrib.auth.models import User
import unittest

from annotation.fake_annotation import get_fake_annotation_version, create_fake_clinvar_data, \
    create_fake_variant_annotation
from annotation.models import HumanProteinAtlasAbundance, HumanProteinAtlasTissueSample, \
    ClinVar, Citation, CitationSource
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from library.django_utils.unittest_utils import URLTestCase
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild


class Test(URLTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls.user = User.objects.get_or_create(username='testuser')[0]
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
        transcript_version = create_fake_transcript_version(cls.grch37)

        cls.gene_id = transcript_version.gene_version.gene_id
        cls.gene_symbol = transcript_version.gene_version.gene_symbol

    def testUrls(self):
        """ No permissions to test """
        URL_NAMES_AND_KWARGS = [
            # Disable annotation as it kicks off cached web tasks
            ("annotation", {}, 200),
            # Don't run annotation versions as it kicks off new versions after loading VEP
            # ("annotation_versions", {}, 200),
            # ("version_diffs", {}, 200),
            ("variant_annotation_runs", {}, 200),
            ("view_annotation_descriptions", {}, 200),
            ("about_new_vep_columns", {}, 200),
            ("view_annotation_version_details", {"annotation_version_id": self.annotation_version_grch37.pk}, 200),
            ("clinvar_citations_tab", {"clinvar_id": self.clinvar_id}, 200),
            ("pubmed_citations_tab", {"pubmed_citations": self.pubmed_citations}, 200),
            ("citations_tab", {"citations_ids_list": self.citations_ids_list}, 200),
            ("citations_json", {"citations_ids_list": self.citations_ids_list_pubmed}, 200),
            # API
            ("api_view_gene_disease_validity", {"gene_symbol": self.gene_symbol}, 200),
            ("api_gene_annotation", {"gene_symbol": self.gene_symbol}, 200),
            ("api_variant_annotation", {"genome_build_name": self.grch37.name, "variant_string": self.variant_string}, 200),
        ]
        self._test_urls(URL_NAMES_AND_KWARGS, self.user)

    def testGridUrls(self):
        build_kwargs = {"genome_build_name": self.grch37.name, "op": "config"}
        human_protein_atlas_version_id = self.annotation_version_grch37.human_protein_atlas_version.pk
        tissue_kwargs = {"human_protein_atlas_version_id": human_protein_atlas_version_id,
                         "tissue_sample_id": self.human_protein_atlas_tissue_sample.pk,
                         "min_abundance": HumanProteinAtlasAbundance.MEDIUM,
                         "op": "config"}
        GRID_LIST_URLS = [
            ("variant_annotation_version_grid", build_kwargs, 200),
            ("annotation_run_grid", build_kwargs, 200),
            ("tissue_gene_grid", tissue_kwargs, 200),
        ]
        self._test_urls(GRID_LIST_URLS, self.user)


if __name__ == "__main__":
    unittest.main()
