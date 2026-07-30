"""
The derivation and matching rules are covered without a database in
annotation/tests/test_gene_annotation_release_matching.py - these cover only what needs real
model rows: that the candidate queryset is scoped to the build and consortium, and that a
match is saved.
"""
from django.test import TestCase
from django.test.utils import override_settings

from annotation.fake_annotation import get_fake_annotation_settings_dict, get_fake_vep_version
from annotation.models import VariantAnnotationVersion
from genes.models import GeneAnnotationImport, GeneAnnotationRelease
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild

RS_2025_08_URL = ("https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/"
                  "GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz")
REFSEQ_GRCH38_RS_2025_08 = "GCF_000001405.40-RS_2025_08 - GCF_000001405.40_GRCh38.p14_genomic.gff"


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class LinkGeneAnnotationReleaseTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")

    def _make_release(self, genome_build, url, version, annotation_consortium=AnnotationConsortium.REFSEQ):
        gene_annotation_import = GeneAnnotationImport.objects.create(annotation_consortium=annotation_consortium,
                                                                    genome_build=genome_build, url=url)
        return GeneAnnotationRelease.objects.create(version=version, annotation_consortium=annotation_consortium,
                                                   genome_build=genome_build,
                                                   gene_annotation_import=gene_annotation_import)

    def _make_vav(self, genome_build, refseq=REFSEQ_GRCH38_RS_2025_08):
        kwargs = get_fake_vep_version(genome_build, AnnotationConsortium.REFSEQ, 2)
        kwargs["refseq"] = refseq
        return VariantAnnotationVersion.objects.create(**kwargs)

    def test_links_and_saves_matching_release(self):
        release = self._make_release(self.grch38, RS_2025_08_URL, "GRCh38_RefSeq_2025_08")
        vav = self._make_vav(self.grch38)

        self.assertEqual(release, vav.link_gene_annotation_release())
        vav.refresh_from_db()
        self.assertEqual(release, vav.gene_annotation_release)

    def test_candidates_are_scoped_to_our_build(self):
        self._make_release(self.grch37, RS_2025_08_URL, "wrong_build")
        vav = self._make_vav(self.grch38)

        self.assertIsNone(vav.link_gene_annotation_release())
        self.assertIsNone(vav.gene_annotation_release)

    def test_candidates_are_scoped_to_our_consortium(self):
        self._make_release(self.grch38, RS_2025_08_URL, "wrong_consortium",
                           annotation_consortium=AnnotationConsortium.ENSEMBL)
        vav = self._make_vav(self.grch38)

        self.assertIsNone(vav.link_gene_annotation_release())

    def test_ambiguous_candidates_stay_unlinked(self):
        """ Two releases built from the same GFF - a human needs to pick """
        self._make_release(self.grch38, RS_2025_08_URL, "first")
        self._make_release(self.grch38, RS_2025_08_URL.replace("https://", "ftp://"), "second")
        vav = self._make_vav(self.grch38)

        self.assertIsNone(vav.link_gene_annotation_release())
        self.assertIsNone(vav.gene_annotation_release)
