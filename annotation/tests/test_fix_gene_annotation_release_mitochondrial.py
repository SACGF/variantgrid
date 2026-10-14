"""
fix_gene_annotation_release_mitochondrial adds cdot's fake MT transcripts to a RefSeq release imported before cdot had
them, rematches a symbol looked up while the release lacked the gene, and writes GeneAnnotation rows for the new
genes into the release's existing GeneAnnotationVersion.
"""
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.fake_data import get_fake_annotation_settings_dict, get_fake_annotation_version
from annotation.management.commands.fix_gene_annotation_release_mitochondrial import (
    add_mitochondrial_genes_to_release,
)
from annotation.models import (
    DBNSFPGeneAnnotation,
    DBNSFPGeneAnnotationVersion,
    GeneAnnotation,
    GeneAnnotationVersion,
)
from annotation.models.models import SubVersionPartition
from genes.fake_data import create_fake_gene_version
from genes.models import (
    GeneAnnotationImport,
    GeneAnnotationRelease,
    GeneSymbol,
    GeneSymbolAlias,
    ReleaseGeneSymbol,
    ReleaseGeneSymbolGene,
    ReleaseGeneVersion,
    ReleaseTranscriptVersion,
    Transcript,
    TranscriptVersion,
)
from genes.models_enums import AnnotationConsortium, GeneSymbolAliasSource
from snpdb.models import GenomeBuild


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class FixGeneAnnotationReleaseMitochondrialTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        annotation_version = get_fake_annotation_version(cls.grch38)
        cls.release = cls._make_refseq_release(cls.grch38)
        dbnsfp_version = DBNSFPGeneAnnotationVersion.objects.create(version="4.test", sha256_hash="not_a_real_hash")
        with SubVersionPartition.defer_new_sub_version():
            cls.gav = GeneAnnotationVersion.objects.create(gene_annotation_release=cls.release,
                                                           ontology_version=annotation_version.ontology_version,
                                                           dbnsfp_gene_version=dbnsfp_version,
                                                           gnomad_import_date=timezone.now())

        # ND4 and its cdot fake transcript, as import_cdot_latest stores them
        gene_version = create_fake_gene_version(cls.grch38, "4538", "ND4", AnnotationConsortium.REFSEQ)
        transcript = Transcript.objects.create(identifier="fake-rna-ND4",
                                               annotation_consortium=AnnotationConsortium.REFSEQ)
        cls.transcript_version = TranscriptVersion.objects.create(transcript=transcript, version=1,
                                                                  gene_version=gene_version,
                                                                  genome_build=cls.grch38,
                                                                  contig=cls.grch38.chrom_contig_mappings["MT"],
                                                                  import_source=gene_version.import_source)
        cls.gene = gene_version.gene

        # HGNC's symbol, which gene lists and dbNSFP use, with RefSeq's as its alias
        mt_nd4 = GeneSymbol.objects.get_or_create(symbol="MT-ND4")[0]
        GeneSymbolAlias.objects.create(alias="ND4", gene_symbol=mt_nd4, source=GeneSymbolAliasSource.HGNC)
        DBNSFPGeneAnnotation.objects.create(version=dbnsfp_version, gene_symbol=mt_nd4)
        # Looked up by a gene list while the release had no MT genes, so it matched nothing
        cls.release_gene_symbol = ReleaseGeneSymbol.objects.create(release=cls.release, gene_symbol=mt_nd4)

    @staticmethod
    def _make_refseq_release(genome_build) -> GeneAnnotationRelease:
        gene_annotation_import = GeneAnnotationImport.objects.create(annotation_consortium=AnnotationConsortium.REFSEQ,
                                                                    genome_build=genome_build,
                                                                    url=f"refseq_gff_{genome_build.name}")
        return GeneAnnotationRelease.objects.create(version="RS_test", annotation_consortium=AnnotationConsortium.REFSEQ,
                                                   genome_build=genome_build,
                                                   gene_annotation_import=gene_annotation_import)

    def test_adds_release_rows_symbol_match_and_gene_annotation(self):
        self.assertEqual({self.gene.pk}, add_mitochondrial_genes_to_release(self.release))

        self.assertTrue(ReleaseGeneVersion.objects.filter(release=self.release,
                                                          gene_version=self.transcript_version.gene_version).exists())
        self.assertTrue(ReleaseTranscriptVersion.objects.filter(release=self.release,
                                                                transcript_version=self.transcript_version).exists())
        self.assertTrue(ReleaseGeneSymbolGene.objects.filter(release_gene_symbol=self.release_gene_symbol,
                                                             gene=self.gene).exists())
        gene_annotation = GeneAnnotation.objects.get(version=self.gav, gene=self.gene)
        self.assertEqual("MT-ND4", gene_annotation.dbnsfp_gene.gene_symbol_id)

    def test_second_run_changes_nothing(self):
        add_mitochondrial_genes_to_release(self.release)

        self.assertEqual(set(), add_mitochondrial_genes_to_release(self.release))
        self.assertEqual(1, GeneAnnotation.objects.filter(version=self.gav, gene=self.gene).count())

    def test_build_without_fake_transcripts_is_left_alone(self):
        release = self._make_refseq_release(GenomeBuild.get_name_or_alias("GRCh37"))

        self.assertEqual(set(), add_mitochondrial_genes_to_release(release))
        self.assertFalse(release.releasegeneversion_set.exists())
