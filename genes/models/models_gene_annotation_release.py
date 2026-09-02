
from django.db import models
from django.db.models import QuerySet
from django.db.models.deletion import CASCADE
from django_extensions.db.models import TimeStampedModel

from genes.models.models_gene import (
    Gene,
    GeneAnnotationImport,
    GeneSymbol,
    GeneVersion,
    TranscriptVersion,
)
from genes.models_enums import AnnotationConsortium
from snpdb.models.models_genome import GenomeBuild


class GeneAnnotationRelease(models.Model):
    """
        Links all the GeneVersion/TranscriptVersion from a particular GTF/GFF file, so they match what a
        VEP build used. Created by the 'import_cdot_gene_annotation_release' management command.

        This release can be set on a VariantAnnotationVersion to be able to get genes/transcripts from a VEP build
    """
    version = models.TextField()  # Needs to support e.g. "109.20190607"
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    gene_annotation_import = models.ForeignKey(GeneAnnotationImport, on_delete=CASCADE)

    class Meta:
        unique_together = ('version', 'annotation_consortium', 'genome_build')

    def get_genes(self):
        return Gene.objects.filter(annotation_consortium=self.annotation_consortium,
                                   geneversion__releasegeneversion__release=self)

    @staticmethod
    def get_for_latest_annotation_versions_for_builds(status=None) -> list['GeneAnnotationRelease']:
        """ status: VariantAnnotationVersion.Status - defaults to ACTIVE (see VariantAnnotationVersion.latest) """
        from annotation.models import VariantAnnotationVersion
        gene_annotation_releases = []
        for genome_build in GenomeBuild.builds_with_annotation().order_by("name"):
            if vav := VariantAnnotationVersion.latest(genome_build, status=status):
                if vav.gene_annotation_release:
                    gene_annotation_releases.append(vav.gene_annotation_release)
        return gene_annotation_releases

    def genes_for_symbols(self, gene_symbols) -> QuerySet:
        rgsg_qs = ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=self,
                                                       release_gene_symbol__gene_symbol__in=gene_symbols)
        return Gene.objects.filter(pk__in=rgsg_qs.values_list("gene_id", flat=True))

    def genes_for_symbol(self, gene_symbol) -> QuerySet:
        return self.genes_for_symbols([gene_symbol])

    def transcript_versions_for_transcript(self, transcript) -> QuerySet[TranscriptVersion]:
        return TranscriptVersion.objects.filter(releasetranscriptversion__release=self,
                                                transcript=transcript)

    def transcript_versions_for_gene(self, gene) -> QuerySet:
        return TranscriptVersion.objects.filter(releasetranscriptversion__release=self,
                                                gene_version__gene=gene)

    def transcript_versions_for_symbol(self, gene_symbol) -> QuerySet:
        return TranscriptVersion.objects.filter(releasetranscriptversion__release=self,
                                                gene_version__gene__in=self.genes_for_symbol(gene_symbol))

    def __str__(self):
        return f"{self.genome_build.slug}/{self.get_annotation_consortium_display()} - v{self.version}"


class ReleaseGeneVersion(models.Model):
    """ Used to be able to reconstruct what versions were in an import """
    release = models.ForeignKey(GeneAnnotationRelease, on_delete=CASCADE)
    gene_version = models.ForeignKey(GeneVersion, on_delete=CASCADE)

    class Meta:
        unique_together = ('release', 'gene_version')


class ReleaseTranscriptVersion(models.Model):
    """ Used to be able to reconstruct what versions were in an import """
    release = models.ForeignKey(GeneAnnotationRelease, on_delete=CASCADE)
    transcript_version = models.ForeignKey(TranscriptVersion, on_delete=CASCADE)

    class Meta:
        unique_together = ('release', 'transcript_version')


class ReleaseGeneSymbol(TimeStampedModel):
    """ Collects (via GeneSymbolReleaseGene) how gene symbol matches genes for a particular release
        Works as global cache, created as symbols appear in gene lists """
    release = models.ForeignKey(GeneAnnotationRelease, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)

    class Meta:
        unique_together = ('release', 'gene_symbol')

    def __str__(self):
        return f"{self.release}: {self.gene_symbol}"


class ReleaseGeneSymbolGene(TimeStampedModel):
    release_gene_symbol = models.ForeignKey(ReleaseGeneSymbol, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    match_info = models.TextField(null=True)

    class Meta:
        unique_together = ('release_gene_symbol', 'gene')
