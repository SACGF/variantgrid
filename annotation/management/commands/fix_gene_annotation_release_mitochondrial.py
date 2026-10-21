"""
Add the 13 mitochondrial coding genes to RefSeq GeneAnnotationReleases built before cdot had them.

A release's ReleaseGeneVersion / ReleaseTranscriptVersion rows decide which genes the genes grid lists, which
genes a gene list symbol matches in analyses and which genes a GeneAnnotationVersion covers. cdot's 'fake-rna-<gene>'
MT transcripts (#2139) arrived after existing releases were imported, so none held an MT gene: a gene list naming
MT-ND4 matched nothing and the MT contig page listed no genes.

RefSeq's NC_012920.1 annotation is identical in every GFF release (GRCh37 105.*, GRCh38 110, RS_2023_10 .. RS_2025_08),
so the build's current fake transcripts fit every release of that build. For each RefSeq release this links them,
matches the release's unmatched HGNC / gene list symbols (the ND4 -> MT-ND4 alias does the work) and writes
GeneAnnotation rows for the new genes into each of the release's existing GeneAnnotationVersions.
Needs the fake transcripts imported first (import_cdot_latest).
"""
from django.core.management.base import BaseCommand

from annotation.management.commands.gene_annotation import Command as GeneAnnotationCommand
from annotation.models import GeneAnnotation, GeneAnnotationVersion
from genes.gene_matching import ReleaseGeneMatcher
from genes.models import (
    GeneAnnotationRelease,
    ReleaseGeneSymbol,
    ReleaseGeneVersion,
    ReleaseTranscriptVersion,
    TranscriptVersion,
)
from genes.models_enums import AnnotationConsortium
from genes.transcript_parts import CDOT_FAKE_TRANSCRIPT_PREFIX


class Command(BaseCommand):
    category = "one-off"

    def handle(self, *args, **options):
        for release in GeneAnnotationRelease.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ):
            add_mitochondrial_genes_to_release(release)


def _fake_transcript_versions(release: GeneAnnotationRelease):
    return TranscriptVersion.objects.filter(genome_build=release.genome_build,
                                            transcript__annotation_consortium=AnnotationConsortium.REFSEQ,
                                            transcript__identifier__startswith=CDOT_FAKE_TRANSCRIPT_PREFIX)


def _match_release_symbols(release: GeneAnnotationRelease):
    """ A symbol that was looked up before the release had its gene (MT-ND4 from a gene list, or HGNC) left a
        ReleaseGeneSymbol with no genes, and match_unmatched_in_hgnc_and_gene_lists skips any symbol with a row -
        so those are matched again explicitly, then symbols never looked up """
    matcher = ReleaseGeneMatcher(release)
    unmatched = list(ReleaseGeneSymbol.objects.filter(release=release, releasegenesymbolgene__isnull=True))
    if unmatched:
        matcher.match_symbols_to_genes(unmatched)
    matcher.match_unmatched_in_hgnc_and_gene_lists()


def add_mitochondrial_genes_to_release(release: GeneAnnotationRelease) -> set[str]:
    """ Returns the gene ids the release gained """
    fake_tv_qs = _fake_transcript_versions(release)
    if not fake_tv_qs.exists():
        print(f"{release}: no cdot fake-rna transcripts for {release.genome_build} - run import_cdot_latest first")
        return set()

    existing_tv_ids = set(release.releasetranscriptversion_set.values_list("transcript_version_id", flat=True))
    existing_gv_ids = set(release.releasegeneversion_set.values_list("gene_version_id", flat=True))
    new_gene_ids = set()
    release_transcript_versions = []
    release_gene_versions = []
    for tv_id, gv_id, gene_id in fake_tv_qs.values_list("pk", "gene_version_id", "gene_version__gene_id"):
        if tv_id not in existing_tv_ids:
            release_transcript_versions.append(ReleaseTranscriptVersion(release=release, transcript_version_id=tv_id))
        if gv_id not in existing_gv_ids:
            existing_gv_ids.add(gv_id)
            release_gene_versions.append(ReleaseGeneVersion(release=release, gene_version_id=gv_id))
            new_gene_ids.add(gene_id)

    if not (release_transcript_versions or release_gene_versions):
        print(f"{release}: already has the mitochondrial genes")
        return set()

    ReleaseTranscriptVersion.objects.bulk_create(release_transcript_versions)
    ReleaseGeneVersion.objects.bulk_create(release_gene_versions)
    print(f"{release}: added {len(release_gene_versions)} gene versions, "
          f"{len(release_transcript_versions)} transcript versions")

    _match_release_symbols(release)

    gene_annotation_command = GeneAnnotationCommand()
    for gav in GeneAnnotationVersion.objects.filter(gene_annotation_release=release, data_archived_date__isnull=True):
        annotated = set(GeneAnnotation.objects.filter(version=gav, gene_id__in=new_gene_ids).values_list("gene_id",
                                                                                                         flat=True))
        if gene_ids := new_gene_ids - annotated:
            num_written = gene_annotation_command.add_gene_annotation_for_genes(gav, gene_ids)
            print(f"{gav}: wrote {num_written} GeneAnnotation rows")
    return new_gene_ids
