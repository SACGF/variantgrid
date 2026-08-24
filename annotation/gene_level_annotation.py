"""
Annotation for gene-level variants (gene fusions), computed locally.

@see snpdb.gene_level_variants for why a fusion is stored as a Variant.

These never reach VEP - their alt is not something VEP can parse, and there is no coordinate to
annotate anyway. What they do need is the same thing VEP output gives everything else: a
VariantAnnotation row so grids and hover cards have a symbol to show, VariantTranscriptAnnotation
rows so an analysis that knows its enrichment kit can swap in that kit's transcript, and
VariantGeneOverlap rows so gene lists and compound-het two-hit detection find them.

It is a pipeline type rather than a one-off at import because VariantGeneOverlap is keyed on the
annotation version and symbol-to-gene resolution is per GeneAnnotationRelease. Written once at
import, fusions would quietly drop out of gene lists at the next annotation version.

Resolution runs HGNC/symbol -> release genes rather than HGNC -> gene: HGNC carries no Entrez ID
while Gene.identifier is the Entrez ID for RefSeq releases, so the symbol is the only bridge. That
also means a later HGNC import improving SEPT14 -> SEPTIN14 improves the mapping without touching a
single Variant, since identity keys on the id and only resolution is versioned.

A fusion sits on no transcript, so there is no VEP 'pick' to inherit. VariantAnnotation is per
(version, variant) and has no sample, so the enrichment kit that would name the lab's transcript is
out of reach - the representative is MANE Select, then RefSeq Select, then the highest version we
hold. The kit's own choice is served by the per-transcript rows, which the grid export swaps in
against a CanonicalTranscriptCollection (@see analysis.grid_export).
"""
import logging
import os
from collections import Counter
from dataclasses import dataclass
from typing import Optional

from django.db import transaction
from django.utils import timezone

from annotation.annotation_version_querysets import get_variants_qs_for_annotation
from annotation.models.models import (
    VariantAnnotation,
    VariantAnnotationVersion,
    VariantGeneOverlap,
    VariantTranscriptAnnotation,
)
from annotation.signals.manual_signals import annotation_run_complete_signal
from genes.models import FusionGeneId, GeneAnnotationRelease, GeneFusion, TranscriptVersion
from library.django_utils.django_partition import temporary_db_table

BULK_INSERT_BATCH_SIZE = 2000


@dataclass(frozen=True)
class ReleaseGeneAnnotation:
    """ What one side of a fusion means in a GeneAnnotationRelease """
    symbol: str
    gene_ids: list[str]  # every gene the symbol maps to - one VariantGeneOverlap each
    gene_id: Optional[str]  # the representative one
    transcript_versions: list[TranscriptVersion]  # one VariantTranscriptAnnotation each
    representative_transcript_version: Optional[TranscriptVersion]


class FusionGeneIdResolver:
    """ FusionGeneId -> the genes and transcripts of a GeneAnnotationRelease. One per run, since it
        caches per-symbol lookups across what is usually a lot of repeated partners. """

    def __init__(self, gene_annotation_release: GeneAnnotationRelease):
        self.gene_annotation_release = gene_annotation_release
        self._cache: dict[int, ReleaseGeneAnnotation] = {}

    def get_release_gene_annotation(self, fusion_gene_id: FusionGeneId) -> ReleaseGeneAnnotation:
        if (cached := self._cache.get(fusion_gene_id.pk)) is not None:
            return cached

        gene_ids = []
        transcript_versions = []
        gene_symbol_id = fusion_gene_id.gene_symbol_id
        if self.gene_annotation_release and gene_symbol_id:
            gene_qs = self.gene_annotation_release.genes_for_symbol(gene_symbol_id)
            gene_ids = sorted(gene_qs.values_list("identifier", flat=True))
            tv_qs = self.gene_annotation_release.transcript_versions_for_symbol(gene_symbol_id)
            transcript_versions = list(tv_qs.select_related("transcript", "gene_version"))

        gene_id = gene_ids[0] if gene_ids else None
        result = ReleaseGeneAnnotation(
            symbol=fusion_gene_id.symbol_str,
            gene_ids=gene_ids,
            gene_id=gene_id,
            transcript_versions=transcript_versions,
            representative_transcript_version=_representative_transcript_version(transcript_versions, gene_id))
        self._cache[fusion_gene_id.pk] = result
        return result


def _representative_transcript_version(transcript_versions: list[TranscriptVersion],
                                       gene_id: Optional[str]) -> Optional[TranscriptVersion]:
    """ MANE Select, then RefSeq Select, then the highest version we hold. canonical_score reads
        cdot's transcript tags, so one ranking serves both annotation consortiums.

        Preference goes to the gene VariantAnnotation.gene will name, so the representative row's
        gene and transcript agree where a symbol maps to more than one gene. """

    def _rank(tv: TranscriptVersion) -> tuple:
        return tv.gene_version.gene_id == gene_id, tv.canonical_score, tv.version

    if transcript_versions:
        return max(transcript_versions, key=_rank)
    return None


def annotate_gene_level_run(annotation_run) -> int:
    """ Peer of dump_and_annotate_variants for the GENE_LEVEL pipeline: write the VariantAnnotation +
        VariantTranscriptAnnotation + VariantGeneOverlap rows for every gene-level variant in this
        run's range, and take the run to FINISHED. Returns how many variants were annotated.

        Unlike the VEP pipelines there is no dump, no subprocess and no separate import lane - the
        rows are written here - so this is also where the run completes. AnnotationRun.get_status()
        reads the timestamp fields, so they are walked through the same states a VEP run passes. """

    variant_annotation_version = annotation_run.variant_annotation_version
    resolver = FusionGeneIdResolver(variant_annotation_version.gene_annotation_release)
    if variant_annotation_version.gene_annotation_release is None:
        logging.warning("%s has no gene_annotation_release - gene-level variants get symbols but no "
                        "gene overlaps, so gene lists will not find them", variant_annotation_version)

    start = timezone.now()
    annotation_run.dump_start = start
    annotation_run.annotation_start = start
    annotation_run.upload_start = start
    annotation_run.save()

    results = Counter()
    variant_annotations = []
    transcript_annotations = []
    gene_overlaps = []
    for gene_fusion in _gene_fusions_for_run(annotation_run).iterator():
        variant_annotation, fusion_transcript_annotations, fusion_gene_overlaps = \
            _build_gene_fusion_annotation(annotation_run, resolver, gene_fusion)
        variant_annotations.append(variant_annotation)
        transcript_annotations.extend(fusion_transcript_annotations)
        gene_overlaps.extend(fusion_gene_overlaps)
        results["annotated"] += 1

    with transaction.atomic():
        _bulk_create_in_partition(variant_annotation_version, VariantAnnotation,
                                  VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION,
                                  variant_annotations)
        _bulk_create_in_partition(variant_annotation_version, VariantTranscriptAnnotation,
                                  VariantAnnotationVersion.TRANSCRIPT_ANNOTATION,
                                  transcript_annotations)
        _bulk_create_in_partition(variant_annotation_version, VariantGeneOverlap,
                                  VariantAnnotationVersion.VARIANT_GENE_OVERLAP,
                                  gene_overlaps)
    results["transcript_annotations"] = len(transcript_annotations)
    results["gene_overlaps"] = len(gene_overlaps)

    annotated_count = results["annotated"]
    end = timezone.now()
    annotation_run.dump_count = annotated_count  # 0 here is what makes get_status() finish an empty run
    annotation_run.annotated_count = annotated_count
    annotation_run.dump_end = end
    annotation_run.annotation_end = end
    annotation_run.upload_end = end
    annotation_run.save()

    annotation_run_complete_signal.send(sender=os.path.basename(__file__),
                                        variant_annotation_version=variant_annotation_version,
                                        pipeline_type=annotation_run.pipeline_type,
                                        annotation_run=annotation_run)

    logging.info("%s gene-level annotation: %s", annotation_run, dict(results))
    return annotated_count


def _bulk_create_in_partition(variant_annotation_version: VariantAnnotationVersion,
                              klass, base_table_name: str, records: list):
    """ Annotation tables are partitioned by inheritance and every annotation queryset is rewritten to
        read the version's partition directly, so a row left in the base table is invisible to all of
        them. @see library.django_utils.django_partition and annotation.annotation_version_querysets """

    if not records:
        return
    partition_table = variant_annotation_version.get_partition_table(base_table_name=base_table_name)
    with temporary_db_table(klass, partition_table):
        klass.objects.bulk_create(records, batch_size=BULK_INSERT_BATCH_SIZE)


def _gene_fusions_for_run(annotation_run):
    annotation_version = annotation_run.annotation_range_lock.version.get_any_annotation_version()
    range_lock = annotation_run.annotation_range_lock
    variant_qs = get_variants_qs_for_annotation(annotation_version,
                                                pipeline_type=annotation_run.pipeline_type,
                                                min_variant_id=range_lock.min_variant_id,
                                                max_variant_id=range_lock.max_variant_id)
    return GeneFusion.objects.filter(variant__in=variant_qs) \
        .select_related("variant", "anchor", "partner")


def _build_gene_fusion_annotation(annotation_run, resolver: FusionGeneIdResolver, gene_fusion: GeneFusion):
    """ (representative annotation, per-transcript annotations, gene overlaps) for one fusion """

    variant_annotation_version = annotation_run.variant_annotation_version
    # The anchor's release annotation is the representative one, matching "the gene this row is about"
    # everywhere else
    release_annotations = [resolver.get_release_gene_annotation(fgi) for fgi in gene_fusion.fusion_gene_ids]
    anchor = release_annotations[0]

    # hgvs_c/hgvs_g both carry the VICC gene-level nomenclature - HGVS defers to VICC for fusions,
    # and a blank g.HGVS reads as broken rather than as "not applicable"
    canonical_str = gene_fusion.canonical_str
    representative_transcript_version = anchor.representative_transcript_version
    variant_annotation = VariantAnnotation(
        version=variant_annotation_version,
        variant=gene_fusion.variant,
        annotation_run=annotation_run,
        symbol=anchor.symbol,
        overlapping_symbols=",".join(sorted({ra.symbol for ra in release_annotations})),
        gene_id=anchor.gene_id,
        transcript_id=representative_transcript_version.transcript_id if representative_transcript_version else None,
        transcript_version=representative_transcript_version,
        hgvs_c=canonical_str,
        hgvs_g=canonical_str,
    )

    # Both partners, so an enrichment kit's canonical transcript on either side has a row to swap in
    transcript_annotations = []
    gene_ids = set()
    seen_transcript_version_ids = set()
    for release_annotation in release_annotations:
        gene_ids.update(release_annotation.gene_ids)
        for transcript_version in release_annotation.transcript_versions:
            if transcript_version.pk in seen_transcript_version_ids:
                continue
            seen_transcript_version_ids.add(transcript_version.pk)
            transcript_annotations.append(VariantTranscriptAnnotation(
                version=variant_annotation_version,
                variant=gene_fusion.variant,
                annotation_run=annotation_run,
                symbol=release_annotation.symbol,
                gene_id=transcript_version.gene_version.gene_id,
                transcript_id=transcript_version.transcript_id,
                transcript_version=transcript_version,
                hgvs_c=canonical_str,
            ))

    gene_overlaps = [
        VariantGeneOverlap(version=variant_annotation_version,
                           annotation_run=annotation_run,
                           variant=gene_fusion.variant,
                           gene_id=gene_id)
        for gene_id in sorted(gene_ids)
    ]
    return variant_annotation, transcript_annotations, gene_overlaps
