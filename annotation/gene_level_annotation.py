"""
Annotation for gene-level variants (gene fusions), computed locally.

@see snpdb.gene_level_variants for why a fusion is stored as a Variant.

These never reach VEP - their alt is not something VEP can parse, and there is no coordinate to
annotate anyway. What they do need is the same thing VEP output gives everything else: a
VariantAnnotation row so grids and hover cards have a symbol to show, and VariantGeneOverlap rows so
gene lists and compound-het two-hit detection find them.

It is a pipeline type rather than a one-off at import because VariantGeneOverlap is keyed on the
annotation version and symbol-to-gene resolution is per GeneAnnotationRelease. Written once at
import, fusions would quietly drop out of gene lists at the next annotation version.

Resolution runs HGNC/symbol -> release genes rather than HGNC -> gene: HGNC carries no Entrez ID
while Gene.identifier is the Entrez ID for RefSeq releases, so the symbol is the only bridge. That
also means a later HGNC import improving SEPT14 -> SEPTIN14 improves the mapping without touching a
single Variant, since identity keys on the id and only resolution is versioned.
"""
import logging
import os
from collections import Counter

from django.db import transaction
from django.utils import timezone

from annotation.annotation_version_querysets import get_variants_qs_for_annotation
from annotation.models.models import VariantAnnotation, VariantGeneOverlap
from annotation.signals.manual_signals import annotation_run_complete_signal
from genes.models import GeneAnnotationRelease, GeneFusion, FusionGeneId


class FusionGeneIdResolver:
    """ FusionGeneId -> the genes and symbols of a GeneAnnotationRelease. One per run, since it
        caches per-symbol lookups across what is usually a lot of repeated partners. """

    def __init__(self, gene_annotation_release: GeneAnnotationRelease):
        self.gene_annotation_release = gene_annotation_release
        self._cache: dict[int, tuple[str, list[str]]] = {}

    def get_symbol_and_gene_ids(self, fusion_gene_id: FusionGeneId) -> tuple[str, list[str]]:
        """ (symbol as it should be displayed, gene ids in this release) """
        if (cached := self._cache.get(fusion_gene_id.pk)) is not None:
            return cached

        gene_ids = []
        if self.gene_annotation_release and fusion_gene_id.gene_symbol_id:
            gene_qs = self.gene_annotation_release.genes_for_symbol(fusion_gene_id.gene_symbol_id)
            gene_ids = list(gene_qs.values_list("identifier", flat=True))

        result = (fusion_gene_id.symbol_str, gene_ids)
        self._cache[fusion_gene_id.pk] = result
        return result


def annotate_gene_level_run(annotation_run) -> int:
    """ Peer of dump_and_annotate_variants for the GENE_LEVEL pipeline: write the VariantAnnotation +
        VariantGeneOverlap rows for every gene-level variant in this run's range, and take the run to
        FINISHED. Returns how many variants were annotated.

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
    for gene_fusion in _gene_fusions_for_run(annotation_run).iterator():
        with transaction.atomic():
            _annotate_gene_fusion(annotation_run, resolver, gene_fusion)
        results["annotated"] += 1

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


def _gene_fusions_for_run(annotation_run):
    annotation_version = annotation_run.annotation_range_lock.version.get_any_annotation_version()
    range_lock = annotation_run.annotation_range_lock
    variant_qs = get_variants_qs_for_annotation(annotation_version,
                                                pipeline_type=annotation_run.pipeline_type,
                                                min_variant_id=range_lock.min_variant_id,
                                                max_variant_id=range_lock.max_variant_id)
    return GeneFusion.objects.filter(variant__in=variant_qs) \
        .select_related("variant", "anchor", "partner")


def _annotate_gene_fusion(annotation_run, resolver: FusionGeneIdResolver, gene_fusion: GeneFusion):
    variant_annotation_version = annotation_run.variant_annotation_version

    symbols = []
    gene_ids = set()
    for fusion_gene_id in gene_fusion.fusion_gene_ids:
        symbol, resolved_gene_ids = resolver.get_symbol_and_gene_ids(fusion_gene_id)
        symbols.append(symbol)
        gene_ids.update(resolved_gene_ids)

    # The anchor's symbol is the representative one, matching "the gene this row is about" everywhere else
    VariantAnnotation.objects.update_or_create(
        version=variant_annotation_version,
        variant=gene_fusion.variant,
        defaults={
            "annotation_run": annotation_run,
            "symbol": symbols[0],
            "overlapping_symbols": ",".join(sorted(set(symbols))),
        })

    overlaps = [
        VariantGeneOverlap(version=variant_annotation_version,
                           annotation_run=annotation_run,
                           variant=gene_fusion.variant,
                           gene_id=gene_id)
        for gene_id in sorted(gene_ids)
    ]
    if overlaps:
        VariantGeneOverlap.objects.bulk_create(overlaps, ignore_conflicts=True)
