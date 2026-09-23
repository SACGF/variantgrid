"""Bring the gene-level strings stored before #1875's wording change into line with what we write now.

Two things changed at once. A classification target such as 'EGFR amplification' used to have its
whitespace stripped like any c.HGVS (tidy_hgvs_whitespace now keeps a gene-level value's spaces), so
the imported value and the c_hgvs evidence read 'EGFRamplification' while the g.HGVS annotation
wrote 'EGFR amplification'. And a splice event's canonical string gained the word 'splice'
(@see genes.gene_splice.display_splice_label), which annotation and the evidence autopopulated
from it still hold without.

The copy number space is put back deterministically - one gene, one word, so 'FGF4amplification'
can only have been 'FGF4 amplification'. A splice or fusion target is left as the lab wrote it,
since 'EGFRvIII' and 'ARV7' are both real spellings and nothing records which one arrived. The
allele info's md5 moves with the value, and a record already holding the spaced value is reported
rather than merged.

Splice annotation is rewritten in place - hgvs_c / hgvs_g on every version's partition, for the
handful of splice variants there are - and the g_hgvs and splice_label evidence that was
autopopulated from it goes through Classification.patch_history so the modification history agrees.
"""
import logging

from django.core.management.base import BaseCommand
from django.db import transaction
from django.db.models import Q

from annotation.models import (
    AnnotationRun,
    VariantAnnotation,
    VariantAnnotationPipelineType,
    VariantAnnotationVersion,
    VariantTranscriptAnnotation,
)
from classification.enums import SpecialEKeys
from classification.models import Classification, ImportedAlleleInfo
from classification.models.classification_utils import PatchMeta
from genes.gene_copy_number import COPY_NUMBER_STRING_PATTERN
from genes.gene_splice import SPLICE_DISPLAY_SUFFIX, get_splice_event_variant
from library.django_utils.django_partition import temporary_db_table
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from library.utils import md5sum_str
from snpdb.models import Variant

SPLICE_ALT_PREFIX = f"<{GeneLevelSymbolicAlt.SPLICE}:"


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report what would change, and change nothing")

    def handle(self, *args, **options):
        dry_run = options['dry_run']
        respaced = respace_copy_number_allele_infos(dry_run)
        relabelled = relabel_splice_variants(dry_run)
        verb = "Would change" if dry_run else "Changed"
        print(f"{verb} {respaced} copy number allele info(s) and {relabelled} splice variant(s)")


def respaced_copy_number_string(value: str):
    """ 'FGF4amplification' -> 'FGF4 amplification'; None for anything else, including a value
        that already has its space """
    if not value or " " in value:
        return None
    if m := COPY_NUMBER_STRING_PATTERN.match(value):
        gene_name, word = m.groups()
        return f"{gene_name} {word}"
    return None


def respace_copy_number_allele_infos(dry_run: bool) -> int:
    changed = 0
    for allele_info in ImportedAlleleInfo.objects.filter(imported_c_hgvs__isnull=False).order_by("pk"):
        old = allele_info.imported_c_hgvs
        new = respaced_copy_number_string(old)
        if new is None:
            continue
        new_md5 = md5sum_str(new)
        clash = ImportedAlleleInfo.objects.filter(
            imported_md5_hash=new_md5, imported_transcript=allele_info.imported_transcript,
            imported_genome_build_patch_version=allele_info.imported_genome_build_patch_version).first()
        if clash:
            print(f"allele info {allele_info.pk}\t{old}\t-> {new} already held by allele info {clash.pk} - "
                  f"left for you to merge")
            continue
        classifications = list(Classification.objects.filter(allele_info=allele_info))
        print(f"allele info {allele_info.pk}\t{old}\t-> {new}\t{len(classifications)} classification(s)")
        if not dry_run:
            _respace_allele_info(allele_info, new, new_md5, classifications)
        changed += 1
    return changed


@transaction.atomic
def _respace_allele_info(allele_info: ImportedAlleleInfo, new: str, new_md5: str,
                         classifications: list[Classification]):
    old = allele_info.imported_c_hgvs
    allele_info.imported_c_hgvs = new
    allele_info.imported_md5_hash = new_md5
    allele_info.save(update_fields=["imported_c_hgvs", "imported_md5_hash"])
    for classification in classifications:
        _rewrite_evidence(classification, {SpecialEKeys.C_HGVS: (old, new)})


def _rewrite_evidence(classification: Classification, replacements: dict[str, tuple[str, str]]):
    """ Each key's old value becomes the new one wherever the history holds it - the evidence, every
        modification's delta, and the published evidence a modification was shared with """
    def patcher(patch_meta: PatchMeta):
        for key, (old, new) in replacements.items():
            if patch_meta.get(key, fallback_existing=False) == old:
                patch_meta.patch_value(key, new)

    if classification.classificationmodification_set.exists():
        classification.patch_history(patcher)
    else:
        patch_meta = PatchMeta(patch=classification.evidence, existing={})
        patcher(patch_meta)
        classification.save(update_fields=["evidence"])


def relabel_splice_variants(dry_run: bool) -> int:
    """ Every splice Variant's annotation and evidence, from the string display_splice_label wrote
        before the suffix to the one it writes now """
    annotated_version_ids = AnnotationRun.objects.filter(
        pipeline_type=VariantAnnotationPipelineType.GENE_LEVEL) \
        .values_list("annotation_range_lock__version", flat=True).distinct()
    versions = list(VariantAnnotationVersion.objects.filter(pk__in=annotated_version_ids))

    changed = 0
    variant_qs = Variant.objects.filter(Variant.get_gene_level_q(), alt__seq__startswith=SPLICE_ALT_PREFIX)
    for variant in variant_qs.select_related("locus", "alt").order_by("pk"):
        event = get_splice_event_variant(variant)
        if event is None:
            continue
        new = event.canonical_str
        if not new.endswith(SPLICE_DISPLAY_SUFFIX):
            continue
        old = new[:-len(SPLICE_DISPLAY_SUFFIX)]
        classifications = list(Classification.objects.filter(
            Q(allele_info__grch37__variant=variant) | Q(allele_info__grch38__variant=variant)).distinct())
        print(f"variant {variant.pk}\t{old}\t-> {new}\t{len(versions)} annotation version(s), "
              f"{len(classifications)} classification(s)")
        if not dry_run:
            _relabel_splice_variant(variant, old, new, versions, classifications)
        changed += 1
    return changed


@transaction.atomic
def _relabel_splice_variant(variant: Variant, old: str, new: str,
                            versions: list[VariantAnnotationVersion], classifications: list[Classification]):
    for version in versions:
        representative_table = version.get_partition_table(
            base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantAnnotation, representative_table):
            VariantAnnotation.objects.filter(variant=variant, hgvs_c=old).update(hgvs_c=new, hgvs_g=new)
        transcript_table = version.get_partition_table(base_table_name=VariantAnnotationVersion.TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantTranscriptAnnotation, transcript_table):
            VariantTranscriptAnnotation.objects.filter(variant=variant, hgvs_c=old).update(hgvs_c=new)
    # The splice label and g.HGVS were autopopulated from the event when the classification was
    # created off the variant; a lab's own imported c.HGVS is its own and stays
    for classification in classifications:
        _rewrite_evidence(classification, {SpecialEKeys.G_HGVS: (old, new),
                                           SpecialEKeys.SPLICE_LABEL: (old, new)})
    logging.info("Relabelled splice variant %s: %s -> %s", variant.pk, old, new)
