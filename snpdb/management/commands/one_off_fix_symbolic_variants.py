"""
Brings stored del/dup/inv variants into line with settings.VARIANT_SYMBOLIC_ALT_SIZE: explicit ones at or
over it become symbolic, symbolic ones under it become explicit. Where both forms of a variant are stored
(#982, #2109), the copy without genotypes is merged into the other - its allele, classifications, tags and
ClinVar records moved across - and deleted.

Run from migrations as a ManualOperation whenever the threshold changes (#1358). --dry-run reports the counts.
"""
from collections import Counter
from functools import cache

from django.conf import settings
from django.core.management.base import BaseCommand
from django.db import transaction
from django.db.models import ProtectedError, Q
from django.db.models.functions import Length

from annotation.annotation_pipeline_routing import EXPLICIT_SYMBOLIC_ALTS, pipeline_type_for_variant
from annotation.models import AnnotationRangeLock, VariantAnnotation
from library.utils import sha256sum_str
from snpdb.models import AlleleConversionTool, GenomeBuild, Locus, Sequence, Variant, VariantAllele


@cache
def _get_sequence(seq: str) -> Sequence:
    """ By hash - seq itself isn't indexed. Created on first sight, so eg <DEL> won't exist on a legacy DB """
    sequence, _ = Sequence.objects.get_or_create(seq_sha256_hash=sha256sum_str(seq), defaults={"seq": seq})
    return sequence


class _NotMerged(Exception):
    pass


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        for genome_build in GenomeBuild.builds_with_annotation():
            stats = Counter()
            for v in self._non_canonical_candidates(genome_build):
                self._canonicalise(dry_run, genome_build, v, stats)
            self._print_stats(f"{genome_build} variants around VARIANT_SYMBOLIC_ALT_SIZE="
                              f"{settings.VARIANT_SYMBOLIC_ALT_SIZE}", stats)

    @staticmethod
    def _non_canonical_candidates(genome_build):
        """ Explicit variants with a sequence long enough to be symbolic, and symbolic ones short enough to
            be explicit. Most explicit candidates are insertions or substitutions that stay as they are """
        size = settings.VARIANT_SYMBOLIC_ALT_SIZE
        long_sequences = Sequence.objects.annotate(seq_length=Length("seq")).filter(seq_length__gt=size)
        q_explicit = Q(svlen__isnull=True) & (Q(locus__ref__in=long_sequences) | Q(alt__in=long_sequences))
        q_short_symbolic = Q(alt__seq__in=EXPLICIT_SYMBOLIC_ALTS, svlen__gt=-size, svlen__lt=size)
        # select_related as v.coordinate reads locus/contig/ref/alt
        qs = Variant.objects.filter(Variant.get_contigs_q(genome_build), q_explicit | q_short_symbolic)
        return qs.select_related("locus__contig", "locus__ref", "alt").order_by("pk").iterator(chunk_size=1000)

    def _canonicalise(self, dry_run: bool, genome_build, v: Variant, stats: Counter):
        vc = v.coordinate
        canonical = vc.as_internal_canonical_form(genome_build)
        if canonical == vc:
            stats["already canonical - no change"] += 1
            return

        direction = "symbolic" if canonical.is_symbolic else "explicit"
        try:
            twin = Variant.get_from_variant_coordinate(canonical, genome_build)
        except Variant.DoesNotExist:
            stats[f"converted to {direction}"] += 1
            if not dry_run:
                self._convert_in_place(v, canonical)
            return

        stats[f"merged with existing {direction} twin"] += 1
        if not dry_run:
            try:
                self._merge_twins(genome_build, v, twin, canonical, stats)
            except _NotMerged as e:
                stats[f"NOT merged - {e}"] += 1

    @staticmethod
    def _convert_in_place(v: Variant, canonical):
        """ Same pk, so genotypes, annotation and range locks stay where they are. A converted variant under
            VariantAnnotationVersion.structural_variant_min_size keeps its STANDARD annotation, as that is still
            the pipeline it routes to """
        v.locus = Locus.objects.get_or_create(contig=v.locus.contig, position=canonical.position,
                                              ref=_get_sequence(canonical.ref))[0]
        v.alt = _get_sequence(canonical.alt)
        v.svlen = canonical.svlen
        v.end = canonical.end  # Stored calculated field - nothing recalcs it on save
        v.save()

    def _merge_twins(self, genome_build, v: Variant, twin: Variant, canonical, stats: Counter):
        """ Keep whichever has genotypes (the canonical twin if neither), move everything else onto it,
            then delete the other. A kept non-canonical variant is converted once the twin is out of the way
            of the (locus, alt, svlen) unique constraint """
        v_has_genotypes = v.cohortgenotype_set.exists()
        if v_has_genotypes and twin.cohortgenotype_set.exists():
            raise _NotMerged("both have cohort genotype data")
        if v_has_genotypes:
            keep, dupe = v, twin
        else:
            keep, dupe = twin, v

        with transaction.atomic():
            self._move_allele(genome_build, keep, dupe)
            self._move_records(keep, dupe)
            AnnotationRangeLock.release_variant(dupe)
            try:
                dupe.delete()
            except ProtectedError as e:
                raise _NotMerged(f"protected ({e.protected_objects.model.__name__})") from e
            if keep == v:
                self._convert_in_place(keep, canonical)

        annotated_by = VariantAnnotation.objects.filter(variant=keep) \
            .values_list("annotation_run__pipeline_type", "version__structural_variant_min_size")
        if any(pipeline_type != pipeline_type_for_variant(keep, sv_min_size)
               for pipeline_type, sv_min_size in annotated_by):
            # The scheduler only annotates new pk ranges, so this stays as it is until a new annotation version
            stats["kept variant annotated by a pipeline it no longer routes to"] += 1

    @staticmethod
    def _move_allele(genome_build, keep: Variant, dupe: Variant):
        dupe_va = VariantAllele.objects.filter(variant=dupe, genome_build=genome_build).first()
        if dupe_va is None:
            return
        keep_va = VariantAllele.objects.filter(variant=keep, genome_build=genome_build).first()
        if keep_va is None:
            dupe_va.variant = keep
            dupe_va.save()
        elif keep_va.allele != dupe_va.allele:
            if not keep_va.allele.merge(AlleleConversionTool.SAME_CONTIG, dupe_va.allele):
                raise _NotMerged("alleles could not be merged (both have ClinGen alleles)")

    @staticmethod
    def _move_records(keep: Variant, dupe: Variant):
        """ What would block the delete (PROTECT) or be lost with it. Rows keyed on the variant alone -
            annotation, genotype counts, caches - are recalculated or go with it """
        if hasattr(dupe, "variantwiki"):
            if hasattr(keep, "variantwiki"):
                raise _NotMerged("both have a variant wiki")
            dupe.variantwiki.variant = keep
            dupe.variantwiki.save()
        dupe.classification_set.update(variant=keep)
        dupe.varianttag_set.update(variant=keep)
        dupe.clinvar_set.update(variant=keep)
        dupe.importedalleleinfo_set.update(matched_variant=keep)
        dupe.resolvedvariantinfo_set.update(variant=keep)
        dupe.modifiedimportedvariant_set.update(variant=keep)
        dupe.createdmanualvariant_set.update(variant=keep)
        dupe.variantcollectionrecord_set.update(variant=keep)
        dupe.candidate_set.update(variant=keep)

    @staticmethod
    def _print_stats(label: str, stats: Counter):
        print(f"{label}:")
        if not stats:
            print("  (nothing)")
        for reason, count in stats.most_common():
            print(f"  {count}\t{reason}")
