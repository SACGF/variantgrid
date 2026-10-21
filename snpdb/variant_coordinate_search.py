""" Matching a VariantCoordinate to existing Variants, falling back to the variants that bcftools normalised it into
    on import (ModifiedImportedVariant). Entry point: get_results_from_variant_coordinate """
from django.db.models import QuerySet

from library.genomics import format_chrom
from snpdb.models import GenomeBuild, Variant, VariantCoordinate
from upload.models import ModifiedImportedVariant


def get_results_from_variant_coordinate(genome_build: GenomeBuild, qs: QuerySet, variant_coordinate: VariantCoordinate,
                                        any_alt: bool = False) -> QuerySet[Variant]:
    """
    :param genome_build: genome build (used for format variant_coordinate.chromosome for variant search
    :param qs: A query set that we'll be searching inside of (including variants found via ModifiedImportedVariant)
    :param variant_coordinate: The variant coordinate to lookup
    :param any_alt: If true, search without using alt and return all matches
    :return: A QuerySet of variants
    """
    chrom = format_chrom(variant_coordinate.chrom, genome_build.reference_fasta_has_chr)
    results = qs.filter(Variant.get_chrom_q(chrom),
                        locus__position=variant_coordinate.position,
                        locus__ref__seq=variant_coordinate.ref)
    if not any_alt:
        results = results.filter(alt__seq=variant_coordinate.alt, svlen=variant_coordinate.svlen)

    if not results.exists():
        if not any_alt:
            return ModifiedImportedVariant.get_variants_for_unnormalized_variant(qs, variant_coordinate)
        else:
            # should we really be searching ModifiedImportVariants with any alt? or should that just happen for
            # the filter of "real" variants
            return ModifiedImportedVariant.get_variants_for_unnormalized_variant_any_alt(qs, variant_coordinate)
    return results
