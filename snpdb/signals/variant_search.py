import re

from django.db.models import QuerySet

from genes.hgvs import HGVSMatcher
from library.genomics import format_chrom
from snpdb.clingen_allele import get_clingen_allele
from snpdb.models import Variant, LOCUS_PATTERN, LOCUS_NO_REF_PATTERN, DbSNP, DBSNP_PATTERN, VariantCoordinate, \
    ClinGenAllele
from snpdb.search2 import search_receiver, SearchInputInstance, SearchExample

COSMIC_PATTERN = re.compile(r"^(COS[VM]).*$", re.IGNORECASE)


@search_receiver(
    search_type=Variant,
    pattern=COSMIC_PATTERN,
    sub_name="COSMIC",
    example=SearchExample(
        note="Provide the COSMIC ID",
        example="COSV53567516"
    )
)
def variant_cosmic_search(search_input: SearchInputInstance):
    for genome_build in search_input.genome_builds:
        variant_qs = search_input.get_visible_variants(genome_build)
        if search_input.match.group(1) == "COSV":
            yield variant_qs.filter(variantannotation__cosmic_id__icontains=search_input.search_string)
        elif search_input.match.group(1) == "COSM":
            yield variant_qs.filter(variantannotation__cosmic_legacy_id__icontains=search_input.search_string)


@search_receiver(
    search_type=Variant,
    pattern=LOCUS_NO_REF_PATTERN,
    sub_name="Locus w/out Ref",
    example=SearchExample(
        note="Provide chromosome and position",
        example="chr1:169519049"
    )
)
def search_variant_locus_no_ref(search_input: SearchInputInstance):
    for genome_build in search_input.genome_builds:
        chrom, position = search_input.match.groups()
        chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        yield search_input.get_visible_variants(genome_build).filter(
            locus__contig__name=chrom,
            locus__position=position
        )


@search_receiver(
    search_type=Variant,
    pattern=LOCUS_PATTERN,
    sub_name="Locus with Ref",
    example=SearchExample(
        note="TO DO provide an example with ref",
        example="chr1:169519049"
    )
)
def search_variant_locus_with_ref(search_input: SearchInputInstance):
    for genome_build in search_input.genome_builds:
        chrom, position, ref = search_input.match.groups()
        chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        yield search_input.get_visible_variants(genome_build).filter(
            locus__contig__name=chrom,
            locus__position=position,
            locus__ref__seq=ref
        )


@search_receiver(
    search_type=Variant,
    pattern=ClinGenAllele.CLINGEN_ALLELE_CODE_PATTERN,
    sub_name="ClinGen Allele ID",
    example=SearchExample(note="ClinGen Allele ID", example="CA285410130")
)
def allele_search(search_input: SearchInputInstance):
    search_string = search_input.search_string
    if ClinGenAllele.looks_like_id(search_string):
        allele = get_clingen_allele(search_string).allele
        for genome_build in search_input.genome_builds:
            try:
                if variant := allele.variant_for_build(genome_build):
                    yield variant
            except ValueError:
                pass # need to re-add support for creating variants
        #yield allele
        # else:
        #     # FIXME none of this will work, neither Variant or CreateManualVariant is previewable
        #     # and have that allele just take you to the variant page
        #     for genome_build in search_input.genome_builds:
        #         variant_qs = variant_qs.filter(variantallele__allele__clingen_allele=clingen_allele,
        #                                        variantallele__genome_build=genome_build)
        #         if variant_qs.exists():
        #             yield variant_qs
        #         else:
        #             if can_create_variants(search_input.user):
        #                 variant_string = clingen_allele.get_variant_string(genome_build)
        #                 variant_string_abbreviated = clingen_allele.get_variant_string(genome_build, abbreviate=True)
        #                 search_message = f"'{clingen_allele}' resolved to '{variant_string_abbreviated}'"
        #                 # FIXME, this isn't previewable
        #                 yield CreateManualVariant(genome_build, variant_string), search_message


def get_results_from_variant_tuples(qs: QuerySet, data: VariantCoordinate, any_alt: bool = False):
    """
    :param qs: A query set that we'll be searching inside of (except for when returning ModifiedImportVariants)
    :param data: The variant coordinate to lookup
    :param any_alt: If true, search without using alt and return all matches
    :return: A QuerySet of variants
    """
    (chrom, position, ref, alt) = data
    position = int(position)

    results = qs.filter(Variant.get_chrom_q(chrom), locus__position=position, locus__ref__seq=ref)
    if not any_alt:
        results = results.filter(alt__seq=alt)

    if not results:
        if not any_alt:
            pass
            # results = ModifiedImportedVariant.get_variants_for_unnormalized_variant(chrom, position, ref, alt)
        else:
            pass
            # should we really be searching ModifiedImportVariants with any alt? or should that just happen for
            # the filter of "real" variants
            # results = ModifiedImportedVariant.get_variants_for_unnormalized_variant_any_alt(chrom, position, ref)
    return results


@search_receiver(
    search_type=Variant,
    pattern=DBSNP_PATTERN,
    sub_name="dbSNP ID",
    example=SearchExample(
        note="Provide the variant's dbSNP ID",
        example="rs6025"
    )
)
def search_variant_db_snp(search_input: SearchInputInstance):
    # Do via API as a full table scan takes way too long with big data
    dbsnp = DbSNP.get(search_input.search_string)

    for genome_build in search_input.genome_builds:
        matcher = HGVSMatcher(genome_build)
        for data in dbsnp.get_alleles_for_genome_build(genome_build):
            if hgvs_string := data.get("hgvs"):
                dbsnp_message = f"dbSNP '{search_input.search_string}' resolved to '{hgvs_string}'"
                variant_tuple = matcher.get_variant_tuple(hgvs_string)
                results = get_results_from_variant_tuples(search_input.get_visible_variants(genome_build), variant_tuple)
                if results.exists():
                    for r in results:
                        yield r, dbsnp_message
                else:
                    variant_string = Variant.format_tuple(*variant_tuple)
                    # FIXME return CreateManualVariant again
                    #search_results.append(SearchResult(CreateManualVariant(genome_build, variant_string),
                    #                                   message=dbsnp_message))