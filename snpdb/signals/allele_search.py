from snpdb.clingen_allele import get_clingen_allele
from snpdb.models import ClinGenAllele, Allele
from snpdb.search2 import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=Allele,
    pattern=ClinGenAllele.CLINGEN_ALLELE_CODE_PATTERN,
    example=SearchExample(note="ClinGen Allele ID", example="CA285410130")
)
def allele_search(search_input: SearchInputInstance):
    search_string = search_input.search_string
    if ClinGenAllele.looks_like_id(search_string):
        clingen_allele = get_clingen_allele(search_string)
        yield clingen_allele.allele

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
