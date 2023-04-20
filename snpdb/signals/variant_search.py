import re

from snpdb.models import Variant
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
