from snpdb.models import VCF
from snpdb.search import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=VCF,
    example=SearchExample(note="Provide part of the VCF's name")
)
def vcf_search(search_input: SearchInputInstance):
    yield VCF.filter_for_user(search_input.user).filter(
        search_input.q_words(),
        genome_build__in=search_input.genome_builds
    )
