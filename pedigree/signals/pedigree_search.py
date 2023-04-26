from pedigree.models import Pedigree
from snpdb.search import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=Pedigree,
    example=SearchExample(note="Search of the pedigree's name")
)
def search_pedigree(search_input: SearchInputInstance):
    return Pedigree.filter_for_user(search_input.user).filter(search_input.q_words(), cohort__genome_build__in=search_input.genome_builds)
