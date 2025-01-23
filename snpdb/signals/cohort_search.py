from snpdb.models import Cohort
from snpdb.search import search_receiver, SearchInputInstance, SearchExample


@search_receiver(
    search_type=Cohort,
    example=SearchExample(
        note="Provide part of the cohort's name"
    )
)
def search_cohort(search_input: SearchInputInstance):
    q = search_input.q_words()
    qs = Cohort.filter_for_user(search_input.user).filter(q,
                                                          genome_build__in=search_input.genome_builds,
                                                          vcf__isnull=True)  # These will show up in VCF search
    yield qs
