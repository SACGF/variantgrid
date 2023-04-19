# def search_cohort(search_string: str, user: User, genome_build: GenomeBuild, **kwargs) -> Iterable[Sample]:
#     return Cohort.filter_for_user(user).filter(name__icontains=search_string,
#                                                genome_build=genome_build)