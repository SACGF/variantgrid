from django.db.models.query_utils import Q
from functools import reduce
import operator

CONTROL_POPULATION_DATABASES = ["gnomad_af", "af_1kg", "af_uk10k", "topmed_af"]


def get_population_af_q(max_allele_frequency, population_databases=None, group_operation=operator.and_):
    """ population_databases : container of strings of databases to use
        (default: CONTROL_POPULATION_DATABASES) """

    if population_databases is None:
        population_databases = CONTROL_POPULATION_DATABASES

    filters = []
    for field in population_databases:
        q_isnull = Q(**{f"variantannotation__{field}__isnull": True})
        q_max_value = Q(**{f"variantannotation__{field}__lte": max_allele_frequency})
        filters.append(q_isnull | q_max_value)

    return reduce(group_operation, filters)
