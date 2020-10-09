from typing import Set

from django.db.models import QuerySet
from django.db.models.query_utils import Q
from functools import reduce
import operator

from annotation.annotation_version_querysets import get_variant_queryset_for_latest_annotation_version
from genes.models import GeneSymbol, Gene
from snpdb.models import VariantZygosityCountCollection
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import VariantAllele


def get_variant_queryset_for_gene_symbol(gene_symbol: GeneSymbol, genome_build: GenomeBuild, traverse_aliases: bool = False):
    """

    :param gene_symbol: Gene Symbol to search on
    :param genome_build: Genome Build
    :param traverse_aliases: If true, will see if the gene_symbol has aliases, will not
    :return:
    """
    # See issue #2449 - people often search for a gene symbol that exists in a different genome build
    # so we need to retrieve all the genes that have ever been associated with the symbol
    genes: Set[Gene]
    if traverse_aliases:
        genes = gene_symbol.alias_meta.genes
    else:
        genes = set(gene_symbol.genes)

    qs = get_variant_queryset_for_latest_annotation_version(genome_build)
    return qs.filter(variantgeneoverlap__gene__in=genes)


def get_has_classifications_q(genome_build):
    va_build_qs = VariantAllele.objects.filter(genome_build=genome_build,
                                               allele__variantallele__variant__classification__isnull=False)
    return Q(variantallele__in=va_build_qs)


def get_has_variant_tags():
    return Q(varianttag__isnull=False)


def variant_qs_filter_has_internal_data(variant_qs: QuerySet, genome_build: GenomeBuild) -> QuerySet:
    qs, count_column = VariantZygosityCountCollection.annotate_global_germline_counts(variant_qs)
    interesting = [Q(**{f"{count_column}__gt": 0}),
                   get_has_classifications_q(genome_build),
                   get_has_variant_tags()]
    q = reduce(operator.or_, interesting)
    return qs.filter(q)
