from collections import defaultdict
from typing import Tuple, List

from analysis.models.nodes.analysis_node import Analysis
from analysis.models.nodes.sources import SampleNode, CohortNode, TrioNode, PedigreeNode


def sort_analyses_by_date_and_merge_details(all_analysis_details) -> List[Tuple[Analysis, str]]:
    sorted_analyses = sorted(all_analysis_details.keys(), key=lambda x: x.created)

    analysis_details = []
    for analysis in sorted_analyses:
        details = ','.join(sorted(all_analysis_details[analysis]))
        analysis_details.append((analysis, details))

    return analysis_details


def get_related_analysis_details_for_samples(user, samples) -> List[Tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = Analysis.filter_for_user(user).values_list("pk", flat=True)

    for sn in SampleNode.objects.filter(analysis__in=analyses_ids,
                                        sample__in=samples).select_related("analysis", "sample"):
        all_analysis_details[sn.analysis].add(sn.sample.name)

    return sort_analyses_by_date_and_merge_details(all_analysis_details)


def get_related_analysis_details_for_cohort(user, cohorts) -> List[Tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = Analysis.filter_for_user(user).values_list("pk", flat=True)

    for cohort_node in CohortNode.objects.filter(analysis__in=analyses_ids,
                                                 cohort__in=cohorts).select_related("analysis", "cohort"):
        all_analysis_details[cohort_node.analysis].add(str(cohort_node.cohort))
    return sort_analyses_by_date_and_merge_details(all_analysis_details)


def get_related_analysis_details_for_trio(user, trios) -> List[Tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = Analysis.filter_for_user(user).values_list("pk", flat=True)
    for trio_node in TrioNode.objects.filter(analysis__in=analyses_ids,
                                             trio__in=trios).select_related("analysis", "trio"):
        all_analysis_details[trio_node.analysis].add(str(trio_node.trio))

    return sort_analyses_by_date_and_merge_details(all_analysis_details)


def get_related_analysis_details_for_pedigree(user, pedigrees) -> List[Tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = Analysis.filter_for_user(user).values_list("pk", flat=True)
    for pedigree_node in PedigreeNode.objects.filter(analysis__in=analyses_ids,
                                                     pedigree__in=pedigrees).select_related("analysis", "pedigree"):
        all_analysis_details[pedigree_node.analysis].add(str(pedigree_node.pedigree))

    return sort_analyses_by_date_and_merge_details(all_analysis_details)
