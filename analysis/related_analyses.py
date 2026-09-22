from collections import defaultdict

from django.db.models import Q

from analysis.models.nodes.analysis_node import Analysis
from analysis.models.nodes.sources import CohortNode, PedigreeNode, SampleNode, TrioNode
from analysis.models.nodes.sources.duo_node import DuoNode
from analysis.models.nodes.sources.quad_node import QuadNode
from patients.models_enums import SampleSourceLevel
from snpdb.models import Sample


def _related_analyses_ids(user):
    """ A template's source node keeps whatever sample it was built from as an example - that doesn't
        make the template related to the sample """
    return Analysis.filter_for_user(user).filter(template_type__isnull=True).values_list("pk", flat=True)


def sort_analyses_by_date_and_merge_details(all_analysis_details) -> list[tuple[Analysis, str]]:
    sorted_analyses = sorted(all_analysis_details.keys(), key=lambda x: x.created)

    analysis_details = []
    for analysis in sorted_analyses:
        details = ','.join(sorted(all_analysis_details[analysis]))
        analysis_details.append((analysis, details))

    return analysis_details


def _sample_node_q(samples) -> Q:
    """ A node on one of the samples, or on a patient / specimen / extraction the samples belong to - a
        group level node leaves sample null and resolves its samples at query time """
    q = Q(source_level=SampleSourceLevel.SAMPLE, sample__in=samples)
    group_ids = defaultdict(set)
    sample_qs = Sample.objects.filter(pk__in=[s.pk for s in samples])
    for extraction_id, specimen_id, patient_id, specimen_patient_id in sample_qs.values_list(
            "extraction", "extraction__specimen", "patient", "extraction__specimen__patient"):
        group_ids["extraction"].add(extraction_id)
        group_ids["specimen"].add(specimen_id)
        group_ids["patient"].update((patient_id, specimen_patient_id))

    for level, field in SampleNode.SOURCE_LEVEL_FIELDS.items():
        if ids := group_ids.get(field, set()) - {None}:
            q |= Q(source_level=level, **{f"{field}__in": ids})
    return q


def get_related_analysis_details_for_samples(user, samples) -> list[tuple[Analysis, str]]:
    """ details are labelled by level, eg "Sample: proband, Patient: Smith, J" """
    names_by_analysis_and_level = defaultdict(lambda: defaultdict(set))
    analyses_ids = _related_analyses_ids(user)

    sample_node_qs = SampleNode.objects.filter(_sample_node_q(samples), analysis__in=analyses_ids)
    for sn in sample_node_qs.select_related("analysis", "sample", "extraction__specimen", "specimen", "patient"):
        source = sn.get_source_object()
        name = source.name if sn.source_level == SampleSourceLevel.SAMPLE else str(source)
        names_by_analysis_and_level[sn.analysis][sn.source_level].add(name)

    analysis_details = []
    for analysis in sorted(names_by_analysis_and_level, key=lambda x: x.created):
        names_by_level = names_by_analysis_and_level[analysis]
        details = ", ".join(f"{SampleSourceLevel(level).label}: {','.join(sorted(names_by_level[level]))}"
                            for level in SampleNode.SOURCE_LEVEL_FIELDS if level in names_by_level)
        analysis_details.append((analysis, details))
    return analysis_details


def get_related_analysis_details_for_cohort(user, cohorts) -> list[tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = _related_analyses_ids(user)

    for cohort_node in CohortNode.objects.filter(analysis__in=analyses_ids,
                                                 cohort__in=cohorts).select_related("analysis", "cohort"):
        all_analysis_details[cohort_node.analysis].add(str(cohort_node.cohort))
    return sort_analyses_by_date_and_merge_details(all_analysis_details)


def get_related_analysis_details_for_trio(user, trios) -> list[tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = _related_analyses_ids(user)
    for trio_node in TrioNode.objects.filter(analysis__in=analyses_ids,
                                             trio__in=trios).select_related("analysis", "trio"):
        all_analysis_details[trio_node.analysis].add(str(trio_node.trio))

    return sort_analyses_by_date_and_merge_details(all_analysis_details)


def get_related_analysis_details_for_quad(user, quads) -> list[tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = _related_analyses_ids(user)
    for quad_node in QuadNode.objects.filter(analysis__in=analyses_ids,
                                             quad__in=quads).select_related("analysis", "quad"):
        all_analysis_details[quad_node.analysis].add(str(quad_node.quad))
    return sort_analyses_by_date_and_merge_details(all_analysis_details)


def get_related_analysis_details_for_duo(user, duos) -> list[tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = _related_analyses_ids(user)
    for duo_node in DuoNode.objects.filter(analysis__in=analyses_ids,
                                           duo__in=duos).select_related("analysis", "duo"):
        all_analysis_details[duo_node.analysis].add(str(duo_node.duo))
    return sort_analyses_by_date_and_merge_details(all_analysis_details)


def get_related_analysis_details_for_pedigree(user, pedigrees) -> list[tuple[Analysis, str]]:
    all_analysis_details = defaultdict(set)
    analyses_ids = _related_analyses_ids(user)
    for pedigree_node in PedigreeNode.objects.filter(analysis__in=analyses_ids,
                                                     pedigree__in=pedigrees).select_related("analysis", "pedigree"):
        all_analysis_details[pedigree_node.analysis].add(str(pedigree_node.pedigree))

    return sort_analyses_by_date_and_merge_details(all_analysis_details)
