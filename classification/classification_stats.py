from collections import Counter, defaultdict
from collections.abc import Iterable
from typing import Any, Optional

import numpy as np
from django.conf import settings
from django.contrib.auth.models import User
from django.db import connection
from django.db.models import Count, F, QuerySet
from django.db.models.fields.json import KT

from classification.enums.classification_enums import CriteriaEvaluation, ShareLevel
from classification.models import EvidenceKeyMap
from classification.models.classification import ClassificationModification
from library.django_utils import get_field_counts
from snpdb.enums import AlleleOriginBucket, ClinicalSignificance
from snpdb.models import Lab

ALLELE_PATH = "classification__allele_info__allele_id"
# Which record represents an allele (or a gene + allele) when several labs have classified it
LATEST_FIRST = ("-created", "-pk")


def get_classification_counts(
        user: User,
        show_unclassified=True,
        unique_alleles=False,
        allele_origin_buckets: Optional[set[AlleleOriginBucket]] = None) -> dict[str, int]:
    qs = get_visible_classifications_qs(user, allele_origin_buckets=allele_origin_buckets)

    field_counts: dict[str, int]
    if unique_alleles:
        keys = set(ClinicalSignificance.LABELS.keys())
        if not show_unclassified:
            keys.remove(None)
        field_counts = get_classification_counts_allele(qs, "clinical_significance", possible_values=keys)
    else:
        field_counts = get_field_counts(qs, "clinical_significance")

    # Add all entries as empty so colors are in the right order in graph.
    classification_counts = {}
    for cs, label in ClinicalSignificance.LABELS.items():
        if cs or show_unclassified:
            classification_counts[label] = 0

    for clinical_significance, clinical_significance_count in field_counts.items():
        clinical_significance_label = ClinicalSignificance.LABELS[clinical_significance]
        classification_counts[clinical_significance_label] = clinical_significance_count
    return classification_counts


def get_classification_counts_allele(qs: QuerySet[ClassificationModification], field: str, possible_values: Iterable[Any]) -> dict[Any, int]:
    """ Counts of field over the latest record per allele, so an allele classified differently by 2 labs counts once """
    latest_per_allele_qs = qs.filter(classification__allele_info__allele__isnull=False) \
        .order_by(ALLELE_PATH, *LATEST_FIRST).distinct(ALLELE_PATH)
    counts_qs = ClassificationModification.objects.filter(pk__in=latest_per_allele_qs.values("pk")) \
        .values_list(field).annotate(count=Count("pk")).order_by()
    field_counts = dict(counts_qs)
    return {value: field_counts.get(value, 0) for value in possible_values}


def get_visible_classifications_qs(user: User, allele_origin_buckets: Optional[set[AlleleOriginBucket]] = None) -> QuerySet[ClassificationModification]:
    # now excludes external labs (don't want to report on Shariant only data within SA Path for example)

    if settings.CLASSIFICATION_STATS_USE_SHARED:
        # Records shared at a discordant level are readable by every user (ShareLevel.has_access), so skip the
        # guardian lookup latest_for_user would do - it is most of the cost of these queries
        qs = ClassificationModification.objects.filter(
            is_last_published=True,
            share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
            classification__withdrawn=False,
            classification__lab__external=False,
            clinical_significance__isnull=False,
        )
        if allele_origin_buckets:
            qs = qs.filter(classification__allele_origin_bucket__in=allele_origin_buckets)
        return qs

    return ClassificationModification.latest_for_user(
        user,
        exclude_withdrawn=True,
        allele_origin_buckets=allele_origin_buckets,
        classification__lab__external=False,
        classification__variant__isnull=False,
    )


def get_grouped_classification_counts(user: User,
                                      field: str,
                                      evidence_key: Optional[str] = None,
                                      field_labels: Optional[dict[str, str]] = None,
                                      max_groups=10,
                                      show_unclassified=True,
                                      norm_factor: Optional[dict[str, float]] = None,
                                      allele_level: bool = False,
                                      allele_origin_buckets: Optional[set[AlleleOriginBucket]] = None) -> list[dict[str, dict]]:
    """ :param user: User used to check visibility of classifications
        :param field: the value we're extracting from evidence to group on (from Classification)
        :param evidence_key: label from ekey lookup
        :param field_labels: labels from dict (falls back on raw value if keys not present)
        mutually exclusive with evidence_key
        :param max_groups: the max size of the returns list
        :param show_unclassified: return counts for clinical_significance=None (as 'Unclassified')
        :param allele_level: if True only the latest record per allele within each group is counted
        :return: data for graphing, a list of dicts with x,y,name and type
    """
    if evidence_key and field_labels:
        raise ValueError("Can't supply both 'evidence_key' and 'field_labels'")

    group_value = KT(f"{field}__{evidence_key}__value") if evidence_key else F(field)
    vc_qs = get_visible_classifications_qs(user, allele_origin_buckets).annotate(group_value=group_value)
    if allele_level:
        rows_sql, params = vc_qs.values_list("pk", "created", "clinical_significance", ALLELE_PATH, "group_value") \
            .order_by().query.sql_with_params()
        # OFFSET 0 keeps the rows subquery from being flattened, so an evidence value is extracted once per record
        # rather than again for each of the DISTINCT ON and GROUP BY that use it
        group_counts_sql = f"""
            SELECT clinical_significance, group_value, COUNT(*) FROM (
                SELECT DISTINCT ON (group_value, allele_id) * FROM ({rows_sql} OFFSET 0)
                    visible (id, created, clinical_significance, allele_id, group_value)
                WHERE allele_id IS NOT NULL ORDER BY group_value, allele_id, created DESC, id DESC
            ) latest_per_allele GROUP BY 1, 2"""
        with connection.cursor() as cursor:
            cursor.execute(group_counts_sql, params)
            group_counts_rows = cursor.fetchall()
    else:
        group_counts_rows = vc_qs.values_list("clinical_significance", "group_value") \
            .annotate(count=Count("pk")).order_by()

    counts = Counter()
    classification_counts = defaultdict(Counter)
    for clinical_significance, value, count in group_counts_rows:
        if value is None:  # filtered here rather than in SQL, where it would extract the evidence value again
            continue
        if field_labels:
            value = field_labels.get(value, value)
        counts[value] += count
        classification_counts[clinical_significance][value] += count

    if norm_factor:
        group_counts = {}
        for g, count in counts.items():
            if nf := norm_factor.get(g):
                group_counts[g] = nf * count
    else:
        group_counts = counts

    top_groups = sorted(group_counts, key=lambda g: (-group_counts[g], g))[:max_groups]

    data = []
    for cs, clinical_significance_label in ClinicalSignificance.LABELS.items():
        if cs or show_unclassified:
            counts = classification_counts[cs]
            if norm_factor:
                x = []
                y = []
                for g in top_groups:
                    if nf := norm_factor.get(g):
                        x.append(g)
                        y.append(nf * counts[g])
            else:
                x = top_groups
                y = [counts[i] for i in top_groups]
            data.append({"x": x,
                         "y": y,
                         "name": clinical_significance_label,
                         "type": 'bar'})
    return data


def get_criteria_counts(
        user: User,
        evidence_field: str,
        critera="acmg",
        allele_origin_buckets: Optional[set[AlleleOriginBucket]] = None
    ) -> dict[str, list[dict]]:
    """ Per clinical significance, the fraction of records where each criterion is met (CriteriaEvaluation.is_met),
        not met (NM, NA or N), or not set """
    acmg_labels = {e.key: e.pretty_label for e in EvidenceKeyMap.instance().criteria_for(critera)}
    criteria_index = {key: i for i, key in enumerate(acmg_labels)}

    vc_qs = get_visible_classifications_qs(user, allele_origin_buckets=allele_origin_buckets)
    vc_qs = vc_qs.filter(clinical_significance__isnull=False).order_by()
    totals = dict(vc_qs.values_list("clinical_significance").annotate(count=Count("pk")))

    # jsonb_each expands each record's evidence once - a -> per criterion would read the whole blob again each time
    rows_sql, params = vc_qs.values_list("clinical_significance", evidence_field).query.sql_with_params()
    value_counts_sql = f"""
        SELECT visible.clinical_significance, evidence.key,
               CASE WHEN jsonb_typeof(evidence.value) = 'object' THEN evidence.value ->> 'value'
                    ELSE evidence.value #>> '{{}}' END AS criteria_value,
               COUNT(*)
        FROM ({rows_sql}) visible (clinical_significance, evidence)
        CROSS JOIN LATERAL jsonb_each(visible.evidence) evidence
        WHERE evidence.key = ANY(%s)
        GROUP BY 1, 2, 3"""
    with connection.cursor() as cursor:
        cursor.execute(value_counts_sql, [*params, list(acmg_labels)])
        value_counts = cursor.fetchall()

    n = len(acmg_labels)
    met_not_met = defaultdict(lambda: (np.zeros(n), np.zeros(n)))
    for clinical_significance, key, value, count in value_counts:
        if value:
            met, not_met = met_not_met[clinical_significance]
            (met if CriteriaEvaluation.is_met(value) else not_met)[criteria_index[key]] += count

    acmg_data = {}
    labels = list(acmg_labels.values())
    for clinical_significance, total in totals.items():
        met, not_met = (counts / total for counts in met_not_met[clinical_significance])
        not_set = 1 - met - not_met

        acmg_data[clinical_significance] = [
            {"x": labels,
             "y": met.tolist(),
             "name": "Met",
             "type": 'bar',
             "marker": {
                 "color": "rgb(80,200,80)"
             }},
            {"x": labels,
             "y": not_met.tolist(),
             "name": "Not met",
             "type": 'bar',
             "marker": {
                 "color": "rgb(140,140,140)"
             }},
            {"x": labels,
             "y": not_set.tolist(),
             "name": "Not set",
             "type": 'bar',
             "marker": {
                 "color": "rgb(220,220,220)"
             }},
        ]

    return acmg_data


def get_lab_gene_counts(user: User, lab: Lab):
    evidence_field = "published_evidence"
    lab_field = "classification__lab"

    gene_symbol_field = f"{evidence_field}__gene_symbol__value"
    kwargs = {gene_symbol_field + "__isnull": False, lab_field: lab}
    vc_qs = get_visible_classifications_qs(user).filter(**kwargs)
    gene_clinical_significance_counts = defaultdict(Counter)

    for gene_symbol, clinical_significance in vc_qs.values_list(gene_symbol_field, "clinical_significance"):
        gene_clinical_significance_counts[gene_symbol][clinical_significance] += 1

    return gene_clinical_significance_counts
