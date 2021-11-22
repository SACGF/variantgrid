import operator
from collections import Counter, defaultdict
from typing import Dict, List, Optional

import numpy as np
from django.conf import settings
from django.contrib.auth.models import User

from classification.enums import ClinicalSignificance
from classification.enums.classification_enums import CriteriaEvaluation
from classification.models import EvidenceKeyMap
from classification.models.classification import Classification, ClassificationModification
from library.django_utils import get_field_counts
from snpdb.models import Lab


def get_classification_counts(user: User, show_unclassified=True) -> Dict[str, int]:
    qs = get_visible_classifications_qs(user)
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


def get_visible_classifications_qs(user: User):
    if settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED:
        qs = ClassificationModification.latest_for_user(
            user,
            published=True,
            exclude_withdrawn=True,
            shared_only=True,
            clinical_significance__isnull=False
        ).exclude(classification__variant__isnull=True)
    else:
        qs = Classification.filter_for_user(user)

    return qs


def get_grouped_classification_counts(user: User,
                                      field: str,
                                      evidence_key: Optional[str] = None,
                                      field_labels: Optional[Dict[str, str]] = None,
                                      max_groups=10,
                                      show_unclassified=True,
                                      norm_factor: Dict[str, float] = None) -> List[Dict[str, Dict]]:
    """ :param user: User used to check visibility of classifications
        :param field: the value we're extracting from evidence to group on (from Classification)
        :param evidence_key: label from ekey lookup
        :param field_labels: labels from dict (falls back on raw value if keys not present)
        mutually exclusive with evidence_key
        :param max_groups: the max size of the returns list
        :param show_unclassified: return counts for clinical_significance=None (as 'Unclassified')
        :return: data for graphing, a list of dicts with x,y,name and type
    """
    if evidence_key and field_labels:
        raise ValueError("Can't supply both 'evidence_key' and 'field_labels'")

    vc_qs = get_visible_classifications_qs(user)
    values_qs = vc_qs.values_list("clinical_significance", field)

    counts = Counter()
    classification_counts = defaultdict(Counter)
    for clinical_significance, field in values_qs:
        if evidence_key:
            # field will either be evidence or published evidence
            value = Classification.get_optional_value_from(field, evidence_key)
        elif field_labels:
            value = field_labels.get(field, field)
        else:
            value = field

        if value is not None:
            counts[value] += 1
            classification_counts[clinical_significance][value] += 1

    if norm_factor:
        group_counts = {}
        for g, count in counts.items():
            if nf := norm_factor.get(g):
                group_counts[g] = nf * count
    else:
        group_counts = counts

    top_groups = [gc[0] for gc in sorted(group_counts.items(), key=operator.itemgetter(1), reverse=True)][:max_groups]

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


def get_criteria_counts(user: User, evidence_field: str) -> Dict[str, List[Dict]]:
    acmg_labels = dict((e.key, e.pretty_label) for e in EvidenceKeyMap.instance().acmg_criteria())
    n = len(acmg_labels)
    acmg_met_not_met_by_significance = defaultdict(lambda: (np.zeros(n), np.zeros(n), np.zeros(n)))
    total_clinical_significance = Counter()

    vc_qs = get_visible_classifications_qs(user)
    vc_qs = vc_qs.filter(clinical_significance__isnull=False)
    num_classifications = 0
    for clinical_significance, evidence in vc_qs.values_list("clinical_significance", evidence_field):
        num_classifications += 1
        total_clinical_significance[clinical_significance] += 1
        met, not_met, not_set = acmg_met_not_met_by_significance[clinical_significance]
        for i, k in enumerate(acmg_labels):
            value = Classification.get_optional_value_from(evidence, k)
            if value:
                if CriteriaEvaluation.is_met(value):
                    met[i] += 1
                else:
                    not_met[i] += 1

            else:
                not_set[i] += 1

    acmg_data = {}
    labels = list(acmg_labels.values())
    for clinical_significance, (met, not_met, not_set) in acmg_met_not_met_by_significance.items():
        total = total_clinical_significance[clinical_significance]
        acmg_met_not_met_by_significance[clinical_significance] = (met/total, not_met/total, not_set/total)

        acmg_data[clinical_significance] = [
            {"x": labels,
             "y": (met / total).tolist(),
             "name": "Met",
             "type": 'bar',
             "marker": {
                 "color": "rgb(80,200,80)"
             }},
            {"x": labels,
             "y": (not_met / total).tolist(),
             "name": "Not met",
             "type": 'bar',
             "marker": {
                 "color": "rgb(140,140,140)"
             }},
            {"x": labels,
             "y": (not_set / total).tolist(),
             "name": "Not set",
             "type": 'bar',
             "marker": {
                 "color": "rgb(220,220,220)"
             }},
        ]

    return acmg_data


def get_lab_gene_counts(user: User, lab: Lab):
    if settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED:
        evidence_field = "published_evidence"
        lab_field = "classification__lab"
    else:
        evidence_field = "evidence"
        lab_field = "lab"

    gene_symbol_field = f"{evidence_field}__gene_symbol__value"
    kwargs = {gene_symbol_field + "__isnull": False, lab_field: lab}
    vc_qs = get_visible_classifications_qs(user).filter(**kwargs)
    gene_clinical_significance_counts = defaultdict(Counter)

    for gene_symbol, clinical_significance in vc_qs.values_list(gene_symbol_field, "clinical_significance"):
        gene_clinical_significance_counts[gene_symbol][clinical_significance] += 1

    return gene_clinical_significance_counts
