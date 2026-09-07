"""
The gene content box on the classification form - for records the create-from-variant page never saw:
made from the Classify & Report tab with nothing to copy, imported, or created before the gene was curated.

It offers the same deduplicated gene candidates as the create page
(classification/models/classification.py:ClassificationConsensus.gene_consensus_groups), and applying one
patches the empty fields as SubmissionSource.CONSENSUS.

Entry points: gene_consensus_panel (the card, its dialog rows, and the apply POST), reached at
classification_gene_consensus.
"""
from dataclasses import dataclass
from datetime import date, datetime
from typing import Any, Optional, Union

from django.core.exceptions import PermissionDenied
from django.http import Http404
from django.http.request import HttpRequest
from django.http.response import HttpResponseBase, JsonResponse
from django.shortcuts import get_object_or_404, render

from classification.enums import AlleleOriginBucket, SpecialEKeys, SubmissionSource
from classification.models import (
    Classification,
    ClassificationConsensus,
    ClassificationModification,
    GeneConsensusGroup,
)
from classification.models.classification import COPY_SCOPES_GENE
from classification.models.evidence_key import EvidenceKeyMap


@dataclass(frozen=True)
class GeneContentValue:
    """ One gene scope key: what this record has now, and what the candidate would bring """
    label: str
    current: Any
    candidate: Any


@dataclass(frozen=True)
class GeneContentComparison:
    """ A candidate row in the dialog - the group, and its content beside this record's own """
    group: GeneConsensusGroup
    values: list[GeneContentValue]


def _as_date(value: Union[date, datetime]) -> date:
    return value.date() if isinstance(value, datetime) else value


def classification_gene_symbol(classification: Classification) -> Optional[str]:
    """ The gene the record is about - what it says itself, else what its allele resolved to """
    if gene_symbol := classification.get(SpecialEKeys.GENE_SYMBOL):
        return gene_symbol
    if allele_info := classification.allele_info:
        if gene_symbols := allele_info.gene_symbols:
            return str(gene_symbols[0])
    return None


def gene_content_applied_at(classification: Classification,
                            allele_origin_bucket: Optional[AlleleOriginBucket]) -> Optional[datetime]:
    """ When gene content was last copied into this record. There is no per-value provenance, so the
        CONSENSUS modification that touched a gene scope key is the record of it """
    gene_scope_keys = set(ClassificationConsensus.gene_scope_keys(allele_origin_bucket))
    modifications = ClassificationModification.objects.filter(classification=classification,
                                                              source=SubmissionSource.CONSENSUS).order_by("-created")
    for modification in modifications:
        if gene_scope_keys.intersection(modification.delta or {}):
            return modification.created
    return None


def newer_than(groups: list[GeneConsensusGroup], applied_at: datetime) -> list[GeneConsensusGroup]:
    """ The gene has been re-curated since this record copied from it """
    applied_date = _as_date(applied_at)
    return [group for group in groups if _as_date(group.representative.curated_date) > applied_date]


def gene_consensus_panel(request: HttpRequest, classification_id: int) -> HttpResponseBase:
    """ The gene content box: the card on the classification form, the rows its dialog shows (?rows=1), and
        the POST that copies a chosen record's gene content into this one """
    classification = get_object_or_404(Classification, pk=classification_id)
    if not classification.can_write(request.user):
        raise PermissionDenied(f"You do not have WRITE permission on classification {classification_id}")

    allele_origin_bucket = AlleleOriginBucket(classification.allele_origin_bucket)
    gene_symbol = classification_gene_symbol(classification)
    if not gene_symbol:
        raise Http404(f"Classification {classification_id} has no gene symbol")

    if request.method == "POST":
        return _apply_gene_consensus(request, classification)

    groups = ClassificationConsensus.gene_consensus_groups(gene_symbol=gene_symbol, user=request.user,
                                                           allele_origin_bucket=allele_origin_bucket,
                                                           exclude_classification=classification)
    context = {
        "classification": classification,
        "gene_symbol": gene_symbol,
        "groups": groups,
        "record_count": sum(group.record_count for group in groups),
    }
    if request.GET.get("rows"):
        context["comparisons"] = _comparisons(classification, groups, allele_origin_bucket)
        return render(request, 'classification/gene_consensus_rows.html', context)

    # The box only takes space while it is useful - once gene content has been copied in it collapses to a
    # line, and only while the gene has been curated again since
    applied_at = gene_content_applied_at(classification, allele_origin_bucket)
    if applied_at:
        groups = newer_than(groups, applied_at)
    context["groups"] = groups
    context["applied_at"] = applied_at
    return render(request, 'classification/gene_consensus_card.html', context)


def _comparisons(classification: Classification, groups: list[GeneConsensusGroup],
                 allele_origin_bucket: Optional[AlleleOriginBucket]) -> list[GeneContentComparison]:
    """ What each candidate would bring, beside what the record already has - only the empty ones are filled """
    e_keys = EvidenceKeyMap.cached()
    gene_scope_keys = ClassificationConsensus.gene_scope_keys(allele_origin_bucket)
    comparisons = []
    for group in groups:
        evidence = group.representative.published_evidence or {}
        values = []
        for key in gene_scope_keys:
            candidate = (evidence.get(key) or {}).get("value")
            if candidate is None or candidate == "" or candidate == []:
                continue
            values.append(GeneContentValue(label=e_keys.get(key).pretty_label,
                                           current=classification.get(key),
                                           candidate=candidate))
        comparisons.append(GeneContentComparison(group=group, values=values))
    return comparisons


def _apply_gene_consensus(request: HttpRequest, classification: Classification) -> HttpResponseBase:
    copy_from_id = request.POST.get("copy_gene_from_vcm_id")
    if not copy_from_id:
        raise Http404("No gene content record was chosen")
    copy_from = get_object_or_404(ClassificationModification, pk=int(copy_from_id))
    copy_from.check_can_view(request.user)
    ClassificationConsensus(modification=copy_from, copy_scopes=COPY_SCOPES_GENE).apply_to(classification,
                                                                                             request.user)
    return JsonResponse({"applied": True})
