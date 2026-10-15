"""
The gene content box on the classification form - for records the create-from-variant page never saw:
made from the Classify & Report tab with nothing to copy, imported, or created before the gene was curated.

It offers the same deduplicated gene candidates as the create page
(classification/models/classification.py:ClassificationConsensus.gene_consensus_groups), and applying one
patches the empty fields as SubmissionSource.CONSENSUS.

Entry points: gene_consensus_panel (the card, its dialog rows, and the apply POST), reached at
classification_gene_consensus.
"""
from dataclasses import replace
from datetime import date, datetime
from typing import Optional, Union

from django.core.exceptions import PermissionDenied
from django.http import Http404
from django.http.request import HttpRequest
from django.http.response import HttpResponseBase, JsonResponse
from django.shortcuts import get_object_or_404, render

from classification.enums import AlleleOriginBucket, SubmissionSource
from classification.models import (
    Classification,
    ClassificationConsensus,
    ClassificationModification,
    ConsensusGroup,
)
from classification.models.classification import COPY_SCOPES_GENE
from classification.models.evidence_key import EvidenceKeyMap
from genes.models.models_gene import AmbiguousGeneSymbolAlias


def _as_date(value: Union[date, datetime]) -> date:
    return value.date() if isinstance(value, datetime) else value


def classification_gene_symbol(classification: Classification) -> Optional[str]:
    """ The gene the record is about, by its current symbol (ClassificationModification.gene_symbol_best), so a
        record submitted with an old name or alias finds the gene's other records. None for an ambiguous alias """
    try:
        if gene_symbol := classification.last_edited_version.gene_symbol_best():
            return gene_symbol.symbol
    except AmbiguousGeneSymbolAlias:
        pass
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


def newer_than(groups: list[ConsensusGroup], applied_at: datetime) -> list[ConsensusGroup]:
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
        context["groups"] = _with_current_values(classification, groups)
        return render(request, 'classification/gene_consensus_rows.html', context)

    # The box only takes space while it is useful - once gene content has been copied in it collapses to a
    # line, and only while the gene has been curated again since
    applied_at = gene_content_applied_at(classification, allele_origin_bucket)
    if applied_at:
        groups = newer_than(groups, applied_at)
    context["groups"] = groups
    context["applied_at"] = applied_at
    return render(request, 'classification/gene_consensus_card.html', context)


def _with_current_values(classification: Classification, groups: list[ConsensusGroup]) -> list[ConsensusGroup]:
    """ The same rows the create page shows, each value beside what this record already has - only the empty
        ones are filled """
    e_keys = EvidenceKeyMap.cached()

    def current(key: str):
        value = classification.get(key)
        return e_keys.get(key).pretty_value(value) if value is not None else None

    return [replace(group, content=[replace(value, current=current(value.key)) for value in group.content])
            for group in groups]


def _apply_gene_consensus(request: HttpRequest, classification: Classification) -> HttpResponseBase:
    copy_from_id = request.POST.get("copy_gene_from_vcm_id")
    if not copy_from_id:
        raise Http404("No gene content record was chosen")
    copy_from = get_object_or_404(ClassificationModification, pk=int(copy_from_id))
    copy_from.check_can_view(request.user)
    ClassificationConsensus(modification=copy_from, copy_scopes=COPY_SCOPES_GENE).apply_to(classification,
                                                                                             request.user)
    return JsonResponse({"applied": True})
