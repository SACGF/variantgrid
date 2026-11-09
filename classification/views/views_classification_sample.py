"""
The classification form's Change sample dialog, opened from the summary card's Sample row: a lab head moves the
record to another sample (classification/sample_reassignment.py).

Entry points: classification_sample_panel (the dialog's preview with ?sample_id=, and the POST that reassigns),
reached at classification_sample.
"""
from django.core.exceptions import PermissionDenied
from django.http import Http404
from django.http.request import HttpRequest
from django.http.response import HttpResponseBase, JsonResponse
from django.shortcuts import get_object_or_404, render

from classification.models.classification import Classification
from classification.sample_reassignment import (
    check_can_reassign_sample,
    reassign_sample,
    sample_has_genotype,
    sample_reassignment_preview,
)
from snpdb.models.models_vcf import Sample


def classification_sample_panel(request: HttpRequest, classification_id: int) -> HttpResponseBase:
    classification = get_object_or_404(Classification, pk=classification_id)
    if not classification.can_write(request.user):
        raise PermissionDenied(f"You do not have WRITE permission on classification {classification_id}")

    if request.method == "POST":
        return _reassign(request, classification)

    sample_id = request.GET.get("sample_id")
    if not sample_id:
        raise Http404("No sample was chosen")
    sample = get_object_or_404(Sample, pk=int(sample_id))
    context = {"classification": classification}
    try:
        context["preview"] = sample_reassignment_preview(classification, sample, request.user)
    except PermissionDenied as pd:
        context["error"] = str(pd)
    return render(request, 'classification/sample_change_preview.html', context)


def _reassign(request: HttpRequest, classification: Classification) -> HttpResponseBase:
    sample_id = request.POST.get("sample_id")
    if not sample_id:
        raise Http404("No sample was chosen")
    sample = get_object_or_404(Sample, pk=int(sample_id))
    check_can_reassign_sample(classification, sample, request.user)
    if not request.POST.get("confirm_no_genotype") and not sample_has_genotype(classification, sample):
        return JsonResponse({"error": f"{sample} has no genotype call for this variant"}, status=400)
    repopulate = bool(request.POST.get("repopulate"))
    reassign_sample(classification, sample, request.user, repopulate=repopulate)
    return JsonResponse({"reassigned": True})
