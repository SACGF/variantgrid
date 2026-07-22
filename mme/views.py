from django.conf import settings
from django.contrib import messages
from django.http import HttpRequest, Http404
from django.shortcuts import render, get_object_or_404, redirect
from django.utils import timezone
from django.views.decorators.http import require_POST

from classification.models.classification import Classification
from mme.models import MMESubmission, MMESubmissionStatus, MMEMatchResult
from mme.serializers.patient_profile import (
    mme_eligible_classifications,
    classification_genomic_feature,
    classification_ontology_slots,
    build_patient_profile,
)
from mme.tasks import submit_mme_submission_task
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


def _external_patient_id(classification) -> str:
    """ Opaque, stable, non-PII id we send as MME patient.id. """
    return f"vg:{classification.pk}"


def _classification_is_eligible(classification) -> bool:
    return mme_eligible_classifications().filter(classification=classification).exists()


def mme_classification_panel(request: HttpRequest, classification_id: int):
    """ MatchMaker Exchange panel for a classification: shows the assembled profile
        (candidate gene/variant, curated disorders, HPO features) and a Submit button
        per configured node. POSTing a node_id creates a draft submission. """
    classification = Classification.get_for_user(request.user, classification_id)
    classification.check_can_view(request.user)

    eligible = _classification_is_eligible(classification)

    if request.method == "POST":
        if not eligible:
            raise Http404("Classification is not shared at ShareLevel.PUBLIC")
        node_id = request.POST.get("node_id")
        if node_id not in settings.MME_NODES:
            raise Http404(f"Unknown MME node '{node_id}'")
        submission, _ = MMESubmission.objects.get_or_create(
            classification=classification,
            node_id=node_id,
            defaults={"external_patient_id": _external_patient_id(classification)},
        )
        return redirect("mme_view_submission", submission_id=submission.pk)

    genomic_features = classification_genomic_feature(classification)
    features, disorders = classification_ontology_slots(classification)

    context = {
        "classification": classification,
        "eligible": eligible,
        "mme_enabled": settings.MME_ENABLED,
        "nodes": sorted(settings.MME_NODES.keys()),
        "genomic_features": genomic_features,
        "features": features,
        "disorders": disorders,
        "submissions": MMESubmission.objects.filter(classification=classification).order_by("node_id"),
    }
    return render(request, "mme/mme_classification_panel.html", context)


def view_mme_submission(request: HttpRequest, submission_id: int):
    """ Show the exact profile that will be POSTed and let the curator confirm/submit. """
    submission = get_object_or_404(MMESubmission, pk=submission_id)
    submission.classification.check_can_view(request.user)

    profile = submission.request_json
    profile_error = None
    if profile is None:
        try:
            profile = build_patient_profile(submission)
        except ValueError as ve:
            profile_error = str(ve)

    context = {
        "submission": submission,
        "profile": profile,
        "profile_error": profile_error,
        "mme_enabled": settings.MME_ENABLED,
        "results": MMEMatchResult.objects.filter(submission=submission).order_by("-score"),
    }
    return render(request, "mme/mme_submission.html", context)


@require_POST
def submit_mme_submission(request: HttpRequest, submission_id: int):
    """ Persist the confirmed profile and fire the worker task. """
    submission = get_object_or_404(MMESubmission, pk=submission_id)
    submission.classification.check_can_write(request.user)

    if not settings.MME_ENABLED:
        messages.error(request, "MatchMaker Exchange is not enabled")
        return redirect("mme_view_submission", submission_id=submission.pk)

    try:
        submission.request_json = build_patient_profile(submission)
    except ValueError as ve:
        messages.error(request, f"Cannot submit: {ve}")
        return redirect("mme_view_submission", submission_id=submission.pk)

    submission.status = MMESubmissionStatus.DRAFT
    submission.error = None
    submission.submitted_by = request.user
    submission.submitted = timezone.now()
    submission.save()

    submit_mme_submission_task.si(submission.pk).apply_async()
    messages.info(request, f"Submission to '{submission.node_id}' queued.")
    return redirect("mme_view_submission", submission_id=submission.pk)


class MMEMatchResultColumns(DatatableConfig[MMEMatchResult]):
    """ Results grid for one submission's returned candidate patients. """

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn('score', orderable=True, default_sort=SortOrder.DESC),
            RichColumn('matched_patient_id', label='Matched patient', orderable=True),
            RichColumn('contact_name', label='Contact', orderable=True),
            RichColumn('contact_href', label='Contact link', orderable=True),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True),
        ]

    def get_initial_queryset(self):
        submission_id = self.get_query_param("submission_id")
        qs = MMEMatchResult.objects.all()
        if submission_id:
            qs = qs.filter(submission_id=submission_id,
                           submission__classification__in=Classification.filter_for_user(self.user))
        else:
            qs = qs.none()
        return qs
