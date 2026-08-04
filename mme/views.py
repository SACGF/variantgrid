from django.conf import settings
from django.contrib import messages
from django.core.exceptions import PermissionDenied
from django.http import HttpRequest, Http404
from django.shortcuts import render, get_object_or_404, redirect
from django.utils import timezone
from django.views.decorators.http import require_POST

from classification.enums.classification_enums import ShareLevel, ClinicalSignificance
from classification.models.classification import Classification, ClassificationModification
from mme.disclaimers import connected_nodes, effective_disclaimer, node_disclaimer
from mme.metrics import get_metrics
from mme.models import MMESubmission, MMESubmissionStatus, MMEMatchResult, MMEInboundMatch
from mme.serializers.patient_profile import (
    mme_eligible_classifications,
    classification_genomic_feature,
    classification_ontology_slots,
    build_patient_profile,
)
from mme.tasks import submit_mme_submission_task
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder

# Card states - which layer of the node/lab/record opt-in a classification fails.
MME_STATE_ELIGIBLE = "eligible"
MME_STATE_NOT_CANDIDATE = "not_candidate"
MME_STATE_NOT_SHARED = "not_shared"
MME_STATE_LAB_NOT_ENABLED = "lab_not_enabled"


def _external_patient_id(classification) -> str:
    """ Opaque, stable, non-PII id we send as MME patient.id. """
    return f"vg:{classification.pk}"


def check_can_view_classification(classification, user) -> None:
    """ Classification only carries a WRITE permission; visibility lives on its
        modifications, whose read permission comes from share level. Same rule as
        ClassificationRef.check_security. """
    if classification.can_write(user):
        return
    if not classification.latest_modification_for_user(user, exclude_withdrawn=False):
        raise PermissionDenied(f"You do not have READ permission to view {classification.pk}")


def _classification_is_eligible(classification) -> bool:
    return mme_eligible_classifications().filter(classification=classification).exists()


def _mme_state(classification) -> dict:
    """ Layer by layer rather than the single `eligible` boolean, so the card can say
        *why* a record does not qualify. """
    if _classification_is_eligible(classification):
        return {"state": MME_STATE_ELIGIBLE, "eligible": True, "clinical_significance": None}

    lab = getattr(classification, "lab", None)
    if not (lab and lab.mme_enabled):
        state = MME_STATE_LAB_NOT_ENABLED
        clinical_significance = None
    else:
        cm = (ClassificationModification.objects
              .filter(classification=classification, is_last_published=True)
              .order_by("-pk").first())
        share_level = cm.share_level if cm else None
        clinical_significance = cm.clinical_significance if cm else None
        if share_level != ShareLevel.PUBLIC.value or classification.withdrawn:
            state = MME_STATE_NOT_SHARED
        else:
            state = MME_STATE_NOT_CANDIDATE

    return {
        "state": state,
        "eligible": False,
        "clinical_significance": ClinicalSignificance.LABELS.get(clinical_significance),
    }


def mme_classification_panel(request: HttpRequest, classification_id: int):
    """ MatchMaker Exchange panel for a classification: shows the assembled profile
        (candidate gene/variant, curated disorders, HPO features), the eligibility state, a
        Submit button per configured node, and any inbound matches peers have made against
        this record. POSTing a node_id creates a draft. Renders in every state, including
        ineligible - it is embedded on the classification page. """
    classification = get_object_or_404(Classification, pk=classification_id)
    check_can_view_classification(classification, request.user)

    state = _mme_state(classification)
    eligible = state["eligible"]

    if request.method == "POST":
        if not eligible:
            raise Http404("Classification is not eligible for MatchMaker Exchange")
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
        "mme_state": state["state"],
        "clinical_significance": state["clinical_significance"],
        "mme_enabled": settings.MME_ENABLED,
        "nodes": sorted(settings.MME_NODES.keys()),
        "genomic_features": genomic_features,
        "features": features,
        "disorders": disorders,
        "submissions": MMESubmission.objects.filter(classification=classification).order_by("node_id"),
        "inbound_matches": (MMEInboundMatch.objects
                            .filter(classification=classification)
                            .select_related("inbound_query")
                            .order_by("-created")),
    }
    # ?fragment=1 skips the base template, for the AJAX-loaded card on the classification page
    template = ("mme/_mme_classification_card.html" if request.GET.get("fragment")
                else "mme/mme_classification_panel.html")
    return render(request, template, context)


def view_mme_submission(request: HttpRequest, submission_id: int):
    """ Show the exact profile that will be POSTed and let the curator confirm/submit. """
    submission = get_object_or_404(MMESubmission, pk=submission_id)
    check_can_view_classification(submission.classification, request.user)

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
        "disclaimer": effective_disclaimer(submission.node_id,
                                           submission.response_disclaimer,
                                           submission.response_terms),
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
        # build_patient_profile re-asserts eligibility - the record may have been downgraded
        # or unshared since the draft was created.
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


def view_mme_inbound_match(request: HttpRequest, inbound_match_id: int):
    """ Everything the querying peer gave us: their patient object, their contact, and the
        score. This is what the depositor notification links to. """
    match = get_object_or_404(MMEInboundMatch, pk=inbound_match_id)
    check_can_view_classification(match.classification, request.user)
    return render(request, "mme/mme_inbound_match.html", {
        "match": match,
        "peer_node_id": match.inbound_query.peer_node_id,
        "peer_disclaimer": node_disclaimer(match.inbound_query.peer_node_id),
    })


def mme_public_metrics(request: HttpRequest):
    """ MME requires metrics be published publicly - no auth, via PUBLIC_PATHS. Shares the
        get_metrics() cache with the API view. """
    return render(request, "mme/mme_metrics.html", {
        "metrics": get_metrics(),
        "disclaimer": settings.MME_DISCLAIMER,
        "terms": settings.MME_TERMS,
    })


def mme_public_disclaimers(request: HttpRequest):
    """ Each connected database's disclaimers and terms. Public, like the metrics page. """
    return render(request, "mme/mme_disclaimers.html", {"nodes": connected_nodes()})


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
