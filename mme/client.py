""" Outbound MatchMaker Exchange client.

Self-contained `requests` client that copies the *shape* of `library/oauth.py:ServerAuth`
(a small object wrapping `requests.post` + an auth header) without importing it or any
`sync/` code. MME auth is a single static per-node header, so there is no OAuth machinery.
"""
import requests
from django.conf import settings
from django.urls import reverse
from django.utils import timezone

from library.constants import MINUTE_SECS
from library.django_utils import get_url_from_view_path
from library.log_utils import AdminNotificationBuilder
from mme.models import MMESubmissionStatus, MMEMatchResult
from mme.notifications import notify_match
from mme.serializers.patient_profile import build_patient_profile


class MMEClient:
    """ Minimal MME node client. One instance per remote node. """

    def __init__(self, node_id: str):
        node = settings.MME_NODES[node_id]
        self.node_id = node_id
        self.base_url = node["base_url"].rstrip("/")
        self.token = node["token"]
        # v1.0 is the version every node in the federation speaks (only GeneMatcher also
        # offers v1.1), so it is the safe default when a node's config doesn't say.
        self.api_version = node.get("api_version", "1.0")

    @property
    def _headers(self) -> dict:
        content_type = f"application/vnd.ga4gh.matchmaker.v{self.api_version}+json"
        return {
            "Content-Type": content_type,
            "Accept": content_type,
            "X-Auth-Token": self.token,
        }

    def match(self, patient_profile: dict, timeout: int = MINUTE_SECS) -> dict:
        # disclaimer/terms are message-level (siblings of `patient`), telling the
        # receiving curator some phenotype terms may be ontology-derived (§6a).
        body = {"patient": patient_profile}
        if settings.MME_DISCLAIMER:
            body["disclaimer"] = settings.MME_DISCLAIMER
        if settings.MME_TERMS:
            body["terms"] = settings.MME_TERMS
        response = requests.post(
            url=f"{self.base_url}/match",
            json=body,
            headers=self._headers,
            timeout=timeout,
        )
        response.raise_for_status()
        return response.json()


def submit(submission) -> None:
    """ Build profile, POST, persist results. Notifies admins on failure
        (mirrors the Alissa upload failure-reporting pattern, no sync dep). """
    try:
        # Re-asserts eligibility - the record may have been downgraded, withdrawn or
        # unshared between the curator clicking submit and this worker running.
        profile = build_patient_profile(submission)
    except ValueError as ve:
        submission.status = MMESubmissionStatus.ERROR
        submission.error = str(ve)
        submission.save()
        raise
    submission.request_json = profile
    try:
        data = MMEClient(submission.node_id).match(profile)
    except Exception as e:
        submission.status = MMESubmissionStatus.ERROR
        submission.error = str(e)
        submission.save()
        nb = AdminNotificationBuilder("MME submission failed")
        nb.add_markdown(f"Classification {submission.classification_id} -> node "
                        f"`{submission.node_id}`: {e}")
        nb.send()
        raise

    submission.response_json = data
    submission.response_disclaimer = data.get("disclaimer") or ""
    submission.response_terms = data.get("terms") or ""
    submission.status = MMESubmissionStatus.SUBMITTED
    submission.submitted = timezone.now()
    submission.save()

    results = data.get("results", [])
    for result in results:
        patient = result.get("patient", {})
        contact = patient.get("contact", {})
        MMEMatchResult.objects.create(
            submission=submission,
            score=result.get("score", {}).get("patient", 0.0),
            matched_patient_id=patient.get("id", ""),
            contact_name=contact.get("name"),
            contact_href=contact.get("href"),
            patient_json=patient,
        )

    if results:
        _notify_curator_of_matches(submission, len(results))


def _notify_curator_of_matches(submission, num_results: int) -> None:
    """ Tell the curator and lab their submission returned candidate matches, via the same
        helper the inbound path uses. """
    url = get_url_from_view_path(
        reverse("mme_view_submission", kwargs={"submission_id": submission.pk}))
    plural = "es" if num_results != 1 else ""
    notify_match(
        classification=submission.classification,
        subject=f"MatchMaker Exchange: {num_results} possible match{plural}",
        body=(f"Your MME submission of classification #{submission.classification_id} "
              f"to node `{submission.node_id}` returned {num_results} possible "
              f"match{plural}. [View matches]({url})."),
    )
