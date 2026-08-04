from drf_spectacular.utils import extend_schema
from rest_framework.exceptions import ParseError
from rest_framework.parsers import JSONParser
from rest_framework.renderers import JSONRenderer
from rest_framework.response import Response
from rest_framework.views import APIView

from mme.auth import authenticate_mme_peer
from mme.disclaimers import mme_response_body
from mme.matching import find_matches
from mme.metrics import get_metrics
from mme.models import MMEInboundQuery
from mme.notifications import record_inbound_matches
from mme.tasks import notify_mme_depositors_task

MME_MEDIA_TYPE = "application/vnd.ga4gh.matchmaker.v1.1+json"


class MMEJSONParser(JSONParser):
    """ MME peers send the GA4GH vendor content-type; parse it as JSON. """
    media_type = MME_MEDIA_TYPE


class MMEJSONRenderer(JSONRenderer):
    """ Respond with the GA4GH vendor content-type. """
    media_type = MME_MEDIA_TYPE


class MMEMatchView(APIView):
    """ Inbound MatchMaker Exchange /match endpoint. Authenticated by the per-peer
        X-Auth-Token we issued that node (mme/auth.py), NOT by a VariantGrid user session
        (its path is exempted from GlobalLoginRequiredMiddleware via PUBLIC_PATHS). """
    authentication_classes = []          # token-header auth, not session/DRF user
    permission_classes = []
    # Accept both the GA4GH vendor media type and plain application/json.
    parser_classes = [MMEJSONParser, JSONParser]
    renderer_classes = [MMEJSONRenderer, JSONRenderer]

    @extend_schema(exclude=True)
    def post(self, request, *args, **kwargs):
        peer_node_id = authenticate_mme_peer(request)

        patient = (request.data or {}).get("patient")
        if not patient or not (patient.get("features") or patient.get("genomicFeatures")):
            raise ParseError("patient with features or genomicFeatures required")

        matches = find_matches(patient)          # simple similarity scoring (mme/matching.py)
        inbound_query = MMEInboundQuery.objects.create(
            peer_node_id=peer_node_id,
            request_json=request.data,
            num_results=len(matches),
        )
        record_inbound_matches(inbound_query, matches, patient)
        # Queued so the peer isn't kept waiting while we send email. Certification traffic
        # carries `"test": true` - audit those, but don't page real curators.
        if matches and not patient.get("test"):
            notify_mme_depositors_task.si(inbound_query.pk).apply_async()

        return Response(mme_response_body({"results": [m.as_result() for m in matches]}))


class MMEMetricsView(APIView):
    """ MME Metrics API - peer-token authenticated per the spec; the same numbers are
        published without auth at /mme/metrics. Caching lives in get_metrics() rather than
        @cache_page, which would serve responses without running the token check. """
    authentication_classes = []
    permission_classes = []
    renderer_classes = [MMEJSONRenderer, JSONRenderer]

    @extend_schema(exclude=True)
    def get(self, request, *args, **kwargs):
        authenticate_mme_peer(request)
        return Response(mme_response_body({"metrics": get_metrics()}))
