from django.conf import settings
from drf_spectacular.utils import extend_schema
from rest_framework import status
from rest_framework.exceptions import ParseError
from rest_framework.parsers import JSONParser
from rest_framework.renderers import JSONRenderer
from rest_framework.response import Response
from rest_framework.views import APIView

from mme.auth import authenticate_mme_peer, authenticate_mme_status_request
from mme.disclaimers import mme_response_body
from mme.matching import find_matches
from mme.metrics import get_metrics
from mme.models import MMEInboundQuery
from mme.notifications import record_inbound_matches
from mme.tasks import notify_mme_depositors_task
from mme.versioning import (
    MME_MEDIA_TYPE,
    MMEContentNegotiation,
    MMEJSONParser,
    MMEJSONRenderer,
    check_request_version,
    software_version,
    supported_media_types,
)


class MMEAPIView(APIView):
    """ Shared config for the inbound MME endpoints. Authenticated by the per-peer
        X-Auth-Token we issued that node (mme/auth.py), NOT by a VariantGrid user session
        (their paths are exempted from GlobalLoginRequiredMiddleware via PUBLIC_PATHS). """
    authentication_classes = []          # token-header auth, not session/DRF user
    permission_classes = []
    content_negotiation_class = MMEContentNegotiation

    def get_authenticate_header(self, request):
        """ MME specifies 401 for a bad API key. DRF downgrades AuthenticationFailed to 403
            unless the view can name a WWW-Authenticate scheme, and listing no
            authentication_classes leaves it with nothing to name - so name ours. """
        return "X-Auth-Token"

    def initial(self, request, *args, **kwargs):
        super().initial(request, *args, **kwargs)
        check_request_version(request)

    def handle_exception(self, exc):
        """ MME requires error bodies carry a human-readable `message`; DRF's default is
            `detail`. """
        response = super().handle_exception(exc)
        if isinstance(response.data, dict) and "detail" in response.data:
            response.data = {"message": str(response.data["detail"])}
        return response

    def finalize_response(self, request, response, *args, **kwargs):
        response = super().finalize_response(request, response, *args, **kwargs)
        if response.status_code == status.HTTP_406_NOT_ACCEPTABLE:
            # Tell the peer which version we do support, so they can retry. The header
            # comes off the renderer, so that is what has to carry it.
            response.accepted_renderer = MMEJSONRenderer()
            response.accepted_media_type = MME_MEDIA_TYPE
        return response


class MMEMatchView(MMEAPIView):
    """ Inbound MatchMaker Exchange /match endpoint. """
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


class MMEMetricsView(MMEAPIView):
    """ MME Metrics API - peer-token authenticated per the spec once MME is on (see
        authenticate_mme_status_request); the same numbers are published without auth at
        /mme/metrics. Caching lives in get_metrics() rather than @cache_page, which would
        serve responses without running the token check. """
    # metrics-api.md specifies application/json, unlike /match's vendor media type - so
    # plain JSON is what a peer sending no Accept header gets.
    renderer_classes = [JSONRenderer, MMEJSONRenderer]

    @extend_schema(exclude=True)
    def get(self, request, *args, **kwargs):
        authenticate_mme_status_request(request)
        return Response(mme_response_body({"metrics": get_metrics()}))


class MMEHeartbeatView(MMEAPIView):
    """ MME Heartbeat API - https://github.com/ga4gh/mme-apis/blob/1.1.1/heartbeat-api.md

        Liveness plus the machine-readable list of protocol versions we serve, so a peer
        can negotiate without trial and error. New in v1.1. """
    # heartbeat-api.md specifies application/json, as for metrics.
    renderer_classes = [JSONRenderer, MMEJSONRenderer]

    @extend_schema(exclude=True)
    def get(self, request, *args, **kwargs):
        authenticate_mme_status_request(request)
        return Response(mme_response_body({
            "heartbeat": {
                "production": not settings.MME_TEST,
                "version": software_version(),
                "accept": supported_media_types(),
            },
        }))
