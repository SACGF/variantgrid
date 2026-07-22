from django.conf import settings
from drf_spectacular.utils import extend_schema
from rest_framework.exceptions import AuthenticationFailed, ParseError
from rest_framework.parsers import JSONParser
from rest_framework.renderers import JSONRenderer
from rest_framework.response import Response
from rest_framework.views import APIView

from mme.matching import find_matches
from mme.models import MMEInboundQuery

MME_MEDIA_TYPE = "application/vnd.ga4gh.matchmaker.v1.1+json"


class MMEJSONParser(JSONParser):
    """ MME peers send the GA4GH vendor content-type; parse it as JSON. """
    media_type = MME_MEDIA_TYPE


class MMEJSONRenderer(JSONRenderer):
    """ Respond with the GA4GH vendor content-type. """
    media_type = MME_MEDIA_TYPE


class MMEMatchView(APIView):
    """ Inbound MatchMaker Exchange /match endpoint. Authenticated by the shared
        X-Auth-Token we issue to peer nodes, NOT by a VariantGrid user session
        (its path is exempted from GlobalLoginRequiredMiddleware via PUBLIC_PATHS). """
    authentication_classes = []          # token-header auth, not session/DRF user
    permission_classes = []
    # Accept both the GA4GH vendor media type and plain application/json.
    parser_classes = [MMEJSONParser, JSONParser]
    renderer_classes = [MMEJSONRenderer, JSONRenderer]

    @extend_schema(exclude=True)
    def post(self, request, *args, **kwargs):
        if not settings.MME_ENABLED:
            raise AuthenticationFailed("MatchMaker Exchange is not enabled")

        token = request.headers.get("X-Auth-Token")
        if not token or not settings.MME_INBOUND_TOKEN or token != settings.MME_INBOUND_TOKEN:
            raise AuthenticationFailed("Invalid MME auth token")

        patient = (request.data or {}).get("patient")
        if not patient or not (patient.get("features") or patient.get("genomicFeatures")):
            raise ParseError("patient with features or genomicFeatures required")

        results = find_matches(patient)          # simple similarity scoring (mme/matching.py)
        MMEInboundQuery.objects.create(
            request_json=request.data,
            num_results=len(results),
        )
        return Response({"results": results})
