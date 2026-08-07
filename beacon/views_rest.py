""" Inbound Beacon v2 REST endpoints (§6): framework (identity/config/map) + the
    g_variants model endpoint. Plain DRF APIViews (no routers), following snpdb/views_rest.

Anonymous requests are allowed (permission AllowAny + the /beacon/ PUBLIC_PATHS entry);
per-tier data scope is still enforced by filter_for_user (anonymous -> public group).
"""
import logging

from django.conf import settings
from django.urls import reverse
from drf_spectacular.types import OpenApiTypes
from drf_spectacular.utils import extend_schema, extend_schema_view, OpenApiParameter
from rest_framework.exceptions import NotFound, ParseError
from rest_framework.response import Response
from rest_framework.views import APIView

from beacon import schema
from beacon.datasets import observations_dataset, classifications_dataset
from beacon.models import BeaconInboundQuery
from beacon.response import clamp_granularity, query_response, info_response, RECORD
from beacon.variant_mapping import (
    BeaconQueryError,
    genome_build_from_assembly,
    variant_qs_for_beacon_query,
)
from snpdb.models import Variant

# Groups every inbound Beacon v2 endpoint under one heading in the Swagger UI.
BEACON_TAGS = ["Beacon v2"]


def _require_enabled():
    if not settings.BEACON_ENABLED:
        raise NotFound("Beacon is not enabled")


class _BeaconFrameworkView(APIView):
    """ Base for the static framework endpoints - anonymous-readable, GET only. """
    permission_classes = []

    def _payload(self, request) -> dict:
        raise NotImplementedError()

    def get(self, request, *args, **kwargs):
        _require_enabled()
        return Response(info_response(self._payload(request)))


@extend_schema_view(get=extend_schema(
    tags=BEACON_TAGS,
    summary="Beacon information (identity & organisation)",
    description="GA4GH Beacon v2 identity/metadata. Anonymous-readable.",
    responses=OpenApiTypes.OBJECT,
))
class BeaconInfoView(_BeaconFrameworkView):
    def _payload(self, request):
        return schema.beacon_info()


@extend_schema_view(get=extend_schema(
    tags=BEACON_TAGS,
    summary="Beacon configuration",
    description="Beacon v2 configuration: entry types, security attributes and maturity.",
    responses=OpenApiTypes.OBJECT,
))
class BeaconConfigurationView(_BeaconFrameworkView):
    def _payload(self, request):
        return schema.configuration()


@extend_schema_view(get=extend_schema(
    tags=BEACON_TAGS,
    summary="Beacon entry types (data model schemas)",
    responses=OpenApiTypes.OBJECT,
))
class BeaconEntryTypesView(_BeaconFrameworkView):
    def _payload(self, request):
        return schema.entry_types()


@extend_schema_view(get=extend_schema(
    tags=BEACON_TAGS,
    summary="Beacon filtering terms",
    responses=OpenApiTypes.OBJECT,
))
class BeaconFilteringTermsView(_BeaconFrameworkView):
    def _payload(self, request):
        return schema.filtering_terms()


@extend_schema_view(get=extend_schema(
    tags=BEACON_TAGS,
    summary="Beacon endpoint map (endpoint discovery)",
    responses=OpenApiTypes.OBJECT,
))
class BeaconMapView(_BeaconFrameworkView):
    def _payload(self, request):
        g_variants_url = request.build_absolute_uri(reverse("beacon_g_variants"))
        return schema.endpoint_map(g_variants_url)


class BeaconServiceInfoView(APIView):
    """ GA4GH service-info profile - returned raw (not wrapped in the Beacon meta). """
    permission_classes = []

    @extend_schema(
        tags=BEACON_TAGS,
        summary="GA4GH service-info",
        description="GA4GH service-info profile, returned raw (not wrapped in the Beacon meta envelope).",
        responses=OpenApiTypes.OBJECT,
    )
    def get(self, request, *args, **kwargs):
        _require_enabled()
        return Response(schema.service_info())


# ------------------------------ g_variants ------------------------------

_PARAM_KEYS = ("referenceName", "start", "referenceBases", "alternateBases", "assemblyId")


def _extract_params(request) -> tuple[dict, str, dict]:
    """ Merge GET query params / POST body into a flat request-parameter dict.
        Returns (params, requested_granularity, raw_request_json_for_audit). """
    if request.method == "POST":
        body = request.data or {}
        query = body.get("query", {}) if isinstance(body, dict) else {}
        params = dict(query.get("requestParameters", {}))
        # tolerate flat top-level params too
        for k in _PARAM_KEYS:
            if k not in params and k in body:
                params[k] = body[k]
        requested = query.get("requestedGranularity") or body.get("requestedGranularity")
        raw = body
    else:
        qp = request.query_params
        params = {k: qp.get(k) for k in _PARAM_KEYS if qp.get(k) is not None}
        requested = qp.get("requestedGranularity")
        raw = dict(qp)
    return params, requested, raw


def _received_request_summary(requested_granularity: str, extra: dict = None) -> dict:
    """ Spec-complete beaconReceivedRequestSummary: apiVersion, requestedSchemas,
        pagination and requestedGranularity are all required.

        `requestParameters` is deliberately omitted: the framework schema constrains its
        values to objects (additionalProperties: {type: object}), so echoing flat g_variants
        coordinates (referenceName/start/... as scalars) would fail validation. It is an
        optional field; the exact request is still captured in the BeaconInboundQuery audit row. """
    config = settings.BEACON_CONFIG
    summary = {
        "apiVersion": config["api_version"],
        "requestedSchemas": [],
        "pagination": {"skip": 0, "limit": 0},
        "requestedGranularity": requested_granularity or config.get("default_granularity", "boolean"),
    }
    if extra:
        summary.update(extra)
    return summary


class BeaconGVariantsView(APIView):
    """ Beacon v2 g_variants query. Exact-coordinate lookup over two datasets
        (observations + public classifications), tiered boolean/count/record by the
        requester's permission. GET (flat params) and POST (Beacon request envelope). """
    permission_classes = []

    @extend_schema(
        tags=BEACON_TAGS,
        summary="Query genomic variants (Beacon v2 g_variants)",
        description="Exact-coordinate lookup over observations + public classifications. "
                    "Anonymous requests are allowed; the response granularity (boolean/count/record) "
                    "is clamped to the requester's permission tier.",
        parameters=[
            OpenApiParameter("referenceName", OpenApiTypes.STR, OpenApiParameter.QUERY,
                             description="Chromosome / reference sequence, eg '17'"),
            OpenApiParameter("start", OpenApiTypes.INT, OpenApiParameter.QUERY,
                             description="0-based start position"),
            OpenApiParameter("referenceBases", OpenApiTypes.STR, OpenApiParameter.QUERY,
                             description="Reference allele bases"),
            OpenApiParameter("alternateBases", OpenApiTypes.STR, OpenApiParameter.QUERY,
                             description="Alternate allele bases"),
            OpenApiParameter("assemblyId", OpenApiTypes.STR, OpenApiParameter.QUERY,
                             description="Genome assembly, eg 'GRCh38'"),
            OpenApiParameter("requestedGranularity", OpenApiTypes.STR, OpenApiParameter.QUERY,
                             description="Requested response detail: 'boolean', 'count' or 'record'"),
        ],
        responses=OpenApiTypes.OBJECT,
    )
    def get(self, request, *args, **kwargs):
        return self._handle(request)

    @extend_schema(
        tags=BEACON_TAGS,
        summary="Query genomic variants (Beacon v2 request envelope)",
        description="POST form of the g_variants query. The body is a Beacon request envelope "
                    "carrying the coordinate fields under query.requestParameters (flat top-level "
                    "params are also tolerated).",
        request=OpenApiTypes.OBJECT,
        responses=OpenApiTypes.OBJECT,
    )
    def post(self, request, *args, **kwargs):
        return self._handle(request)

    def _handle(self, request):
        _require_enabled()
        params, requested_granularity, raw = _extract_params(request)

        authenticated = bool(request.user and request.user.is_authenticated)

        missing = [k for k in _PARAM_KEYS if not params.get(k)]
        if missing:
            # No coordinate params at all => a bare "list" query: return a valid, empty
            # result envelope (200) so spec clients / the beacon-verifier get a conformant
            # response. A partial coordinate is genuinely malformed => 400.
            if len(missing) == len(_PARAM_KEYS):
                return Response(self._empty_response(requested_granularity, authenticated))
            raise ParseError(f"Missing required Beacon g_variants parameters: {', '.join(missing)}")

        try:
            genome_build = genome_build_from_assembly(params["assemblyId"])
            start = int(params["start"])
        except BeaconQueryError as e:
            raise ParseError(str(e)) from e
        except (TypeError, ValueError) as e:
            raise ParseError(f"Invalid 'start' coordinate: {params.get('start')!r}") from e

        variant = variant_qs_for_beacon_query(
            params["referenceName"], start, params["referenceBases"],
            params["alternateBases"], genome_build).first()
        allele = variant.allele if variant is not None else None

        granularity = clamp_granularity(requested_granularity, authenticated)
        is_record = granularity == RECORD

        obs = observations_dataset(request.user, variant, genome_build)
        cls = classifications_dataset(request.user, allele, is_record)

        exists = obs.exists or cls.exists
        num_total_results = (obs.reportable_count or 0) + (cls.reportable_count or 0)
        result_sets = [obs.result_set(), cls.result_set()]

        response = query_response(
            granularity,
            _received_request_summary(requested_granularity),
            exists=exists,
            num_total_results=num_total_results,
            result_sets=result_sets,
        )

        self._audit(raw, granularity, authenticated, obs, cls)
        return Response(response)

    @staticmethod
    def _empty_response(requested_granularity, authenticated):
        """ Valid, empty g_variants envelope for a bare (no-coordinate) query. """
        granularity = clamp_granularity(requested_granularity, authenticated)
        return query_response(
            granularity,
            _received_request_summary(requested_granularity),
            exists=False, num_total_results=0, result_sets=[])

    @staticmethod
    def _audit(request_json, granularity, authenticated, obs, cls):
        try:
            BeaconInboundQuery.objects.create(
                request_json=request_json,
                granularity=granularity,
                authenticated=authenticated,
                observations_exists=obs.exists,
                observations_count=obs.count,
                classifications_exists=cls.exists,
                classifications_count=cls.count,
            )
        except Exception:  # audit must never break the response
            logging.exception("Beacon: failed to write BeaconInboundQuery audit row")
        logging.info("Beacon g_variants: granularity=%s authenticated=%s observations(exists=%s,count=%s) "
                     "classifications(exists=%s,count=%s)", granularity, authenticated,
                     obs.exists, obs.count, cls.exists, cls.count)


class BeaconGVariantByIdView(APIView):
    """ GET /g_variants/{id} - single variant detail by our variantInternalId (Variant pk).
        Runs the same two-dataset lookup as the query, clamped by the requester's tier. """
    permission_classes = []

    @extend_schema(
        tags=BEACON_TAGS,
        summary="Retrieve a single genomic variant by internal id",
        description="Single-variant detail by variantInternalId (the VariantGrid Variant pk). "
                    "Runs the same two-dataset lookup as the query, clamped to the requester's tier.",
        parameters=[
            OpenApiParameter("variant_id", OpenApiTypes.INT, OpenApiParameter.PATH,
                             description="variantInternalId (VariantGrid Variant pk)"),
        ],
        responses=OpenApiTypes.OBJECT,
    )
    def get(self, request, variant_id, *args, **kwargs):
        _require_enabled()
        variant = Variant.objects.filter(pk=variant_id).first()
        if variant is None:
            raise NotFound(f"No variant with internal id {variant_id}")
        genome_build = variant.any_genome_build

        authenticated = bool(request.user and request.user.is_authenticated)
        granularity = clamp_granularity(RECORD, authenticated)
        is_record = granularity == RECORD

        obs = observations_dataset(request.user, variant, genome_build)
        cls = classifications_dataset(request.user, variant.allele, is_record)

        exists = obs.exists or cls.exists
        num_total_results = (obs.reportable_count or 0) + (cls.reportable_count or 0)
        summary = _received_request_summary(RECORD, extra={"variantInternalId": str(variant_id)})
        return Response(query_response(
            granularity, summary, exists=exists, num_total_results=num_total_results,
            result_sets=[obs.result_set(), cls.result_set()]))
