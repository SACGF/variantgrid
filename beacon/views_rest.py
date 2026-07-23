""" Inbound Beacon v2 REST endpoints (§6): framework (identity/config/map) + the
    g_variants model endpoint. Plain DRF APIViews (no routers), following snpdb/views_rest.

Anonymous requests are allowed (permission AllowAny + the /beacon/ PUBLIC_PATHS entry);
per-tier data scope is still enforced by filter_for_user (anonymous -> public group).
"""
import logging

from django.conf import settings
from django.urls import reverse
from drf_spectacular.utils import extend_schema
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


def _require_enabled():
    if not settings.BEACON_ENABLED:
        raise NotFound("Beacon is not enabled")


class _BeaconFrameworkView(APIView):
    """ Base for the static framework endpoints - anonymous-readable, GET only. """
    permission_classes = []

    def _payload(self, request) -> dict:
        raise NotImplementedError()

    @extend_schema(exclude=True)
    def get(self, request, *args, **kwargs):
        _require_enabled()
        return Response(info_response(self._payload(request)))


class BeaconInfoView(_BeaconFrameworkView):
    def _payload(self, request):
        return schema.beacon_info()


class BeaconConfigurationView(_BeaconFrameworkView):
    def _payload(self, request):
        return schema.configuration()


class BeaconEntryTypesView(_BeaconFrameworkView):
    def _payload(self, request):
        return schema.entry_types()


class BeaconFilteringTermsView(_BeaconFrameworkView):
    def _payload(self, request):
        return schema.filtering_terms()


class BeaconMapView(_BeaconFrameworkView):
    def _payload(self, request):
        g_variants_url = request.build_absolute_uri(reverse("beacon_g_variants"))
        return schema.endpoint_map(g_variants_url)


class BeaconServiceInfoView(APIView):
    """ GA4GH service-info profile - returned raw (not wrapped in the Beacon meta). """
    permission_classes = []

    @extend_schema(exclude=True)
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


def _received_request_summary(params: dict, requested_granularity: str) -> dict:
    return {
        "requestedGranularity": requested_granularity,
        "requestParameters": {k: params.get(k) for k in _PARAM_KEYS},
    }


class BeaconGVariantsView(APIView):
    """ Beacon v2 g_variants query. Exact-coordinate lookup over two datasets
        (observations + public classifications), tiered boolean/count/record by the
        requester's permission. GET (flat params) and POST (Beacon request envelope). """
    permission_classes = []

    @extend_schema(exclude=True)
    def get(self, request, *args, **kwargs):
        return self._handle(request)

    @extend_schema(exclude=True)
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
                return Response(self._empty_response(params, requested_granularity, authenticated))
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
            _received_request_summary(params, requested_granularity),
            exists=exists,
            num_total_results=num_total_results,
            result_sets=result_sets,
        )

        self._audit(raw, granularity, authenticated, obs, cls)
        return Response(response)

    @staticmethod
    def _empty_response(params, requested_granularity, authenticated):
        """ Valid, empty g_variants envelope for a bare (no-coordinate) query. """
        granularity = clamp_granularity(requested_granularity, authenticated)
        return query_response(
            granularity,
            _received_request_summary(params, requested_granularity),
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

    @extend_schema(exclude=True)
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
        summary = {"requestedGranularity": RECORD, "variantInternalId": str(variant_id)}
        return Response(query_response(
            granularity, summary, exists=exists, num_total_results=num_total_results,
            result_sets=[obs.result_set(), cls.result_set()]))
