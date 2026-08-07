""" MME API versioning - https://github.com/ga4gh/mme-apis/blob/master/search-api.md

The API version rides in the media type (`application/vnd.ga4gh.matchmaker.v1.1+json`),
not in the URL, so there is no path component to version. A request naming a *major*
version we serve is compatible - minor versions are backwards compatible - and the response
echoes back the version the peer asked for. A major version we cannot serve is 406 (not
415), carrying our latest supported version in the response Content-Type.

This matters in practice: most of the deployed federation (DECIPHER, Monarch, MyGene2,
matchbox, PhenomeCentral, PubCaseFinder) speaks v1.0, and only GeneMatcher speaks v1.1.
"""
import re

from rest_framework.exceptions import NotAcceptable
from rest_framework.negotiation import DefaultContentNegotiation
from rest_framework.parsers import JSONParser
from rest_framework.renderers import JSONRenderer

from manual.models import Deployment


def _media_type(version: str) -> str:
    return f"application/vnd.ga4gh.matchmaker.v{version}+json"


# The version of the spec we implement, and what we advertise as our latest.
MME_API_VERSION = "1.1"
MME_MEDIA_TYPE = _media_type(MME_API_VERSION)

# Versions we serve, newest last - published verbatim as the /heartbeat `accept` array.
# Peers upgrade at their own pace (v1.0 and v1.1 have coexisted in the federation for
# years), so serving several at one base URL is the norm, and is how a future v2 stays
# backwards compatible with v1 peers.
MME_SUPPORTED_VERSIONS = ("1.0", "1.1")

# Derived: the compatibility test is on the major, since minors are cross-compatible.
# Adding a major is not just a matter of listing it - majors are by definition not
# wire-compatible, so a new one also needs its own profile builder + matching semantics,
# selected on `request_version(request)`. The list is the switch, not the work.
MME_SUPPORTED_MAJOR_VERSIONS = frozenset(v.split(".")[0] for v in MME_SUPPORTED_VERSIONS)


def supported_media_types() -> list[str]:
    """ Every version we accept, as media types - the /heartbeat `accept` array. """
    return [_media_type(v) for v in MME_SUPPORTED_VERSIONS]


def software_version() -> str:
    """ OUR release version for /heartbeat `version` - the running build, not the API
        version. VariantGrid has no semantic release number, so the last recorded
        deployment's git hash is what identifies a build (GeneMatcher returns e.g. "6.2"). """
    deployment = Deployment.objects.order_by("-created").first()
    if deployment and deployment.git_hash:
        return deployment.git_hash[:8]
    return "unknown"

_MME_MEDIA_TYPE_RE = re.compile(
    r"application/vnd\.ga4gh\.matchmaker\.v(?P<major>\d+)\.(?P<minor>\d+)\+json")


class MMEJSONParser(JSONParser):
    """ MME peers send the GA4GH vendor content-type; parse it as JSON. Any v1.x lands
        here via MMEContentNegotiation, regardless of this declared media_type. """
    media_type = MME_MEDIA_TYPE


class MMEJSONRenderer(JSONRenderer):
    """ Respond with the GA4GH vendor content-type. """
    media_type = MME_MEDIA_TYPE


def parse_media_type_version(value: str | None) -> tuple[str, str] | None:
    """ (major, minor) if `value` is an MME vendor media type, else None. """
    media_type = (value or "").split(";")[0].strip()
    if m := _MME_MEDIA_TYPE_RE.fullmatch(media_type):
        return m.group("major"), m.group("minor")
    return None


def compatible_media_type(value: str | None) -> str | None:
    """ The bare media type if `value` names an MME version we can serve, else None. """
    if version := parse_media_type_version(value):
        if version[0] in MME_SUPPORTED_MAJOR_VERSIONS:
            return (value or "").split(";")[0].strip()
    return None


def request_version(request) -> tuple[str, str]:
    """ The (major, minor) a peer is speaking, defaulting to what we implement when they
        send plain JSON. The seam for serving a second major: branch on the major here
        rather than adding a versioned URL, so peers keep one base_url across upgrades. """
    version = parse_media_type_version(request.content_type)
    if version is None:
        version = parse_media_type_version(request.headers.get("Accept"))
    return version or tuple(MME_API_VERSION.split("."))


def check_request_version(request) -> None:
    """ A request carrying a body names its API version in Content-Type. A bodyless GET
        (/metrics) has no meaningful content type, and plain application/json is let
        through - see parser_classes - but a version we cannot serve is 406 rather than an
        unsupported-media-type 415. Accept-header versions are negotiated separately, by
        MMEContentNegotiation. """
    if request.method not in ("POST", "PUT", "PATCH"):
        return
    media_type = (request.content_type or "").split(";")[0].strip()
    if not media_type or media_type == "application/json":
        return
    if compatible_media_type(media_type):
        return
    raise NotAcceptable(f"Unsupported MatchMaker Exchange API version '{media_type}' - "
                        f"this server supports v{MME_API_VERSION}")


class MMEContentNegotiation(DefaultContentNegotiation):
    """ Route any compatible v1.x media type to the MME parser/renderer. DRF matches media
        types by exact sub-type, so a v1.0 peer would otherwise get 415/406 off our v1.1
        declarations. """

    def select_parser(self, request, parsers):
        if compatible_media_type(request.content_type):
            for parser in parsers:
                if isinstance(parser, MMEJSONParser):
                    return parser
        return super().select_parser(request, parsers)

    def select_renderer(self, request, renderers, format_suffix=None):
        for accepted in request.headers.get("Accept", "").split(","):
            # Echo their version back, per the spec - not necessarily the one we declare.
            if media_type := compatible_media_type(accepted):
                for renderer in renderers:
                    if isinstance(renderer, MMEJSONRenderer):
                        # Response.rendered_content takes the Content-Type from
                        # renderer.media_type rather than the negotiated type, so pin the
                        # peer's version onto the instance (one per request) to echo it.
                        renderer.media_type = media_type
                        return renderer, media_type
        return super().select_renderer(request, renderers, format_suffix)
