""" Static framework metadata for the Beacon v2 framework endpoints, built from
    settings.BEACON_CONFIG (§6). These are mostly-static identity/config documents. """
from django.conf import settings

from beacon.response import GRANULARITIES

# Canonical Beacon v2 schema URIs. The configuration/map `response` objects require a
# `$schema` property (beaconMapSchema / beaconConfigurationSchema), used by clients for
# schema discovery; omitting it fails spec validation (EGA beacon-verifier).
_SCHEMA_BASE = "https://raw.githubusercontent.com/ga4gh-beacon/beacon-v2/main/framework/json/configuration"
BEACON_MAP_SCHEMA = f"{_SCHEMA_BASE}/beaconMapSchema.json"
BEACON_CONFIGURATION_SCHEMA = f"{_SCHEMA_BASE}/beaconConfigurationSchema.json"

# The two datasets we serve (§5.5), advertised in service-info / entry types.
OBSERVATIONS_DATASET_ID = "variantgrid_observations"
CLASSIFICATIONS_DATASET_ID = "variantgrid_classifications"

DATASETS = [
    {
        "id": OBSERVATIONS_DATASET_ID,
        "name": "VariantGrid sample observations",
        "description": "Presence and zygosity counts of the variant across sample genotypes, "
                       "scoped to the datasets the requester may read.",
    },
    {
        "id": CLASSIFICATIONS_DATASET_ID,
        "name": "VariantGrid public classifications",
        "description": "Variants carrying a published, publicly-shared ACMG classification "
                       "(the same consented set shared to ClinVar / MatchMaker Exchange).",
    },
]


def _organization() -> dict:
    org = settings.BEACON_CONFIG.get("organization", {})
    return {
        "id": org.get("id", ""),
        "name": org.get("name", ""),
        "welcomeUrl": org.get("welcome_url", ""),
        "contactUrl": org.get("contact_url", ""),
    }


def beacon_info() -> dict:
    """ GET / and /info - Beacon identity/metadata. """
    config = settings.BEACON_CONFIG
    return {
        "id": config["beacon_id"],
        "name": config["name"],
        "apiVersion": config["api_version"],
        "environment": config.get("environment", "prod"),
        "organization": _organization(),
        "datasets": DATASETS,
    }


def service_info() -> dict:
    """ GET /service-info - GA4GH service-info profile. """
    config = settings.BEACON_CONFIG
    org = _organization()
    return {
        "id": config["beacon_id"],
        "name": config["name"],
        "type": {
            "group": "org.ga4gh",
            "artifact": "beacon",
            "version": config["api_version"],
        },
        "description": "VariantGrid GA4GH Beacon v2 genomic data-sharing endpoint.",
        "organization": {"name": org["name"], "url": org["welcomeUrl"]},
        "contactUrl": org["contactUrl"],
        "version": config["api_version"],
        "environment": config.get("environment", "prod"),
    }


def configuration() -> dict:
    """ GET /configuration - Beacon configuration object. """
    config = settings.BEACON_CONFIG
    return {
        "$schema": BEACON_CONFIGURATION_SCHEMA,
        "maturityAttributes": {"productionStatus": "DEV"},
        "securityAttributes": {
            "defaultGranularity": config.get("default_granularity", "boolean"),
            "securityLevels": ["PUBLIC"],
        },
        "entryTypes": _entry_types(),
    }


def _entry_types() -> dict:
    # Phase 1: only the genomicVariant entry type is served.
    return {
        "genomicVariant": {
            "id": "genomicVariant",
            "name": "Genomic Variant",
            "ontologyTermForThisType": {"id": "SO:0001059", "label": "sequence_alteration"},
            "partOfSpecification": "Beacon v2.0.0",
            "defaultSchema": {
                "id": "ga4gh-beacon-variant-v2.0.0",
                "name": "Default schema for a genomic variant",
            },
        }
    }


def entry_types() -> dict:
    """ GET /entry_types - supported entry types. """
    return {"entryTypes": _entry_types()}


def filtering_terms() -> dict:
    """ GET /filtering_terms - supported ontology filters (phase 1: minimal / none). """
    return {"filteringTerms": [], "resources": []}


def endpoint_map(g_variants_url: str) -> dict:
    """ GET /map - endpoint map. rootUrl/singleEntryUrl must be absolute URIs (spec
        `format: uri`); relative paths crash strict clients (e.g. the EGA beacon-verifier
        whose URL parser rejects a base-less relative URL). """
    return {
        "$schema": BEACON_MAP_SCHEMA,
        "endpointSets": {
            "genomicVariant": {
                "entryType": "genomicVariant",
                "rootUrl": g_variants_url,
                "singleEntryUrl": f"{g_variants_url}/{{id}}",
            }
        }
    }


def supported_granularities() -> tuple:
    return GRANULARITIES
