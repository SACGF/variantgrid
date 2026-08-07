""" MME Metrics API - https://github.com/ga4gh/mme-apis/blob/master/metrics-api.md

Counts cover `mme_eligible_classifications()`, so published metrics, inbound matching and
outbound submission all describe one dataset. Building them resolves every eligible
classification's genomic + ontology profile, so they are cached and refreshed nightly
rather than computed per request.
"""
from django.conf import settings
from django.core.cache import cache
from django.utils import timezone

from library.constants import DAY_SECS, HOUR_SECS
from mme.genes import canonical_gene_key
from mme.models import MMEInboundMatch, MMEInboundQuery
from mme.serializers.patient_profile import (
    classification_genomic_feature,
    classification_ontology_slots,
    mme_eligible_classifications,
)

MME_METRICS_CACHE_KEY = "mme_metrics"
# Longer than the nightly refresh, so a missed beat run serves stale metrics rather than
# recomputing inside a request.
MME_METRICS_CACHE_SECS = DAY_SECS + HOUR_SECS


def build_metrics() -> dict:
    """ Compute the metrics payload from scratch. Called by the nightly task. """
    num_cases = 0
    submitters: set = set()
    gene_total = 0
    unique_genes: set[str] = set()
    variant_total = 0
    unique_variants: set[tuple] = set()
    feature_total = 0
    unique_features: set[str] = set()
    feature_sets = 0
    cases_with_diagnosis = 0

    for cm in mme_eligible_classifications().select_related("classification__lab").iterator():
        classification = cm.classification
        num_cases += 1
        submitters.add(classification.lab_id)   # lab, the unit that opts in and is contacted

        for gf in (classification_genomic_feature(classification) or []):
            if gene := gf.get("gene"):
                gene_total += 1
                if gene_key := canonical_gene_key(gene):
                    unique_genes.add(gene_key)
            if variant := gf.get("variant"):
                variant_total += 1
                unique_variants.add((variant.get("assembly"), variant.get("referenceName"),
                                     variant.get("start"), variant.get("referenceBases"),
                                     variant.get("alternateBases")))

        features, disorders = classification_ontology_slots(classification)
        if features:
            feature_sets += 1
            feature_total += len(features)
            unique_features.update(f["id"] for f in features)
        if disorders:
            cases_with_diagnosis += 1

    return {
        "numberOfCases": num_cases,
        "numberOfSubmitters": len(submitters),
        "numberOfGenes": gene_total,
        "numberOfUniqueGenes": len(unique_genes),
        "numberOfVariants": variant_total,
        "numberOfUniqueVariants": len(unique_variants),
        "numberOfFeatures": feature_total,
        "numberOfUniqueFeatures": len(unique_features),
        "numberOfFeatureSets": feature_sets,
        "numberOfUniqueGenesMatched": _unique_genes_matched(),
        "numberOfCasesWithDiagnosis": cases_with_diagnosis,
        "numberOfRequestsReceived": MMEInboundQuery.objects.count(),
        "numberOfPotentialMatchesSent": MMEInboundMatch.objects.count(),
        "dateGenerated": timezone.now().isoformat(),
    }


def _unique_genes_matched() -> int:
    """ Distinct genes among our classifications we have actually returned to a peer. """
    unique = set()
    for classification in {m.classification for m in
                           MMEInboundMatch.objects.select_related("classification")}:
        for gf in (classification_genomic_feature(classification) or []):
            if (gene := gf.get("gene")) and (gene_key := canonical_gene_key(gene)):
                unique.add(gene_key)
    return len(unique)


def get_metrics(refresh: bool = False) -> dict:
    """ Cached metrics payload. `refresh=True` recomputes (nightly task). """
    if not settings.MME_ENABLED:
        # Eligibility is empty, so this is zeros over an empty queryset - cheap to build,
        # and skipping the cache means enabling MME doesn't serve stale zeros for a day.
        return build_metrics()
    if not refresh:
        if (cached := cache.get(MME_METRICS_CACHE_KEY)) is not None:
            return cached
    metrics = build_metrics()
    cache.set(MME_METRICS_CACHE_KEY, metrics, MME_METRICS_CACHE_SECS)
    return metrics
