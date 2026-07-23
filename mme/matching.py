""" Inbound similarity scoring for MatchMaker Exchange /match.

v1 is deliberately simple: overlap of candidate genes, disorders and HPO features
between the querying patient and each of our ShareLevel.PUBLIC classifications — the
same consented set we are willing to submit outward, so inbound and outbound expose a
consistent dataset. Richer OntologySnake semantic similarity can be layered on later.
"""
import logging

from django.conf import settings

from mme.contact import mme_contact_for_classification
from mme.serializers.patient_profile import (
    mme_eligible_classifications,
    classification_ontology_slots,
    classification_genomic_feature,
)

# v1 iterates our PUBLIC classifications in Python; cap the scan so an inbound query
# can't run unbounded. Truncation is logged (never silently dropped).
MAX_CLASSIFICATIONS_SCANNED = 2000
# Only return candidates with at least this combined similarity.
MIN_SCORE = 0.01
MAX_RESULTS = 20

# Relative weight of each evidence type in the combined score.
_GENE_WEIGHT = 0.5
_DISORDER_WEIGHT = 0.25
_FEATURE_WEIGHT = 0.25


def _genes_from_profile(genomic_features) -> set[str]:
    genes = set()
    for gf in (genomic_features or []):
        gene = (gf.get("gene") or {}).get("id")
        if gene:
            genes.add(gene.upper())
    return genes


def _jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def _extract_query(patient: dict) -> tuple[set[str], set[str], set[str]]:
    genes = _genes_from_profile(patient.get("genomicFeatures"))
    features = {f.get("id") for f in patient.get("features", []) if f.get("id")}
    disorders = {d.get("id") for d in patient.get("disorders", []) if d.get("id")}
    return genes, features, disorders


def _our_patient_object(classification, genomic_features, features, disorders) -> dict:
    """ Non-PII patient object describing one of our classifications, in MME shape.
        Contact follows the classification's lab; if that can't be resolved we fall
        back to the server contact as-is rather than drop the match. """
    patient = {
        "id": f"vg:{classification.pk}",   # opaque, stable, non-PII
        "contact": mme_contact_for_classification(classification) or (settings.MME_CONTACT or {}),
        "species": "NCBITaxon:9606",
    }
    if genomic_features:
        patient["genomicFeatures"] = genomic_features
    if disorders:
        patient["disorders"] = disorders
    if features:
        patient["features"] = features
    return patient


def find_matches(patient: dict) -> list[dict]:
    """ Score the inbound `patient` against our PUBLIC classifications.
        Returns MME result objects: {"score": {"patient": x}, "patient": {...}}. """
    query_genes, query_features, query_disorders = _extract_query(patient)

    qs = mme_eligible_classifications()
    scored = []
    scanned = 0
    for cm in qs.iterator():
        if scanned >= MAX_CLASSIFICATIONS_SCANNED:
            logging.warning("MME find_matches: scan capped at %d classifications; "
                            "some candidates were not evaluated", MAX_CLASSIFICATIONS_SCANNED)
            break
        scanned += 1
        classification = cm.classification

        genomic_features = classification_genomic_feature(classification)
        features, disorders = classification_ontology_slots(classification)

        our_genes = _genes_from_profile(genomic_features)
        our_features = {f["id"] for f in features}
        our_disorders = {d["id"] for d in disorders}

        score = (_GENE_WEIGHT * _jaccard(query_genes, our_genes)
                 + _DISORDER_WEIGHT * _jaccard(query_disorders, our_disorders)
                 + _FEATURE_WEIGHT * _jaccard(query_features, our_features))
        if score >= MIN_SCORE:
            scored.append((score, _our_patient_object(
                classification, genomic_features, features, disorders)))

    scored.sort(key=lambda t: t[0], reverse=True)
    return [{"score": {"patient": round(score, 4)}, "patient": obj}
            for score, obj in scored[:MAX_RESULTS]]
