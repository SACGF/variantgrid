""" Outbound: build a MatchMaker Exchange `patient` profile from a Classification.

This is *not* a DRF serializer (that framework shapes inbound HTTP); it is a plain
builder in the `ExportRow` spirit, producing a dict from a Classification. The MME
`patient` object has three information-bearing parts, populated from the classification:

    genomicFeatures  candidate variant + gene   <- classification's allele  (§6b)
    disorders        curated disease            <- condition terms (OMIM/Orphanet)  (§6a)
    features         observed phenotype (HPO)   <- condition HPO terms + linked patient (§6a)

A classification may populate any subset; MME requires at least one of features /
genomicFeatures.
"""
from django.conf import settings

from classification.enums.classification_enums import SpecialEKeys, ShareLevel
from classification.models.classification import ClassificationModification
from mme.contact import mme_contact_for_classification
from ontology.models import OntologyTerm, OntologyService, OntologySnake, OntologyTermRelation


def classification_ontology_slots(classification, include_patient_phenotype: bool = True) -> tuple[list[dict], list[dict]]:
    """ Route condition + linked-patient terms into MME (features_hpo, disorders).
        Each condition term is routed by its ontology service, not assumed to be disease.
        Derived terms carry a label suffix + `_derivedFrom` id for provenance.

        include_patient_phenotype gates the linked-patient HPO block (part 2). MME always
        includes it (the submission is curator-driven); Beacon passes False unless the
        linked Patient is itself shared publicly (§5.5 patient extra gate). """
    features: dict[str, dict] = {}   # keyed by HPO id to de-dup
    disorders: list[dict] = []
    disorder_ids: set[str] = set()

    def add_feature(term, observed="yes", derived_from=None):
        entry = {"id": term.pk, "label": term.name, "observed": observed}
        if derived_from:
            entry["label"] = f"{term.name} (via {derived_from.pk})"
            entry["_derivedFrom"] = derived_from.pk   # underscore = MME extension field
        features.setdefault(term.pk, entry)

    def add_disorder(mim_id: str):
        if mim_id not in disorder_ids:
            disorder_ids.add(mim_id)
            disorders.append({"id": mim_id})

    # 1. Condition under curation (curator-resolved; any ontology)
    resolved = classification.condition_resolution_obj          # cached_property
    for term in (resolved.terms if resolved else []):
        if term.ontology_service == OntologyService.HPO:
            add_feature(term)                                   # observed diagnosis term
        elif term.ontology_service == OntologyService.OMIM:
            add_disorder(f"MIM:{term.index}")
        elif term.ontology_service == OntologyService.MONDO:
            # EXACT alias to OMIM (lossless); MME has no MONDO slot
            if settings.MME_ONTOLOGY_SNAKE_EXACT and (omim := OntologyTermRelation.as_omim(term)):
                add_disorder(f"MIM:{omim.index}")
        # Orphanet maps straight through once that ontology service is supported

        # Approximate: expand a disease term to its typical HPO phenotypes (opt-in)
        if settings.MME_ONTOLOGY_PHENOTYPE_EXPANSION and term.ontology_service in (
                OntologyService.OMIM, OntologyService.MONDO):
            for snake in OntologySnake.snake_from(term, OntologyService.HPO):
                add_feature(snake.leaf_term, observed="no", derived_from=term)

    # 2. Linked patient phenotype (auto-matched; HPO only, genuinely observed)
    sample = classification.sample if include_patient_phenotype else None  # nullable FK
    patient = getattr(sample, "patient", None) if sample else None
    ptp = getattr(patient, "patient_text_phenotype", None) if patient else None
    if ptp and ptp.phenotype_description:
        for term in OntologyTerm.objects.filter(
                pk__in=ptp.phenotype_description.get_ontology_term_ids(),
                ontology_service=OntologyService.HPO):
            add_feature(term)

    return list(features.values()), disorders


def classification_genomic_feature(classification) -> list[dict] | None:
    """ Build the MME genomicFeatures entry from a classification's allele.
        MME variant coords are 0-based; ours are 1-based -> subtract 1.
        Returns None when there is neither a resolvable gene symbol nor variant. """
    gene_symbol = classification.get(SpecialEKeys.GENE_SYMBOL)

    variant = None
    try:
        genome_build = classification.get_genome_build()
        variant = classification.get_variant_for_build(genome_build)
    except ValueError:
        genome_build = None

    if not gene_symbol and variant is None:
        return None

    feature: dict = {}
    if gene_symbol:
        feature["gene"] = {"id": str(gene_symbol)}
    if variant is not None:
        vc = variant.coordinate                                # VariantCoordinate (1-based)
        if vc.ref and vc.alt and not vc.is_symbolic:
            feature["variant"] = {
                "assembly": genome_build.name,                 # "GRCh37" / "GRCh38"
                "referenceName": vc.chrom.replace("chr", ""),
                "start": vc.position - 1,                       # 1-based -> 0-based
                "referenceBases": vc.ref,
                "alternateBases": vc.alt,
            }
    if not feature:
        return None
    return [feature]


def mme_eligible_classifications():
    """ Latest PUBLIC ('3rd Party Databases') published, non-withdrawn modifications
        whose owning lab has opted into MME. Eligibility is a three-layer AND: node
        (settings.MME_ENABLED) x lab (Lab.mme_enabled) x record (share_level=PUBLIC,
        not withdrawn). share_level=PUBLIC is the exact consent signal for an external
        DB like MME; ALL_USERS is internal-only and must be excluded. """
    return (ClassificationModification.objects
            .filter(is_last_published=True,
                    share_level=ShareLevel.PUBLIC.value,
                    classification__withdrawn=False,
                    classification__lab__mme_enabled=True)
            .select_related("classification"))


def build_patient_profile(submission) -> dict:
    """ Assemble the MME `patient` object from the submission's classification.
        The request-level `disclaimer`/`terms` are siblings of `patient` in the
        envelope, so they are attached by the client (§7), not here. """
    classification = submission.classification
    contact = mme_contact_for_classification(classification)
    if not (contact.get("name") and contact.get("href")):
        lab = getattr(classification, "lab", None)
        lab_name = getattr(lab, "name", None) or getattr(classification, "lab_id", None)
        raise ValueError(
            f"Cannot submit to MME: no resolvable contact for lab '{lab_name}' "
            f"and settings.MME_CONTACT (name + href) is not configured")

    genomic_features = classification_genomic_feature(classification)      # §6b
    features, disorders = classification_ontology_slots(classification)    # §6a

    profile = {
        "id": submission.external_patient_id,   # opaque, stable, non-PII
        "contact": contact,
        "species": "NCBITaxon:9606",
    }
    if genomic_features:
        profile["genomicFeatures"] = genomic_features
    if disorders:
        profile["disorders"] = disorders
    if features:
        profile["features"] = features
    if not settings.MME_ENABLED_PRODUCTION_SUBMIT:   # keep test-mode until certified
        profile["test"] = True

    # MME requires at least one of features / genomicFeatures
    if "features" not in profile and "genomicFeatures" not in profile:
        raise ValueError("Cannot submit to MME: no HPO features or genomic features")
    return profile
