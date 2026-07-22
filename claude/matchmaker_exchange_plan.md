# MatchMaker Exchange (MME) Integration Plan

Implementation plan for GitHub issue **#1662** — integrate MatchMaker Exchange for
rare-disease patient matching. Split from #40.

## 1. Background

MatchMaker Exchange is a **patient-centric, phenotype-aware** federation protocol
(distinct from Beacon, which is variant-centric). You POST a *patient profile*
(candidate gene(s)/variant(s) + HPO terms + a contact) to a remote node; each node
runs its own similarity matching and returns ranked candidate patients plus a contact
for follow-up.

Active nodes: DECIPHER, GeneMatcher, seqr, PhenomeCentral, RD-Connect, MyGene2,
PatientMatcher, IRUD.

- Spec: [`ga4gh/mme-apis`](https://github.com/ga4gh/mme-apis) (v1.1, Nov 2017). Live
  endpoint list on the repo wiki.
- **Connection is bilateral, not centrally gated.** Per the API repo: *"to query most
  MME services, you need to request an authentication token from that service"* — nodes
  exchange tokens directly. Becoming a *recognised, discoverable MME Service* (and
  serving inbound reciprocal matches) is the layer with central oversight — the Steering
  Committee's Service Requirements + agreements (§10).

### Protocol summary

`POST <node>/match`

Headers:
- `Content-Type: application/vnd.ga4gh.matchmaker.v1.1+json`
- `Accept: application/vnd.ga4gh.matchmaker.v1.1+json`
- `X-Auth-Token: <per-node token>`

Request body:
```json
{
  "patient": {
    "id": "unique-within-our-node",
    "label": "optional, no PII",
    "contact": { "name": "Lab X", "href": "mailto:curator@lab.org", "institution": "..." },
    "species": "NCBITaxon:9606",
    "sex": "FEMALE|MALE|OTHER",
    "disorders": [
      { "id": "MIM:143100" }
    ],
    "features": [
      { "id": "HP:0001250", "label": "Seizure", "observed": "yes" }
    ],
    "genomicFeatures": [
      {
        "gene": { "id": "BRCA1" },
        "variant": {
          "assembly": "GRCh38",
          "referenceName": "17",
          "start": 43044294,
          "referenceBases": "A",
          "alternateBases": "G"
        },
        "zygosity": 1,
        "type": { "id": "SO:0001583", "label": "missense_variant" }
      }
    ]
  }
}
```

Rules:
- At least one of `features` (HPO) or `genomicFeatures` is mandatory (both preferred).
- `variant.start` / `end` are **0-based** (BED-style). Our internal coordinates are
  1-based — convert on the way out.
- Test queries carry `"test": true` in the patient object.

Response body:
```json
{
  "results": [
    { "score": { "patient": 0.83 }, "patient": { /* same shape as request */ } }
  ]
}
```

## 2. Design principles for this build

- **Standalone `mme/` Django app.** Do **not** build on `sync/`'s `SyncRunner`
  framework — that is Shariant-classification-specific and not worth generalising.
  Do **not** depend on any `sync/alissa/*` code — Alissa is abandoned and will be
  removed. We *copy the good patterns* (config-in-settings, `requests`-based client,
  admin notification on failure) into self-contained MME code.
- **Reuse read-only helpers that are already general:** ontology term access,
  HGVS/ClinGen variant formatting, `Patient`/phenotype models, DRF `APIView` for the
  inbound endpoint.
- **Consent = `share_level = ShareLevel.PUBLIC`.** MME shares variant + condition data
  outward, so the submission unit is a **publicly-shared Classification**. `PUBLIC` is
  labelled *"3rd Party Databases"* and is already the bar for ClinVar submission — a
  classification set to PUBLIC is, by the lab's own decision, cleared to go to external
  databases, so no additional patient-level consent flag is needed.
  (`Patient.research_consent` is for *internal* database mining and is **not** the MME
  gate.)
- **Explicit and manual by default.** First release is a manual "submit to MME" action
  on a PUBLIC classification, not an automatic background sweep. Auto-submit can come
  later behind a setting.

## 3. New app: `mme/`

```
mme/
  __init__.py
  apps.py
  models.py
  admin.py
  client.py                # outbound HTTP client (copies ServerAuth pattern, no sync dep)
  serializers/
    __init__.py
    patient_profile.py     # Patient -> MME patient JSON (features + genomicFeatures)
    match_request.py       # DRF serializers for inbound /match
  matching.py              # inbound similarity scoring
  views_rest.py            # inbound POST /match  (DRF APIView)
  views.py                 # UI: submit patient, view match results
  urls.py
  migrations/
  templates/mme/
  tests/
    test_patient_profile.py
    test_client.py
    test_inbound_match.py
```

Register in `settings.INSTALLED_APPS` (via the appropriate
`variantgrid/settings/components/default_settings.py` list), and add
`mme.urls` to the project URL conf.

## 4. Configuration & credentials

MME uses a **static per-node auth token** (not OAuth), so this is simpler than
`library/oauth.py:ServerAuth`. Only the **auth tokens** are secret and go through
`get_secret()` (see `default_settings.py:215` — `CLINGEN_ALLELE_REGISTRY_*`); everything
else (enabled flag, our public contact details, node URLs, ontology toggles, disclaimer
text) is a plain Python setting.

`variantgrid/settings/components/default_settings.py`:
```python
# --- MatchMaker Exchange -------------------------------------------------
MME_ENABLED = False

# Our own contact details, sent as the `contact` block on every outbound patient.
# Public info, not secret.
MME_CONTACT = {
    "name": "SA Pathology Genetics",
    "href": "mailto:mme@sapath.org",
    "institution": "SA Pathology",
}

# Remote nodes we can submit to: public base_url + api_version as plain config; the
# token THEY issued to us is the only secret, pulled via get_secret().
MME_NODES = {
    "decipher": {"base_url": "https://...", "api_version": "1.1",
                 "token": get_secret("MME.decipher_token")},
    "genematcher": {"base_url": "https://...", "api_version": "1.1",
                    "token": get_secret("MME.genematcher_token")},
}

# Token that remote nodes must present to call OUR inbound /match endpoint. Secret.
MME_INBOUND_TOKEN = get_secret("MME.inbound_token")

# Follow EXACT OntologyTermRelations to alias each condition term into the ontology
# form MME expects (e.g. a MONDO diagnosis -> its EXACT-equivalent OMIM for `disorders`).
# These are lossless, same-concept aliases. On by default; labs that curate in MONDO
# need this to produce any `disorders`. Off = send only terms already in MME's ontology.
MME_ONTOLOGY_SNAKE_EXACT = True

# Expand a curated disease term to its *typical* HPO phenotypes (disease -> phenotype is
# NOT an equivalence, so this is approximate and not "observed" in the patient). Off by
# default; turn on to widen the phenotype match surface, accepting looser matches.
MME_ONTOLOGY_PHENOTYPE_EXPANSION = False

# Standing disclaimer sent as the request-level `disclaimer`, so receiving curators
# know some phenotype terms may be ontology-derived and matches may be approximate.
MME_DISCLAIMER = (
    "Some phenotype (HPO) terms may be derived from the curated diagnosis via ontology "
    "cross-references and were not necessarily directly observed; terms carry a "
    "`_derivedFrom` field where so derived. Please contact us to confirm before acting "
    "on a match."
)
```

Deployments override the non-secret values in their per-host env settings file; the
`config/settings_config.json` template gains only the secret **token** keys
(`MME.*_token`, `MME.inbound_token`).

## 5. Models (`mme/models.py`)

We track what we submitted and what came back, for audit and follow-up. `Patient`
already extends `ExternallyManagedModel` / `GuardianPermissionsMixin`
(`patients/models.py:122`); MME records hang off it.

```python
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.utils import timezone

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Patient


class MMESubmissionStatus(models.TextChoices):
    DRAFT = "D", "Draft"
    SUBMITTED = "S", "Submitted"
    ERROR = "E", "Error"


class MMESubmission(models.Model):
    """ A single outbound submission of one Classification to one remote MME node.
        The MME 'patient' profile is synthesised from the classification: candidate
        variant from its allele, disorders from its condition-under-curation, and HPO
        features from a linked patient (if any). The unit is the classification, not
        the Patient, because classifications are the consented (share_level=PUBLIC),
        structured unit and are frequently not linked to a Patient at all. """
    classification = models.ForeignKey("classification.Classification", on_delete=models.CASCADE)
    node_id = models.CharField(max_length=64)          # key into settings.MME_NODES
    # Opaque, stable, non-PII id we send as patient.id (e.g. the classification's
    # lab record id / UUID).
    external_patient_id = models.CharField(max_length=255)
    status = models.CharField(max_length=1, choices=MMESubmissionStatus.choices,
                              default=MMESubmissionStatus.DRAFT)
    request_json = models.JSONField(null=True)          # exact profile we POSTed
    response_json = models.JSONField(null=True)         # raw /match response
    error = models.TextField(null=True, blank=True)
    created = models.DateTimeField(default=timezone.now)
    submitted = models.DateTimeField(null=True, blank=True)
    submitted_by = models.ForeignKey(User, null=True, on_delete=models.SET_NULL)

    class Meta:
        unique_together = ("classification", "node_id")


class MMEMatchResult(models.Model):
    """ One candidate patient returned in a /match response (inbound or outbound). """
    submission = models.ForeignKey(MMESubmission, null=True, on_delete=models.CASCADE)
    score = models.FloatField()
    matched_patient_id = models.CharField(max_length=255)   # remote node's patient id
    contact_name = models.TextField(null=True)
    contact_href = models.TextField(null=True)
    patient_json = models.JSONField()                       # full returned patient obj
    created = models.DateTimeField(default=timezone.now)
```

If we serve inbound matches, add an `MMEInboundQuery` audit row (who queried, when,
what profile, how many of our patients we returned) — same shape, keeps a compliance
trail.

## 6. Outbound: classification → MME profile serializer

This is the core mapping. It is *not* a DRF serializer (that framework is for shaping
inbound HTTP); it is a plain builder in the `ExportRow` spirit, producing a dict from a
**classification**. The MME `patient` object has three information-bearing parts:

| MME field | Meaning | Ontology | VariantGrid source |
| --- | --- | --- | --- |
| `genomicFeatures` | candidate variant + gene | — | classification's allele (§6b) |
| `disorders` | curated disease | OMIM / Orphanet | condition terms of disease type (§6a) |
| `features` | observed phenotype | HPO | condition HPO terms **+** linked patient's phenotype (§6a) |

A classification may populate any subset. The condition under curation is the primary
ontology source and is often present **even without a linked sample/patient**.

### 6a. Ontology terms → `features` (HPO) + `disorders` (disease)

**The condition under curation is not disease-only.** A classification's condition
accepts *any* ontology — structured or free text, resolved and linked to real terms — so
`ConditionResolved.terms` can hold HPO phenotype terms, MONDO/OMIM/Orphanet disease
terms, or a mix. The mapping therefore **routes each term by its ontology service** into
the correct MME slot, rather than assuming disease.

Condition terms live in `Classification.condition_resolution`
(`classification/models/classification.py:514`, a `ConditionResolvedDict` JSON blob),
surfaced as a `ConditionResolved` object via `condition_resolution_obj`
(`classification.py:637`); `ConditionResolved.terms` is a `List[OntologyTerm]`
(`classification.py:235`).

**Routing** each condition term by ontology service, using `OntologyTermRelation` /
`OntologySnake` (`ontology/models/models_ontology.py`):
- **HPO term** → MME `features` (`HP:#######`, direct).
- **OMIM term** → MME `disorders` (`MIM:######`, direct).
- **MONDO term** → MME `disorders` via the EXACT alias crosswalk `OntologySnake.as_omim()`
  (`models_ontology.py:697`, which walks an `EXACT` `OntologyTermRelation`). A MONDO
  disease and its EXACT-linked OMIM are the *same concept*, so this is lossless. Gated
  on **`MME_ONTOLOGY_SNAKE_EXACT`** (default on) — with it off, a MONDO-only condition
  yields no `disorders` (MME has no MONDO slot).

Two clearly-separated ontology behaviours, each with its own setting:
- **EXACT aliasing** (`MME_ONTOLOGY_SNAKE_EXACT`, default **on**): lossless same-concept
  crosswalks so curated terms reach the ontology MME expects (above).
- **Phenotype expansion** (`MME_ONTOLOGY_PHENOTYPE_EXPANSION`, default **off**):
  approximate. A disease term does **not** map to HPO by equivalence — disease and
  phenotype are different concept types — so to add phenotype `features` from a diagnosis
  we expand via `OntologySnake.snake_from(term, OntologyService.HPO)`
  (`models_ontology.py:1100`), taking the leaf HPO terms. These are the disease's
  *typical* phenotypes, not observed in this patient — hence off by default.

**Provenance on derived terms.** MME has no per-feature comment field, and `observed`
accepts only `"yes"`/`"no"` — no "inferred". So an expanded HPO term is honestly *not
observed*. We annotate every derived term so the receiving curator/algorithm can weight
or ignore it:
- per-feature `label` suffix, e.g. `"Seizure (via MONDO:0005027)"`;
- an underscore-extension `"_derivedFrom": "MONDO:0005027"` (spec-sanctioned custom
  field for tracking/reference);
- and a request-level `disclaimer` (`settings.MME_DISCLAIMER`, added in §6d).

Patient phenotype supplements the HPO `features` directly (these *were* observed):
`Classification.sample` (`classification.py:469`, nullable) → `Sample.patient` → the
NLP-matched free-text phenotype pipeline (`annotation/models/models_phenotype_match.py`);
`PhenotypeDescription.get_ontology_term_ids()` (`models_phenotype_match.py:59`) returns
matched term ids (filter to HPO).

```python
from django.conf import settings

from ontology.models import OntologyTerm, OntologyService, OntologySnake


def classification_ontology_slots(classification) -> tuple[list[dict], list[dict]]:
    """ Route condition + patient terms into MME (features_hpo, disorders).
        Derived terms carry a label suffix + `_derivedFrom` for provenance. """
    features: dict[str, dict] = {}   # keyed by HPO id to de-dup
    disorders: list[dict] = []

    def add_feature(term, observed="yes", derived_from=None):
        entry = {"id": term.pk, "label": term.name, "observed": observed}
        if derived_from:
            entry["label"] = f"{term.name} (via {derived_from.pk})"
            entry["_derivedFrom"] = derived_from.pk   # underscore = MME extension field
        features.setdefault(term.pk, entry)

    # 1. Condition under curation (curator-resolved; any ontology)
    resolved = classification.condition_resolution_obj          # classification.py:637
    for term in (resolved.terms if resolved else []):
        if term.ontology_service == OntologyService.HPO:
            add_feature(term)                                   # observed diagnosis term
        elif term.ontology_service == OntologyService.OMIM:
            disorders.append({"id": f"MIM:{term.index}"})
        elif term.ontology_service == OntologyService.MONDO:
            # EXACT alias to OMIM (lossless); MME has no MONDO slot
            if settings.MME_ONTOLOGY_SNAKE_EXACT and (omim := OntologySnake.as_omim(term)):
                disorders.append({"id": f"MIM:{omim.index}"})
        # Orphanet maps straight through once that ontology service is supported

        # Approximate: expand a disease term to its typical HPO phenotypes (opt-in)
        if settings.MME_ONTOLOGY_PHENOTYPE_EXPANSION and term.ontology_service in (
                OntologyService.OMIM, OntologyService.MONDO):
            for snake in OntologySnake.snake_from(term, OntologyService.HPO):
                add_feature(snake.leaf_term, derived_from=term)

    # 2. Linked patient phenotype (auto-matched; HPO only, genuinely observed)
    sample = classification.sample                              # classification.py:469
    patient = getattr(sample, "patient", None) if sample else None
    ptp = getattr(patient, "patient_text_phenotype", None) if patient else None
    if ptp and ptp.phenotype_description:
        for term in OntologyTerm.objects.filter(
                pk__in=ptp.phenotype_description.get_ontology_term_ids(),
                ontology_service=OntologyService.HPO):
            add_feature(term)

    return list(features.values()), disorders
```

> **Design note — term confirmation.** Curator-resolved condition terms and observed
> patient HPO are high confidence; `snake_from` phenotype expansions are approximate and
> carry `_derivedFrom`. The submit UI (§9) shows the assembled `features` + `disorders`,
> flags each by provenance, and lets the curator prune before POSTing. Persist the
> confirmed set on `MMESubmission.request_json`.

### 6b. Genomic features (candidate variant from the classification's allele)

**Are candidate variants required?** No. Per the MME v1.1 spec a profile needs **at
least one** of `features` (HPO) or `genomicFeatures`; both are preferred. But the
candidate variant is the strongest signal — several major nodes (GeneMatcher in
particular) match primarily on **candidate gene** — and since the submission unit *is* a
classification, its variant is always present.

The classification's allele → gene + coordinates, reusing existing HGVS/coordinate/
ClinGen machinery (`snpdb/models/models_variant.py`, `genes/hgvs/`,
`snpdb/models/models_clingen_allele.py`). The ClinVar export formatter
(`classification/views/exports/classification_export_formatter_clinvar.py`) is the
closest existing "curated variant → external genomics API" template.

```python
def classification_genomic_feature(classification) -> dict | None:
    """ Build the MME genomicFeatures entry from a classification's allele.
        MME variant coords are 0-based; ours are 1-based -> subtract 1. """
    genome_build = classification.get_genome_build()
    variant = classification.variant                       # classification.py
    if variant is None:
        return None
    vc = variant.coordinate                                # VariantCoordinate (1-based)
    gene_symbol = classification.get_gene_symbol()         # e.g. via c_hgvs

    feature = {"gene": {"id": str(gene_symbol)}}
    if vc.ref and vc.alt and not _is_symbolic(vc.alt):
        feature["variant"] = {
            "assembly": genome_build.name,                 # "GRCh37" / "GRCh38"
            "referenceName": vc.chrom.replace("chr", ""),
            "start": vc.position - 1,                       # 1-based -> 0-based
            "referenceBases": vc.ref,
            "alternateBases": vc.alt,
        }
    return [feature]
```

### 6c. Eligibility gate — only `ShareLevel.PUBLIC` classifications may be sent

The submission unit is a classification, and exposing it to an external federation is
governed **solely** by its share level — the same bar as ClinVar submission. No
patient-level consent flag is involved (`Patient.research_consent` is for internal DB
mining and does **not** gate MME).

`ShareLevel.PUBLIC` is defined with the label **"3rd Party Databases"**
(`classification/enums/classification_enums.py:344`, index 4 — the highest level). MME
*is* a third-party database, so `PUBLIC` is the semantically exact — and only
appropriate — gate:

- **Eligible ⇔ `share_level = ShareLevel.PUBLIC`.** A classification set to PUBLIC is,
  by the lab's own decision, already cleared to flow to external databases.
- **Not** `ShareLevel.ALL_USERS` ("App Users" = logged-in VariantGrid users on this
  instance; still internal). Note `shared_only=True` on the existing modification query
  filters `DISCORDANT_LEVEL_KEYS = {ALL_USERS, PUBLIC}`
  (`classification/models/classification.py:2419`, `classification_enums.py:347`) —
  that is **too broad** for MME. Filter to `PUBLIC` explicitly.
- Only **published** modifications are visible outside the owning lab (share_level lives
  on `ClassificationModification`), so query the latest *published* modification.

```python
from classification.enums.classification_enums import ShareLevel
from classification.models.classification import ClassificationModification


def mme_eligible_classifications():
    """ Latest PUBLIC ('3rd Party Databases') published, non-withdrawn modifications.
        share_level=PUBLIC is the exact consent signal for an external DB like MME;
        ALL_USERS is internal-only and must be excluded. """
    return (ClassificationModification.objects
            .filter(is_last_published=True,
                    share_level=ShareLevel.PUBLIC.key,
                    classification__withdrawn=False)
            .select_related("classification"))
```

### 6d. Assemble the profile

```python
from django.conf import settings


def build_patient_profile(submission) -> dict:
    classification = submission.classification
    genomic_features = classification_genomic_feature(classification)  # §6b
    features, disorders = classification_ontology_slots(classification)  # §6a

    profile = {
        "id": submission.external_patient_id,   # opaque, stable, non-PII
        "contact": settings.MME_CONTACT,
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
```

The request-level `disclaimer` / `terms` are **siblings of `patient`** in the envelope,
not part of the patient object, so they are attached where the request body is assembled
(the client, §7) rather than here.

## 7. Outbound client (`mme/client.py`)

Self-contained `requests` client — copies the *shape* of `library/oauth.py:ServerAuth`
(a small object wrapping `requests.get/post` + auth header) without importing it or any
`sync/` code. MME auth is a single static header, so no OAuth machinery.

```python
import requests
from django.conf import settings
from django.utils import timezone

from library.constants import MINUTE_SECS
from library.log_utils import AdminNotificationBuilder
from mme.models import MMESubmissionStatus, MMEMatchResult
from mme.serializers.patient_profile import build_patient_profile


class MMEClient:
    """ Minimal MME node client. One instance per remote node. """

    def __init__(self, node_id: str):
        node = settings.MME_NODES[node_id]
        self.node_id = node_id
        self.base_url = node["base_url"].rstrip("/")
        self.token = node["token"]
        self.api_version = node.get("api_version", "1.1")

    @property
    def _headers(self) -> dict:
        content_type = f"application/vnd.ga4gh.matchmaker.v{self.api_version}+json"
        return {
            "Content-Type": content_type,
            "Accept": content_type,
            "X-Auth-Token": self.token,
        }

    def match(self, patient_profile: dict, timeout: int = MINUTE_SECS) -> dict:
        # disclaimer/terms are message-level (siblings of `patient`), telling the
        # receiving curator some phenotype terms may be ontology-derived (§6a).
        body = {"patient": patient_profile}
        if settings.MME_DISCLAIMER:
            body["disclaimer"] = settings.MME_DISCLAIMER
        response = requests.post(
            url=f"{self.base_url}/match",
            json=body,
            headers=self._headers,
            timeout=timeout,
        )
        response.raise_for_status()
        return response.json()


def submit(submission) -> None:
    """ Build profile, POST, persist results. Notifies admins on failure
        (mirrors the AlissaUploadSyncer failure-reporting pattern, no sync dep). """
    profile = build_patient_profile(submission)
    submission.request_json = profile
    try:
        data = MMEClient(submission.node_id).match(profile)
    except Exception as e:
        submission.status = MMESubmissionStatus.ERROR
        submission.error = str(e)
        submission.save()
        nb = AdminNotificationBuilder("MME submission failed")
        nb.add_markdown(f"Classification {submission.classification_id} -> node "
                        f"`{submission.node_id}`: {e}")
        nb.send()
        raise

    submission.response_json = data
    submission.status = MMESubmissionStatus.SUBMITTED
    submission.submitted = timezone.now()
    submission.save()

    for result in data.get("results", []):
        contact = result.get("patient", {}).get("contact", {})
        MMEMatchResult.objects.create(
            submission=submission,
            score=result.get("score", {}).get("patient", 0.0),
            matched_patient_id=result.get("patient", {}).get("id", ""),
            contact_name=contact.get("name"),
            contact_href=contact.get("href"),
            patient_json=result.get("patient", {}),
        )
```

Run submissions on a worker queue (`@app.task(queue='web_workers')`) so the UI POST
returns immediately; the client call is network-bound.

## 8. Inbound: serving `/match` (reciprocity requirement)

MME membership requires we also answer queries. Build with DRF `APIView` — same pattern
as `ontology/views_rest.py`.

`mme/views_rest.py`:
```python
from django.conf import settings
from rest_framework.exceptions import AuthenticationFailed, ParseError
from rest_framework.response import Response
from rest_framework.views import APIView

from mme.matching import find_matches
from mme.models import MMEInboundQuery


class MMEMatchView(APIView):
    """ Inbound MME /match endpoint. Authenticated by the shared X-Auth-Token
        we issue to peer nodes, NOT by a VariantGrid user session. """
    authentication_classes = []          # token-header auth, not session/DRF user
    permission_classes = []

    def post(self, request, *args, **kwargs):
        token = request.headers.get("X-Auth-Token")
        if not token or token != settings.MME_INBOUND_TOKEN:
            raise AuthenticationFailed("Invalid MME auth token")

        patient = (request.data or {}).get("patient")
        if not patient or not (patient.get("features") or patient.get("genomicFeatures")):
            raise ParseError("patient with features or genomicFeatures required")

        results = find_matches(patient)          # §8b similarity scoring
        MMEInboundQuery.objects.create(
            request_json=request.data,
            num_results=len(results),
        )
        return Response({"results": results})
```

> **Auth carve-out note.** This endpoint is authenticated by the MME shared token, not
> a logged-in VariantGrid user, so it must be exempted from
> `GlobalLoginRequiredMiddleware` (add its path to the public-URL allowlist) — this is
> intentional, like other machine endpoints. The static-token check above is the
> access control. Rate-limit and log every call (`MMEInboundQuery`).

`mme/matching.py` — v1 scoring can be deliberately simple (`OntologySnake` semantic
similarity over HPO terms + shared candidate genes/disorders), matching against our
`ShareLevel.PUBLIC` classifications only (§6c) — the same set we are willing to submit
outward, so inbound and outbound expose a consistent, consented dataset.
`OntologySnake` (`ontology/models/models_ontology.py:951`) already does gene↔term
traversal for phenotype-driven scoring.

## 9. UI (`mme/views.py` + templates)

- **Classification page:** an "MatchMaker Exchange" panel, shown only when the
  classification is `ShareLevel.PUBLIC` (§6c). It displays the assembled profile —
  candidate gene/variant, curated disorders (from condition under curation), and HPO
  features (from a linked patient, if any) — and a **Submit to MME** button per
  configured node.
- **Submit flow:** GET builds a draft `MMESubmission`, shows the exact profile (curator
  confirms/prunes disorders + HPO; the variant comes from the classification), POST
  fires the worker task.
- **Results view:** table of `MMEMatchResult` rows — score, matched patient id, contact
  (name + `href`), and the returned phenotype/genes. Contacting the other curator is a
  manual human step (that is the whole point of MME).
- Use **DataTables** (`DatatableConfig` + `RichColumn`, `snpdb/views/datatable_view.py`)
  for the results grid, per project convention.

## 10. Privacy, consent & governance (non-code, blocking)

These gate go-live and are larger than the code. **The onboarding weight depends on
scope** — there are two paths, and the technical connection is bilateral either way:

**Path A — bilateral connection with a specific partner node** (lightest). The protocol
is decentralized: request an auth token from each node we want to query, issue them one
in return. Good for querying one or two friendly/collaborator nodes. No central
Steering Committee process, though the partner sets their own terms.

**Path B — become a recognised, discoverable MME Service** (needed to be listed and to
serve inbound reciprocal matches to the established network). The official *"I have a
database"* page: *"To become a Matchmaker Exchange Service, each new site must achieve
the Service Requirements set forth by the Steering Committee."* In practice the big
nodes (GeneMatcher, DECIPHER…) expect this before issuing a token. Involves:
- Contact **api@matchmakerexchange.org** to start.
- Meet/sign the three governance documents: **Informed Consent Policy**, **End User
  Agreement**, **Service Requirements** (2026 revision).
- Steering Committee review/approval.

Common to both paths:
- Exchange auth tokens with each remote node we query (fills `MME_NODES`) and issue our
  inbound token to peers.
- Consent is the classification's **`share_level = ShareLevel.PUBLIC`** — the same bar
  already used for ClinVar submission; no separate patient consent flag
  (`Patient.research_consent` is unrelated / internal-only). The MME `patient.id` and
  `label` we send must be non-PII (use the classification's opaque id).
- Keep `"test": true` until certified against a target node's test instance.
- Legal/ethics sign-off on serving inbound queries (we expose PUBLIC classifications'
  gene/variant + condition, and linked-patient HPO where present) — heavier for Path B.

A separate hand-off checklist for the lab head / research assistant lives in
`claude/matchmaker_exchange_onboarding.md`.

## 11. Rollout order

1. `mme/` app skeleton, models, migrations, settings/secrets block (`MME_ENABLED=False`).
2. Outbound profile serializer (§6) + unit tests with fake patients
   (`annotation/tests/test_data_fake_genes.py`, `snpdb/tests/utils/`). No network.
3. `MMEClient` + `submit()` task (§7); test against a node's **test** instance in
   test-mode.
4. Submit/confirm UI + results view (§9).
5. Inbound `/match` endpoint + simple matcher (§8); test with a peer's test query.
6. Governance/onboarding (§10) in parallel — the real critical path.

## 12. Testing

- `test_patient_profile.py` — condition-term **routing**: an HPO condition term lands in
  `features`, an OMIM term in `disorders` (`MIM:`), a MONDO term is crosswalked to OMIM
  via `as_omim`; a mixed-ontology condition splits across both slots; patient HPO merges
  into `features` and de-dups against condition HPO; classification with no linked
  patient still submits (condition terms + variant); coordinate 0-based conversion;
  "no features and no genomicFeatures" raises; `test` flag present until production
  submit enabled.
- `test_ontology_snake.py` — `MME_ONTOLOGY_SNAKE_EXACT=True` aliases a MONDO condition to
  an OMIM `disorders` entry; `=False` yields no `disorders` for a MONDO-only condition.
  `MME_ONTOLOGY_PHENOTYPE_EXPANSION=True` on a disease-only condition adds HPO `features`,
  each carrying a `_derivedFrom` id and a "(via …)" label suffix; `=False` (default) adds
  none. Observed patient HPO never carries `_derivedFrom`; the request envelope includes
  `disclaimer` when set.
- `test_eligibility.py` — `mme_eligible_classifications()` includes only
  `share_level=PUBLIC` + `is_last_published` + non-withdrawn; excludes `ALL_USERS`,
  `LAB`, unpublished.
- `test_client.py` — mock `requests.post`; header/content-type correctness; token sent;
  `raise_for_status` → `MMESubmissionStatus.ERROR` + admin notification; results
  persisted to `MMEMatchResult`.
- `test_inbound_match.py` — reject missing/incorrect `X-Auth-Token` (401); reject body
  with neither features nor genomicFeatures (400); only PUBLIC classifications returned;
  `MMEInboundQuery` audit row written. Extend `URLTestCase`
  (`library/django_utils/unittest_utils.py`) for URL status coverage.

## 13. Key existing references

| Concern | Existing code to copy / reuse |
| --- | --- |
| Outbound HTTP client shape | `library/oauth.py:ServerAuth` (copy pattern, static-token only) |
| Failure → admin notification | `library/log_utils.AdminNotificationBuilder` |
| Consent / share level | `classification/enums/classification_enums.py:243` (`ShareLevel.PUBLIC` = "3rd Party Databases") |
| Eligible classifications query | `ClassificationModification` + `is_last_published`/`share_level` (`classification/models/classification.py:2419`) |
| Condition terms (any ontology) | `Classification.condition_resolution_obj.terms` (`classification.py:637`, `ConditionResolved` `:234`) |
| MONDO → OMIM crosswalk (EXACT) | `OntologySnake.as_omim()` (`ontology/models/models_ontology.py:697`) |
| Cross-ontology enrichment (e.g. disease → HPO) | `OntologySnake.snake_from(term, to_ontology)` (`models_ontology.py:1100`); paths via `OntologyTermRelation` |
| Classification → patient → HPO | `Classification.sample` (`classification.py:469`) → `Sample.patient`; `PhenotypeDescription.get_ontology_term_ids()` (`models_phenotype_match.py:59`) |
| Variant → coords / HGVS / CAID | `snpdb/models/models_variant.py` (`Allele`, `VariantCoordinate`), `genes/hgvs/`, `snpdb/models/models_clingen_allele.py` |
| "Curated variant → external API" template | `classification/views/exports/classification_export_formatter_clinvar.py` |
| Inbound DRF endpoint | `ontology/views_rest.py` (`APIView` + `@extend_schema`) |
| Results grid | `snpdb/views/datatable_view.py` (`DatatableConfig`, `RichColumn`) |
| Secret tokens only | `variantgrid/settings/components/default_settings.py:215` (`get_secret`); other MME config is plain settings |

**Explicitly avoided dependencies:** `sync/sync_runner.py` / `SyncRunner`
(Shariant-specific) and everything under `sync/alissa/` (abandoned, slated for removal).
Patterns are copied into `mme/`, not imported.
