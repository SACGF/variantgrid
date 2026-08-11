# TSO 500 Phase 4 — how identifiers cross the lab boundary

Design for Phase 4 of the TSO 500 order of work
([`tso500_overall_plan.md`](tso500_overall_plan.md) § "Phase 4"), covering
[#1707](https://github.com/SACGF/variantgrid/issues/1707) (patient/specimen/extraction API), the
seqauto extraction link, and [#1559](https://github.com/SACGF/variantgrid/issues/1559)
(specimen-level measures).

The overall plan is the spec for *what* and *why*; this is the *how*. Where the two disagree, the
overall plan wins and this file is wrong.

One subject in four steps:

1. **Resolve** — one helper that turns a client's identifier into a `Patient` / `Specimen` /
   `Extraction`, used by everything below (§1).
2. **Create** — #1707's API, so the client can accession before or alongside posting a run (§2).
3. **Name** — two ways for a VCF to say which extraction it belongs to: a seqauto link call (§3) and
   an upload metadata key (§4).
4. **Reconcile** — absorb the two arriving in either order (§5), then hang measures off the result (§6).

#431's decision governs the phase: **upload the files separately and join later**. So everything here
is assignment on top of the per-file upload that already exists — no bulk-import path, no run loader.

---

## What already exists (read before starting)

| Thing | Where |
|---|---|
| `Patient` / `Specimen` / `Extraction`, `ExternalPK`, `ExternallyManagedModel` | `patients/models.py:65,85,131,310,398` |
| `Sample.extraction`, `Sample.variants_type` | `snpdb/models/models_vcf.py:341,343` |
| `SequencingSample.extraction` (Phase 1, nothing sets it) | `seqauto/models/models_seqauto.py:286` |
| `SequencingSampleLookupSerializer`, `SampleSheetSerializer._create_sequencing_samples` | `seqauto/serializers/sequencing_serializers.py:206,250` |
| `create_backend_vcf_links` → `link_samples_and_vcfs_to_sequencing` | `upload/vcf/vcf_import.py:405,436` |
| `get_samples_by_sequencing_sample` (prefix match, VCF name fallback) | `seqauto/models/models_seqauto.py:1008` |
| `FileUpload.metadata` blob, key vocabulary, upload-time validation | `upload/upload_metadata.py`, `upload/models/models.py:70` |
| `GenomeBuildMismatchException` — the "declared contradicts detected" precedent | `upload/vcf/vcf_import.py:557` |
| `SimpleVCFImportInfo.add_message_count` — surfacing counts on the pipeline page | `upload/models/models.py:648` |
| Existing REST viewsets + router registration | `seqauto/views_rest.py`, `seqauto/urls.py:155` |
| Guardian assignment on patient create | `patients/import_records.py:184` |

`patients` has **no serializers and no `views_rest.py` today** — this phase creates both.

---

## §1 — Naming an externally managed object

The link call (§3), the metadata key (§4) and #1707's own payloads (§2) all quote the same kind of
identifier, and all three models are `ExternallyManagedModel`, so it is written once.

### The reference payload

A bare string is the **local reference** — `Extraction.reference_id`, `Specimen.reference_id`,
`Patient.patient_code`. Reaching an `ExternalPK` takes an object, and takes **both `code` and
`external_type`**: a code on its own does not identify anything, because `external_type` names the
identifier scheme it belongs to (a deployment carries several — HelixID, SAPOrderNumber — and
`external_pk_autocomplete_form_factory` at `patients/forms.py:166` is already built per type). Both
forms are JSON, so the metadata blob nests either.

```jsonc
"2600000001C"                                                   // the local reference
{"code": "H12345", "external_type": "HelixID"}                  // an ExternalPK - code AND type
{"code": "H12345", "external_type": "HelixID",
 "external_manager": "HELIX"}                                   // narrowed to one manager's namespace
{"reference_id": "2600000001C"}                                 // the bare form, written out
{"code": "H12345", "external_type": "HelixID",
 "reference_id": "2600000001C"}                                 // both - narrows to one object
```

A `code` without an `external_type` is rejected on shape, while the client is still connected, rather
than becoming a row nobody can resolve.

Nothing needs pinning to a VG-side vocabulary beyond that. `external_pk` is a OneToOne *from* each
model, so `Extraction.objects.filter(external_pk__code=...)` can only return extractions — the model
being queried already scopes the lookup. `external_manager` stays optional narrowing: the uniqueness
is `(code, external_type, external_manager)`, so on a deployment where two managers share a code and
type the pair alone is genuinely ambiguous, and that parks as Needs attention rather than picking one.

This is deliberately stricter than `search_external_pk` (`patients/signals/external_pk_search.py:19`),
which does filter on code alone. Search is a human typing something into a box and reading the
results; this is a machine asserting a link that nobody will re-read.

### The only per-model knowledge

Which field holds the local reference. `Specimen` and `Extraction` call it `reference_id`;
`Patient` calls it `patient_code`. A class attribute on the abstract base, overridden once:

```python
# patients/models.py
class ExternallyManagedModel(TimeStampedModel):
    external_pk = models.OneToOneField(ExternalPK, null=True, on_delete=CASCADE)
    # The field holding the local reference beside external_pk - @see patients.external_references
    LOCAL_REFERENCE_FIELD = "reference_id"


class Patient(GuardianPermissionsMixin, HasPhenotypeDescriptionMixin, ExternallyManagedModel, PreviewModelMixin):
    LOCAL_REFERENCE_FIELD = "patient_code"
```

### The resolver

```python
# patients/external_references.py
"""
How a client names a Patient / Specimen / Extraction, and how VG turns that into a row.

An externally managed object is found by its ExternalPK, which is unique on
(code, external_type, external_manager) - so a code needs its type to mean anything, and reaching one
takes both. An object that is not externally managed falls back to the local reference beside it,
which is what a bare string means. The caller supplies whichever it has, and each narrows the query
where present, so no deployment is forced into a canonical form it does not run.

Resolution has three outcomes rather than two. Nothing found is Pending, not an error: an extraction
legitimately arrives after the VCF that names it (@see patients.tasks.extraction_matching_tasks).
More than one found is Needs attention - reference_id is unique per parent rather than globally, so a
bare reference can match several rows, and that is ambiguity rather than a precedence rule.
"""
from dataclasses import dataclass
from typing import Optional

from django.contrib.auth.models import User
from django.db.models import Model, Q, QuerySet

from patients.models_enums import MatchStatus

REFERENCE_KEYS = frozenset({"code", "external_type", "external_manager", "reference_id"})


class ExternalReferenceError(ValueError):
    """ A reference we can reject on shape alone, before any lookup """


@dataclass(frozen=True)
class ExternalReference:
    reference_id: Optional[str] = None     # the local reference, whatever the model calls it
    code: Optional[str] = None             # ExternalPK.code - meaningless without external_type
    external_type: Optional[str] = None    # which identifier scheme the code belongs to
    external_manager: Optional[str] = None  # ExternalModelManager.name
    # Inferred from a sample name rather than asserted by a client (@see the sample-name fallback in
    # §5). Takes no part in the query - it only says how much to trust the row
    derived: bool = False

    @classmethod
    def from_data(cls, data) -> 'ExternalReference':
        """ Accepts the bare-string and object forms - shape only, no database access """
        if isinstance(data, str):
            data = {"reference_id": data}
        if not isinstance(data, dict):
            raise ExternalReferenceError(f"Reference must be a string or an object, got {type(data).__name__}")
        if unknown := set(data) - REFERENCE_KEYS:
            raise ExternalReferenceError(f"Unknown reference key(s): {', '.join(sorted(unknown))}. "
                                         f"Accepted keys: {', '.join(sorted(REFERENCE_KEYS))}")
        reference = cls(**{k: str(v) for k, v in data.items() if v is not None})
        if not (reference.reference_id or reference.code):
            raise ExternalReferenceError("Reference must supply 'reference_id' or 'code'")
        if bool(reference.code) != bool(reference.external_type):
            # ExternalPK is unique on (code, external_type, external_manager) - one without the other
            # names nothing, and is worth saying now rather than parking a row nobody can resolve
            raise ExternalReferenceError("'code' and 'external_type' must be supplied together")
        return reference

    @property
    def has_external_pk(self) -> bool:
        return bool(self.code)

    def as_json(self) -> dict:
        """ What gets parked on a row - round-trips back through from_data """
        return {k: v for k, v in self.__dict__.items() if v is not None}

    def q(self, model: type[Model]) -> Q:
        q = Q()
        if self.code:
            q &= Q(external_pk__code=self.code, external_pk__external_type=self.external_type)
        if self.external_manager:
            q &= Q(external_pk__external_manager=self.external_manager)
        if self.reference_id:
            q &= Q(**{model.LOCAL_REFERENCE_FIELD: self.reference_id})
        return q

    def halves(self) -> list['ExternalReference']:
        """ The external and local halves on their own, for telling 'neither matched' apart from
            'each matched something different' """
        parts = []
        if self.has_external_pk:
            parts.append(ExternalReference(code=self.code, external_type=self.external_type,
                                           external_manager=self.external_manager))
        if self.reference_id:
            parts.append(ExternalReference(reference_id=self.reference_id))
        return parts

    def __str__(self):
        return ", ".join(f"{k}={v}" for k, v in sorted(self.as_json().items()))


@dataclass
class ResolvedReference:
    reference: ExternalReference
    status: str                     # MatchStatus
    obj: Optional[Model] = None
    error: Optional[str] = None

    @property
    def matched(self) -> bool:
        return self.status == MatchStatus.MATCHED


def resolve_reference(model: type[Model], reference: ExternalReference,
                      user: Optional[User] = None) -> ResolvedReference:
    """ user scopes the search through filter_for_user - a client can only attach objects it can see,
        and one it cannot reads as Pending rather than Needs attention """
    qs: QuerySet = model.filter_for_user(user) if user else model.objects.all()
    matches = list(qs.filter(reference.q(model))[:2])
    if len(matches) == 1:
        return ResolvedReference(reference, MatchStatus.MATCHED, obj=matches[0])
    if len(matches) > 1:
        return ResolvedReference(reference, MatchStatus.NEEDS_ATTENTION,
                                 error=f"{reference} matches more than one {model.__name__}")

    halves = reference.halves()
    if len(halves) > 1:
        # An ExternalPK and a local reference that each resolve, but to different rows. The combined
        # query just returns nothing, which reads as Pending - but this will never settle on its own
        found = {h: qs.filter(h.q(model)).first() for h in halves}
        distinct = {obj.pk for obj in found.values() if obj}
        if len(distinct) > 1:
            described = "; ".join(f"{h} -> {obj}" for h, obj in found.items() if obj)
            return ResolvedReference(reference, MatchStatus.NEEDS_ATTENTION,
                                     error=f"{reference} names more than one {model.__name__}: {described}")
    return ResolvedReference(reference, MatchStatus.PENDING,
                             error=f"No {model.__name__} found for {reference}")
```

```python
# patients/models_enums.py
class MatchStatus(models.TextChoices):
    """ An identifier's journey from posted to resolved - independent feeds arrive in any order, so a
        reference we cannot resolve yet is parked rather than rejected """
    MATCHED = 'M', 'Matched'
    PENDING = 'P', 'Pending'
    NEEDS_ATTENTION = 'N', 'Needs attention'
```

### Where the search is scoped by user

`resolve_reference` takes an optional user and filters through `filter_for_user`. The API passes
`request.user`; the VCF import passes `vcf.user`; the reconcile task passes `sample.vcf.user`. An
extraction that exists but is not visible therefore parks as Pending and promotes to Needs attention
after the window (§5) — the same path as one that does not exist yet, which is the safe reading.

`SequencingSample` resolves **unscoped** (`user=None`). It has no user of its own, and seqauto only
runs on intranet deployments where records are global — so there is nothing to scope against and no
confidentiality boundary being crossed. Worth revisiting only if seqauto ever runs somewhere records
are not shared.

---

## §2 — #1707: the patient / specimen / extraction API

First, because everything after it quotes the identifiers this creates. Every route into a `Specimen`
today is a human one (the sample form, the patient CSV, the admin).

### Serializers

```python
# patients/serializers.py
from django.db import transaction
from rest_framework import serializers

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.external_references import ExternalReference, ExternalReferenceError, resolve_reference
from patients.models import Extraction, ExternalModelManager, ExternalPK, Patient, Specimen, SpecimenMeasure


class ExternalReferenceField(serializers.Field):
    """ The bare-string or object form of an identifier - @see patients.external_references """

    def to_internal_value(self, data) -> ExternalReference:
        try:
            return ExternalReference.from_data(data)
        except ExternalReferenceError as e:
            raise serializers.ValidationError(str(e)) from e

    def to_representation(self, value) -> dict:
        return value.as_json() if isinstance(value, ExternalReference) else value


class ExternalPKSerializer(serializers.Serializer):
    """ All three of code / external_type / external_manager are required, because that triple is what
        ExternalPK is unique on - anything less names a row we could not look up again. """
    code = serializers.CharField()
    external_type = serializers.CharField()
    external_manager = serializers.CharField()

    def validate_external_manager(self, value) -> ExternalModelManager:
        """ On an intranet deployment the set of managers is known, so an unfamiliar name is a typo
            rather than a new system - and silently creating one would decide can_modify (ie whether
            VG users may edit those records) as a side effect of a misspelling. A superuser creating
            one through the API is a deliberate act, so that is allowed either way """
        if manager := ExternalModelManager.objects.filter(name=value).first():
            return manager

        user = self.context["request"].user
        if settings.PATIENTS_API_EXTERNAL_MANAGER_CREATE_ADMIN_ONLY and not user.is_superuser:
            known = ", ".join(ExternalModelManager.objects.values_list("name", flat=True)) or "(none configured)"
            raise serializers.ValidationError(f"Unknown external_manager '{value}'. Known: {known}")
        return ExternalModelManager.objects.create(name=value)

    def get_or_create(self) -> ExternalPK:
        external_pk, _ = ExternalPK.objects.get_or_create(**self.validated_data)
        return external_pk
```

Each of the three model serializers is an **upsert keyed on the identifiers**, so a client that
re-posts a run gets the same rows rather than duplicates:

```python
class PatientSerializer(serializers.ModelSerializer):
    external_pk = ExternalPKSerializer(required=False)

    class Meta:
        model = Patient
        fields = ("id", "patient_code", "family_code", "first_name", "last_name", "date_of_birth",
                  "date_of_death", "sex", "affected", "external_pk")
        read_only_fields = ("id",)

    def create(self, validated_data):
        external_pk_data = validated_data.pop("external_pk", None)
        user = self.context["request"].user
        with transaction.atomic():
            patient = self._get_existing(validated_data, external_pk_data, user)
            if patient is None:
                patient = Patient.objects.create(**validated_data)
                assign_permission_to_user_and_groups(user, patient)
            else:
                patient.check_can_write(user)
                for field, value in validated_data.items():
                    setattr(patient, field, value)
            if external_pk_data:
                patient.external_pk = ExternalPKSerializer(data=external_pk_data).get_or_create()
            patient.save()
        return patient
```

`SpecimenSerializer` and `ExtractionSerializer` take the same shape, with their parent named by an
`ExternalReferenceField` instead of a PK:

```python
class SpecimenSerializer(serializers.ModelSerializer):
    patient = ExternalReferenceField()
    external_pk = ExternalPKSerializer(required=False)

    class Meta:
        model = Specimen
        fields = ("id", "patient", "reference_id", "description", "collected_by", "collection_date",
                  "received_date", "tissue_status", "external_pk")
        read_only_fields = ("id",)

    def validate_patient(self, reference: ExternalReference) -> Patient:
        """ A specimen has nowhere to live without its patient, so an unresolvable one is a 400 rather
            than a parked row - unlike a VCF naming an extraction, which has a row to park on """
        resolved = resolve_reference(Patient, reference, self.context["request"].user)
        if not resolved.matched:
            raise serializers.ValidationError(resolved.error)
        return resolved.obj
```

`tissue_status` serializes as its `TissueStatus` value (Phase 3's field — Reference / Affected /
Unknown, never null). The API is where the DNA and RNA arms of one block get told apart:
`nucleic_acid_source` is on `Extraction`, so both arms post against one specimen.

### Viewsets

```python
# patients/views_rest.py
class PatientViewSet(viewsets.ModelViewSet):
    serializer_class = PatientSerializer

    def get_queryset(self):
        return Patient.filter_for_user(self.request.user)


class SpecimenViewSet(viewsets.ModelViewSet):
    serializer_class = SpecimenSerializer

    def get_queryset(self):
        return Specimen.filter_for_user(self.request.user)


class ExtractionViewSet(viewsets.ModelViewSet):
    serializer_class = ExtractionSerializer

    def get_queryset(self):
        return Extraction.filter_for_user(self.request.user)

    def perform_create(self, serializer):
        super().perform_create(serializer)
        # New upstream data - re-run anything parked waiting for exactly this
        reconcile_pending_extractions.delay()
```

Guardian gives all three `filter_for_user` (Phase 3 made `Specimen` and `Extraction` delegate to
their patient), so no hand-rolled filters and no per-view permission classes — DRF's
`DEFAULT_PERMISSION_CLASSES = [IsAuthenticated]` is global.

### URLs

```python
# patients/urls.py
router = routers.DefaultRouter()
router.register(r'api/v1/patient', PatientViewSet, basename='api_patient')
router.register(r'api/v1/specimen', SpecimenViewSet, basename='api_specimen')
router.register(r'api/v1/extraction', ExtractionViewSet, basename='api_extraction')
router.register(r'api/v1/specimen_measure', SpecimenMeasureViewSet, basename='api_specimen_measure')

urlpatterns = [
    ...
    path('', include(router.urls), name='patients_apis'),
    path('api/v1/specimen_measure/bulk_create', SpecimenMeasureBulkCreateView.as_view(),
         name='api_specimen_measure_bulk_create'),
]
```

Add the new names to `shariantcommon.py`'s `URLS_NAME_REGISTER` overrides beside `view_specimen` /
`view_extraction` — Shariant has no patients menu and should not grow a patient API.

---

## §3 — Route 1: link the `SequencingSample` to the `Extraction`

By the time a file arrives the client has already made its seqauto calls and the §2 patient calls.
`Sample` ↔ `SequencingSample` links itself by path (`create_backend_vcf_links` →
`link_samples_and_vcfs_to_sequencing`); the only thing missing on that path is which extraction the
`SequencingSample` is.

Its own endpoint rather than fields on `SampleSheetSerializer`: most sample sheets have no extraction
to name, and most samples that have an extraction never had a sample sheet.

```python
# seqauto/serializers/sequencing_serializers.py
class SequencingSampleExtractionLinkSerializer(serializers.Serializer):
    """ One call per sequencing sample. All of that arm's VCFs inherit the extraction from it, since
        link_samples_and_vcfs_to_sequencing carries it down to every Sample it creates """
    sequencing_sample = SequencingSampleLookupSerializer()
    extraction = ExtractionReferenceField()

    def save(self, **kwargs):
        sequencing_sample = SequencingSampleLookupSerializer.get_object(self.validated_data["sequencing_sample"])
        resolved = resolve_reference(Extraction, self.validated_data["extraction"])
        sequencing_sample.apply_extraction_match(resolved)
        return sequencing_sample
```

```python
# seqauto/views_rest.py
class SequencingSampleExtractionLinkView(APIView):
    """ Name the extraction a sequencing sample was made from.

        The sequencing sample is the anchor: without it there is nothing to park a claim on, so an
        unknown one is a 400. An extraction that hasn't been created yet is a 202 - it is parked on
        the row and re-resolved by reconcile_pending_extractions. """

    @extend_schema(
        summary="Link a sequencing sample to the extraction it was made from",
        request=SequencingSampleExtractionLinkSerializer,
        responses=OpenApiTypes.OBJECT,
    )
    def post(self, request, *args, **kwargs):
        serializer = SequencingSampleExtractionLinkSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        sequencing_sample = serializer.save()
        response = {
            "sequencing_sample": str(sequencing_sample),
            "match_status": sequencing_sample.get_extraction_match_status_display(),
            "match_error": sequencing_sample.extraction_match_error,
            "extraction": str(sequencing_sample.extraction) if sequencing_sample.extraction else None,
        }
        status_code = status.HTTP_200_OK if sequencing_sample.extraction else status.HTTP_202_ACCEPTED
        return Response(response, status=status_code)
```

Registered as `path('api/v1/sequencing_sample/link_extraction', ...,
name='api_sequencing_sample_link_extraction')` beside the other bulk-create paths in
`seqauto/urls.py`.

### Carrying it down to `Sample`

`link_samples_and_vcfs_to_sequencing` already loops over the resolved pairs with a `modified_sample`
flag (`upload/vcf/vcf_import.py:482`):

```python
        for sequencing_sample, sample in samples_by_sequencing_sample.items():
            modified_sample = False
            # Set sample.variants_type from EnrichmentKit (eg to Mixed somatic for somatic panels)
            if ek := sequencing_sample.enrichment_kit:
                ...

            # One link call per sequencing sample reaches all of that arm's VCFs through here
            if sequencing_sample.extraction_id and not sample.extraction_id:
                sample.extraction = sequencing_sample.extraction
                sample.extraction_match_status = MatchStatus.MATCHED
                modified_sample = True
            elif sequencing_sample.extraction_reference and not sample.extraction_reference:
                # Still parked upstream - carry the claim so one reconcile pass settles both rows
                sample.extraction_reference = sequencing_sample.extraction_reference
                sample.extraction_match_status = sequencing_sample.extraction_match_status
                sample.extraction_match_date = timezone.now()
                modified_sample = True
```

### A re-sent sheet

`_create_sequencing_samples` builds fresh `SequencingSample` rows, so an existing link is carried
across from the previous current sheet's row with the same `sample_id` rather than making the client
re-link:

```python
# seqauto/serializers/sequencing_serializers.py
def carry_extractions_to_new_sample_sheet(previous: SampleSheet, current: SampleSheet):
    """ A re-sent sheet builds new SequencingSample rows, and the extraction is a property of the
        library rather than of the sheet that described it """
    current_by_sample_id = {ss.sample_id: ss for ss in current.sequencingsample_set.all()}
    for old in previous.sequencingsample_set.filter(Q(extraction__isnull=False) |
                                                    Q(extraction_reference__isnull=False)):
        new = current_by_sample_id.get(old.sample_id)
        if new and not (new.extraction_id or new.extraction_reference):
            new.extraction = old.extraction
            new.extraction_reference = old.extraction_reference
            new.extraction_match_status = old.extraction_match_status
            new.extraction_match_error = old.extraction_match_error
            new.extraction_match_date = old.extraction_match_date
            new.save()


class SampleSheetSerializer(serializers.ModelSerializer):
    def create(self, validated_data):
        sequencing_samples_data = validated_data.pop('sequencingsample_set')
        sequencing_run = validated_data["sequencing_run"]
        # Captured before set_as_current_sample_sheet moves the pointer
        previous_sheet = SampleSheet.objects.filter(
            sequencingruncurrentsamplesheet__sequencing_run=sequencing_run).first()
        sample_sheet, created = SampleSheet.objects.update_or_create(...)
        self._create_sequencing_samples(sample_sheet, sequencing_samples_data)
        if previous_sheet and previous_sheet != sample_sheet:
            carry_extractions_to_new_sample_sheet(previous_sheet, sample_sheet)
        sample_sheet.set_as_current_sample_sheet(sequencing_run, created)
        return sample_sheet
```

**Keep this on the API rather than the sample-sheet signal.** The SA Pathology app hooks
`sequencing_run_sample_sheet_created_signal` (`seqauto/models/models_seqauto.py:248`) to assign
patients from HelixIDs in the sheet, and the same hook could reach extractions. It stays where it is:
it only fires where seqauto runs, it is per-deployment glue rather than anything an API contract can
state, and route 2 has to exist regardless.

---

## §4 — Route 2: an `extraction` key in upload metadata

The general mechanism, and the only one available off the seqauto path — a deployment not running
seqauto, or a hand-uploaded file anywhere. Phase 2 built the blob and the per-file-type key
vocabulary, so this is two more keys.

```python
# upload/upload_metadata.py
GENOME_BUILD = "genome_build"
SOURCE = "source"
EXTRACTION = "extraction"                  # one extraction for every sample in the file
SAMPLE_EXTRACTIONS = "sample_extractions"  # VCF sample name -> extraction, for a multi-sample VCF

VCF_METADATA_KEYS = frozenset({GENOME_BUILD, SOURCE, EXTRACTION, SAMPLE_EXTRACTIONS})


def _as_json_value(key: str, value):
    """ Upload query params arrive as strings, so a JSON value comes through as its own text """
    if isinstance(value, str) and value.lstrip()[:1] in "{[":
        try:
            return json.loads(value)
        except json.JSONDecodeError as e:
            raise UploadMetadataError(f"'{key}' is not valid JSON: {e}") from e
    return value


def _validate_extraction(value):
    """ Shape only - the reference resolves at import.

        Phase 2 validates genome_build while the client is still connected precisely so a typo is a
        400; existence-checking an extraction here would instead reject the ordering race that
        patients.tasks.extraction_matching_tasks exists to absorb. """
    try:
        return ExternalReference.from_data(_as_json_value(EXTRACTION, value)).as_json()
    except ExternalReferenceError as e:
        raise UploadMetadataError(f"'{EXTRACTION}': {e}") from e


def _validate_sample_extractions(value):
    """ Keyed on VCF sample name, so nothing but sample names lives at this level """
    value = _as_json_value(SAMPLE_EXTRACTIONS, value)
    if not isinstance(value, dict) or not value:
        raise UploadMetadataError(f"'{SAMPLE_EXTRACTIONS}' must be a non-empty object keyed on "
                                  f"VCF sample name")
    try:
        return {name: ExternalReference.from_data(ref).as_json() for name, ref in value.items()}
    except ExternalReferenceError as e:
        raise UploadMetadataError(f"'{SAMPLE_EXTRACTIONS}': {e}") from e
```

Sending both is an `UploadMetadataError` — one says the whole file is an extraction and the other
says it is per-sample, and guessing which wins is the kind of silent decision Phase 2's upload-time
validation exists to avoid. `validate_upload_metadata` grows the check, since it is the only place
that sees the keys together:

```python
    if EXTRACTION in metadata and SAMPLE_EXTRACTIONS in metadata:
        raise UploadMetadataError(f"Send '{EXTRACTION}' or '{SAMPLE_EXTRACTIONS}', not both")
```

Single-sample files, which is every file in the TSO 500 set, just send the bare string:

```
POST /upload/api/v1/file?source=DRAGEN%20TSO500%20CNV&extraction=2600000001C
```

#### Why two keys rather than one

`Sample.extraction` is per-sample, and a VCF can carry several sample columns — a tumour/normal file
is two libraries off two extractions in one VCF. So metadata has to be able to say either "this whole
file is extraction X" or "column TUMOUR is X, column NORMAL is Y".

Those are two different types, and a single key carrying both would be ambiguous, since each form is
a JSON object with no structural tell:

```jsonc
{"code": "H12345", "external_type": "HelixID"}      // one reference for the file
{"TUMOUR": "2600000001C", "NORMAL": "2600000002C"}  // a per-sample map
```

Telling those apart means asking whether the keys look like reference keys, which quietly breaks the
day a VCF has a sample column named `code` — failing silently and wrongly by attaching the whole file
to an extraction nobody named. It is the same hazard as sample names sharing a namespace with the
metadata vocabulary at the top level, where a sample called `genome_build` would collide with a key
that already means something.

Two keys make each one a single type, and confine sample names to a level where nothing else lives:

```jsonc
{"extraction": "2600000001C"}
{"sample_extractions": {"TUMOUR": "2600000001C",
                        "NORMAL": {"code": "H12345", "external_type": "HelixID"}}}
```

No nesting, no reserved word inside a value, and an unknown top-level key is already a 400 from
Phase 2's vocabulary check.

### Applying it at import

`create_vcf_from_vcf` creates samples (`configure_vcf_from_header`), then does the seqauto linking.
A third step runs after both, so it can see what route 1 already set.

**Not in the sample-creation loop itself**, for two reasons. `configure_vcf_from_header` is also
called on VCF reload (its own docstring says so), so putting `extraction` in the
`Sample.objects.update_or_create` defaults would overwrite a hand-assigned extraction every time
someone reloads a VCF. And the seqauto carry-down runs *after* sample creation, so a disagreement
between the two routes isn't visible from inside that loop. Reading the metadata once, after both,
keeps one function deciding — and `apply_extraction_match` leaves an already-set FK alone, so the
reload case is safe either way round.

```python
# upload/vcf/vcf_import.py
def create_vcf_from_vcf(upload_step, vcf_reader) -> VCF:
    ...
    backend_vcf = create_backend_vcf_links(uploaded_vcf)
    if backend_vcf:
        logging.info("Handle backend VCF")
        link_samples_and_vcfs_to_sequencing(backend_vcf, upload_step=upload_step)
        backend_vcf_import_start_signal.send(sender=os.path.basename(__file__), backend_vcf=backend_vcf)

    assign_sample_extractions(vcf, upload_step)
    return vcf


def assign_sample_extractions(vcf: VCF, upload_step=None):
    """ Route 2 - the extraction upload metadata keys, which need no seqauto records at all.

        Where a file has both this and a seqauto link and they name different extractions we fail the
        import, the same rule a declared genome_build the header contradicts gets. """
    samples = list(vcf.sample_set.all())
    declared = get_metadata_extractions(_get_file_upload(vcf), [s.vcf_sample_name for s in samples])
    if unknown := set(declared) - {s.vcf_sample_name for s in samples}:
        raise ExtractionMismatchException(f"'{SAMPLE_EXTRACTIONS}' names sample(s) not in this VCF: "
                                          f"{', '.join(sorted(unknown))}")

    for sample in samples:
        reference = declared.get(sample.vcf_sample_name) or _derive_extraction_reference(sample.name)
        if reference is None:
            continue
        resolved = resolve_reference(Extraction, reference, vcf.user)
        if sample.extraction_id and resolved.matched and resolved.obj != sample.extraction:
            msg = f"Sample '{sample.vcf_sample_name}': declared extraction '{reference}' disagrees " \
                  f"with '{sample.extraction}' linked through the sequencing sample"
            raise ExtractionMismatchException(msg)
        sample.apply_extraction_match(resolved)
        if upload_step and not resolved.matched:
            SimpleVCFImportInfo.add_message_count(1, resolved.error, upload_step)
```

`get_metadata_extractions(file_upload)` returns `{vcf_sample_name: ExternalReference}` from whichever
key was sent — `extraction` expands to every sample in the file (in practice one), `sample_extractions`
maps by name and is left as-is:

```python
# upload/upload_metadata.py
def get_metadata_extractions(file_upload, sample_names) -> dict[str, ExternalReference]:
    """ Keyed on VCF sample name either way, so the caller has one shape to handle.

        A name in sample_extractions that the VCF doesn't have is the client's mistake to see, so it
        comes back too rather than being dropped - assign_sample_extractions reports it. """
    metadata = (file_upload.metadata or {}) if file_upload else {}
    if per_sample := metadata.get(SAMPLE_EXTRACTIONS):
        return {name: ExternalReference.from_data(ref) for name, ref in per_sample.items()}
    if declared := metadata.get(EXTRACTION):
        reference = ExternalReference.from_data(declared)
        return {name: reference for name in sample_names}
    return {}
```

`ExtractionMismatchException` sits beside `GenomeBuildMismatchException` and propagates the same way,
failing the upload step. A `sample_extractions` name matching no VCF sample fails the import too — it
means the client is naming a sample that isn't in the file it uploaded, which no later arrival fixes.

Parked references show on the upload pipeline page through `SimpleVCFImportInfo`, the way Phase 2's
copy-neutral skips do — so "this VCF is not attached to anything yet" is visible at import rather
than only on the sample page.

---

## §5 — Reconciling what arrives out of order

Mocha hit this and the answer is worth copying wholesale (SACGF/mocha#126, commits `ed1d0df`,
`982022e`). Here the parked claim always has a row to live on, so it needs no new table.

### The columns

Both `Sample` and `SequencingSample` already carry `extraction`, so the state goes beside it. One
abstract mixin so the two models, the reconcile task and the templates share the behaviour:

```python
# patients/models.py
class ExtractionMatchMixin(models.Model):
    """ A claim about which Extraction a row belongs to, which may not be resolvable yet.

        Neither Sample nor SequencingSample is a TimeStampedModel, so extraction_match_date carries
        when the claim was parked - which is what promotes a stale Pending to Needs attention. """
    extraction = models.ForeignKey(Extraction, null=True, blank=True, on_delete=SET_NULL)
    extraction_reference = models.JSONField(null=True, blank=True)
    extraction_match_status = models.CharField(max_length=1, choices=MatchStatus.choices,
                                               null=True, blank=True)
    extraction_match_error = models.TextField(null=True, blank=True)
    extraction_match_date = models.DateTimeField(null=True, blank=True)

    class Meta:
        abstract = True

    def apply_extraction_match(self, resolved: 'ResolvedReference', save=True) -> bool:
        """ Leaves a confirmed link alone, so a settled row never flaps back """
        if self.extraction_id:
            return False
        self.extraction_reference = resolved.reference.as_json()
        self.extraction_match_status = resolved.status
        self.extraction_match_error = resolved.error
        self.extraction_match_date = timezone.now()
        if resolved.matched:
            self.extraction = resolved.obj
        if save:
            self.save()
        return True
```

`Sample` and `SequencingSample` swap their `extraction` field declaration for the mixin:

```python
class Sample(GuardianPermissionsMixin, SortByPKMixin, PreviewModelMixin, ExtractionMatchMixin, models.Model):
```

The `extraction` field itself is identical, so the migration only adds the four new columns.

**Idempotency is keyed on the posted identifiers, never on the nullable FK.** Nothing here gains a
`unique_together` spanning `extraction` — Postgres treats nulls as distinct, so a parked resend would
duplicate. The sequencing-sample link upserts on the `SequencingSample` itself, and route 2 upserts
on `(vcf, vcf_sample_name)` as sample creation already does.

### The task

```python
# patients/tasks/extraction_matching_tasks.py
"""
Independent feeds arrive in any order: the client may post a run before it accessions the specimen, or
accession it a day later. A row whose extraction cannot be resolved yet is parked rather than rejected,
and this re-resolves the parked ones - on a schedule, and again whenever new extractions land.
"""
import celery
from django.conf import settings
from django.db.models import Q
from django.utils import timezone

from patients.external_references import ExternalReference, resolve_reference
from patients.models import Extraction, MatchStatus
from seqauto.models import SequencingSample
from snpdb.models import Sample

PENDING_STATES = [MatchStatus.PENDING, MatchStatus.NEEDS_ATTENTION]


@celery.shared_task(queue='db_workers')
def reconcile_pending_extractions() -> dict:
    counts = {"matched": 0, "still_pending": 0, "needs_attention": 0, "from_sequencing_sample": 0}

    for model in (SequencingSample, Sample):
        qs = model.objects.filter(extraction__isnull=True,
                                  extraction_reference__isnull=False,
                                  extraction_match_status__in=PENDING_STATES)
        for row in qs.iterator():
            user = getattr(getattr(row, "vcf", None), "user", None)
            reference = ExternalReference.from_data(row.extraction_reference)
            resolved = resolve_reference(Extraction, reference, user)
            if resolved.status == MatchStatus.PENDING and _past_pending_window(row):
                # Past the window this is a real mismatch rather than the load race, and wants a human
                resolved.status = MatchStatus.NEEDS_ATTENTION
            row.apply_extraction_match(resolved)
            counts[_count_key(resolved.status)] += 1

    # Route 1 arriving after the VCF: the link call set SequencingSample.extraction, but
    # link_samples_and_vcfs_to_sequencing had already run and had nothing to carry down
    unlinked = Sample.objects.filter(extraction__isnull=True,
                                     samplefromsequencingsample__sequencing_sample__extraction__isnull=False)
    for sample in unlinked.select_related("samplefromsequencingsample__sequencing_sample").iterator():
        sample.extraction = sample.samplefromsequencingsample.sequencing_sample.extraction
        sample.extraction_match_status = MatchStatus.MATCHED
        sample.extraction_match_date = timezone.now()
        sample.save()
        counts["from_sequencing_sample"] += 1

    return counts


def _past_pending_window(row) -> bool:
    parked = row.extraction_match_date
    if parked is None:
        return False
    days = settings.PATIENT_EXTRACTION_MATCH_PENDING_DAYS
    return timezone.now() - parked > timezone.timedelta(days=days)
```

The second half is the part that is easy to miss: a link call made *after* the VCF imported has no
path to `Sample` at all, because `link_samples_and_vcfs_to_sequencing` runs once, at import.

### Wiring

```python
# variantgrid/settings/components/celery_settings.py
CELERY_IMPORTS = (
    ...
    'patients.tasks.extraction_matching_tasks',
)
CELERY_TASK_ROUTES = {
    ...
    'patients.tasks.extraction_matching_tasks.reconcile_pending_extractions': DB_WORKERS,
}

# variantgrid/celery.py
app.conf.beat_schedule['reconcile-pending-extractions'] = {
    'task': 'patients.tasks.extraction_matching_tasks.reconcile_pending_extractions',
    'schedule': HOUR_SECS,
}

# variantgrid/settings/components/default_settings.py
# How long a parked extraction reference stays Pending before it needs a human - past this it is a
# real mismatch rather than feeds arriving out of order. Same window Mocha settled on
PATIENT_EXTRACTION_MATCH_PENDING_DAYS = 3

# An external_manager the API doesn't recognise is a typo on an intranet deployment, where the set of
# tracking systems is known - so only a superuser creates one via the API. A public server taking
# records from systems it has never seen would set this False
PATIENTS_API_EXTERNAL_MANAGER_CREATE_ADMIN_ONLY = True
```

Fired on schedule and again from `ExtractionViewSet.perform_create` (§2), which is when new upstream
data actually lands.

### The sample-name fallback

Everything above assumes something upstream quotes an identifier — a seqauto link call, or an
`extraction` key on the upload. A deployment with no tracking system posts neither, and for those the
sample name is the only identifier that exists. The TSO 500 names carry it:
`ExampleSample_DNA_2600000001C` contains `2600000001C`, which is exactly
`Extraction.reference_id` — a ten-digit specimen accession plus the container suffix naming the arm.

This is a fallback, not the mechanism. It is guesswork anywhere a client could have told us, so it is
off by default and configured per deployment:

```python
# variantgrid/settings/components/default_settings.py
# Deployments with no tracking system to quote identifiers: a regex read against a VCF sample name,
# whose 'extraction' named group is the extraction's reference_id. Only consulted where nothing was
# posted, so it can never override a client
PATIENT_EXTRACTION_SAMPLE_NAME_REGEX = None  # eg r"(?P<extraction>\d{10}[A-Z])$"
```

**It runs at import, in `assign_sample_extractions`, and only where nothing was posted.** Three
things follow from that, which is why the earlier "import or reconcile?" question dissolves:

- *It cannot fight a real reference.* The guard is `not sample.extraction_reference` — a posted
  reference has already been parked by the time the fallback is reached, so there is nothing to give
  first refusal to.
- *Late arrival is already handled.* The derived reference is parked exactly like a posted one, so if
  the extraction is created afterwards, `reconcile_pending_extractions` picks it up with no extra
  code. Running it at import therefore costs nothing and attaches the sample immediately in the
  common case where the extraction already exists.
- *It never creates anything.* Deriving a reference is a guess about naming; creating a `Specimen` or
  `Extraction` from it would be inventing patient records. An unresolvable derived reference parks
  and eventually asks a human, the same as any other.

The derived reference is marked as such, so anyone looking at a Needs attention row can tell whether
the lab asserted the link or VG inferred it. `ExternalReference` carries a flag that takes no part in
the query and round-trips through the parked JSON:

```python
@dataclass(frozen=True)
class ExternalReference:
    ...
    derived: bool = False   # inferred from a sample name rather than asserted by a client
```

```python
# upload/vcf/vcf_import.py - consulted per sample in assign_sample_extractions, only where the
# metadata named nothing for that sample
def _derive_extraction_reference(sample_name: str) -> Optional[ExternalReference]:
    """ Only for deployments with nothing upstream to quote an identifier - @see the setting """
    if pattern := settings.PATIENT_EXTRACTION_SAMPLE_NAME_REGEX:
        if m := re.search(pattern, sample_name):
            return ExternalReference(reference_id=m.group("extraction"), derived=True)
    return None
```

One limitation worth naming: `Extraction.reference_id` is unique per specimen rather than globally,
so a derived reference can in principle match rows under two specimens. That parks as Needs attention
like any other ambiguity. In practice a TSO 500 reference embeds its specimen (`2600000001C` starts
with `2600000001`), so it does not arise — a deployment whose naming does not carry the specimen
would need the regex to yield a parent too, which is not built here.

---

## §6 — #1559: specimen-level measures

Small, and it proves Phase 1's identifiers work across the LIMS boundary. **These arrive by API from
the lab client, not by VG parsing pipeline output** — the measures are scattered across the run (TMB
summary in one place, MSI in a separate JSON, GIS with tumour fraction and ploidy in another written
by a different tool) and none of that structure is worth VG knowing about. When Illumina moves a file
between releases it becomes a client change, not a VG release.

### The model

```python
# patients/models_enums.py
class SpecimenMeasureType(models.TextChoices):
    """ Vendor-neutral - other panels push the same shape """
    TMB = 'T', 'Tumour mutational burden'
    MSI = 'M', 'Microsatellite instability'
    GIS = 'G', 'Genomic instability score'
    TUMOUR_FRACTION = 'F', 'Tumour fraction'
    PLOIDY = 'P', 'Ploidy'
```

```python
# patients/models.py
class SpecimenMeasure(TimeStampedModel):
    """ A scalar measured on the material rather than on any one variant - TMB, MSI, GIS.

        Both the score and the call are stored. HRD gives a genomic instability score and no
        positive/negative call: the threshold that turns it into one is lab policy, not vendor output,
        so without the raw value plus the threshold applied and by whom a later re-interpretation is
        unreconstructable. Mirrors somatic:tmb_value beside somatic:tmb_status in the evidence keys. """
    specimen = models.ForeignKey(Specimen, on_delete=CASCADE)
    # Which arm produced the number - enrichment, since the measure describes the specimen
    extraction = models.ForeignKey(Extraction, null=True, blank=True, on_delete=SET_NULL)
    measure_type = models.CharField(max_length=1, choices=SpecimenMeasureType.choices)
    value = models.FloatField(null=True, blank=True)
    unit = models.TextField(null=True, blank=True)      # eg 'mut/Mb', '%'
    call = models.TextField(null=True, blank=True)      # the lab's call - 'High', 'Stable'
    threshold = models.TextField(null=True, blank=True)  # the threshold that produced the call
    threshold_source = models.TextField(null=True, blank=True)  # whose policy set it
    method = models.TextField(blank=True)               # tool and version the client transcribed from
    source_payload = models.JSONField(default=dict, blank=True)  # raw, so 'which file' stays answerable
    measured_date = models.DateTimeField(null=True, blank=True)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)

    class Meta:
        # One current value per measure - the report pulls these together and wants a single MSI, so a
        # resend replaces rather than accumulating. Not keyed on the nullable extraction FK: Postgres
        # treats nulls as distinct, so that would let a resend duplicate instead of updating
        unique_together = ("specimen", "measure_type")

    @classmethod
    def get_permission_class(cls):
        return Patient

    def get_permission_object(self):
        return self.specimen.patient
```

Four things held to, from the issue:

- **The client transcribes, it never computes.** Lifting a value out of vendor output is fine;
  deriving one moves a clinical calculation into an unversioned client nobody reviews. The API
  accepts what it is given and stores `source_payload` beside it.
- **Raw payload beside the parsed fields**, so "which file and which pipeline version produced this
  number" is answerable at report time. `HelixNGSOrder.raw_data` in Mocha is the existing pattern.
- **Keyed on the specimen/extraction external identifiers from #1704**, not a vendor pair-ID string.
- **Vendor-neutral schema** rather than a TSO500-shaped table.

Per-gene absolute and minor copy number are **not** in scope — per-gene rather than per-specimen, and
they arrive as a VCF, so they are ingestion. They land whenever Phase 0's request produces
`abcn_annotated.vcf`.

### The API

```python
# patients/serializers.py
class SpecimenMeasureSerializer(serializers.ModelSerializer):
    specimen = ExternalReferenceField()
    extraction = ExternalReferenceField(required=False)

    class Meta:
        model = SpecimenMeasure
        fields = ("id", "specimen", "extraction", "measure_type", "value", "unit", "call", "threshold",
                  "threshold_source", "method", "source_payload", "measured_date")
        read_only_fields = ("id",)

    def create(self, validated_data):
        """ A measure has no row to live on without its specimen, so an unresolvable one is a 400 -
            unlike a VCF naming an extraction, which parks on the sample it already created.

            A re-post replaces: there is one current MSI for a specimen, and the report wants that
            rather than a history to choose from. TimeStampedModel.modified records when it changed """
        ...
        measure, _ = SpecimenMeasure.objects.update_or_create(
            specimen=specimen, measure_type=validated_data["measure_type"], defaults=defaults)
        return measure
```

`SpecimenMeasureBulkCreateView` takes `{"specimen": <reference>, "measures": [...]}` so a run's TMB,
MSI, GIS, tumour fraction and ploidy post in one call, following
`SequencingFilesBulkCreateView`'s shape (`seqauto/views_rest.py:135`).

### Showing them

`view_specimen` gains a measures table under the extractions one:

```python
# patients/views.py
    context = {"specimen": specimen,
               "form": form,
               "extractions": specimen.extraction_set.order_by("pk").prefetch_related(visible_samples),
               "measures": specimen.specimenmeasure_set.order_by("measure_type", "-measured_date"),
               "has_write_permission": has_write_permission}
```

```html
{# patients/templates/patients/view_specimen.html #}
<hr/>
<h4>Measures</h4>
{% if measures %}
    <table class="table">
        <thead><tr><th>Measure<th>Value<th>Call<th>Method<th>Measured</tr></thead>
        <tbody>
        {% for measure in measures %}
            <tr>
                <td>{{ measure.get_measure_type_display }}
                <td>{{ measure.value|default_if_none:"-" }} {{ measure.unit|default:"" }}
                <td>{{ measure.call|default:"-" }}
                <td>{{ measure.method|default:"-" }}
                <td>{{ measure.measured_date|default_if_none:"-" }}
            </tr>
        {% endfor %}
        </tbody>
    </table>
{% else %}
    <p>No measures recorded for this specimen.</p>
{% endif %}
```

---

## §7 — Surfacing the match state

`status` and `match_error` read-only on the serializers, and unmatched rows findable:

- **`view_extraction` / `view_specimen`** — unchanged for matched rows; nothing new needed, since a
  matched sample already appears.
- **`view_sample`** — where `extraction` is null and `extraction_reference` is set, show the parked
  reference, the status and the error rather than an empty Extraction row, so the person looking at
  "why is this sample not attached" gets the answer on the page.
- **Admin** — `list_display` and `list_filter` on `extraction_match_status` for `Sample` and
  `SequencingSample`, which is the cheapest "show me everything needing attention" for now.
- **Health check** — the count of unmatched rows joins the existing overall stats, so it shows up in
  the daily server-status notification without anyone having to go looking:

```python
# patients/signals/extraction_match_health_check.py
@receiver(signal=health_check_overall_stats_signal)
def extraction_match_health_check(sender, **kwargs):
    """ Parked extraction references are only a problem once they stop clearing on their own, so
        report Needs attention as the number and Pending as context """
    counts = Counter()
    for model in (Sample, SequencingSample):
        qs = model.objects.filter(extraction__isnull=True, extraction_reference__isnull=False)
        counts.update(qs.values_list("extraction_match_status", flat=True))

    needs_attention = counts[MatchStatus.NEEDS_ATTENTION]
    if not (needs_attention or counts[MatchStatus.PENDING]):
        return None
    return HealthCheckTotalAmount(
        emoji=":test_tube:",
        amount=needs_attention,
        name="unmatched extractions",
        extra=f"{counts[MatchStatus.PENDING]} still pending")
```

Registered in `patients/apps.py:12` beside the search receivers. It reports nothing when everything
is matched, so a healthy deployment gains no line.
- **Sample grid** — the extraction column and its link are Phase 6's grid pass
  (`snpdb/grids.py:127`, already showing `extraction__specimen__reference_id`); this phase leaves the
  data ready rather than reworking the columns.

---

## §8 — Migrations

| App | Contents |
|---|---|
| `patients/0017` | `SpecimenMeasure`; no field changes to `Specimen`/`Extraction` |
| `snpdb/0205` | `Sample.extraction_reference`, `_match_status`, `_match_error`, `_match_date` |
| `seqauto/0047` | the same four on `SequencingSample` |

Moving `extraction` onto `ExtractionMatchMixin` is a no-op — abstract base fields are inlined into
the concrete model, so `makemigrations` sees only the four additions. Confirm with
`manage.py makemigrations --check --dry-run` before writing anything by hand.

No `ManualOperation` — nothing existing needs backfilling, since every row's `extraction` is either
already set or was never claimed.

---

## §9 — Tests

New `patients/tests/test_external_references.py`:

- A bare string resolves against the local reference, and does **not** match a row whose `ExternalPK`
  code happens to equal it.
- `code` without `external_type` (and vice versa) is an `ExternalReferenceError` on shape, before any
  query runs.
- `code` + `external_type` resolves; the same code under a different type does not.
- `external_manager` narrows where two managers share a code and type; without it that pair is Needs
  attention.
- `code`/`external_type` + `reference_id` naming different objects → Needs attention with both named,
  not Pending and not a precedence winner.
- A bare `reference_id` matching two patients' specimens → Needs attention (`reference_id` is unique
  per parent, not globally).
- Nothing found → Pending, with the reference in the error.
- `Patient` resolves on `patient_code` where the other two resolve on `reference_id`.

New `patients/tests/test_patient_api.py` (`APITestCase`):

- Creating patient → specimen → extraction by external identifiers, then re-posting the same payload,
  yields one row of each.
- A specimen naming an unknown patient is a 400.
- A user cannot see or attach another user's patient.
- `tissue_status` round-trips.
- An unknown `external_manager` is a 400 naming the configured ones; a superuser posting the same
  payload creates it; with `PATIENTS_API_EXTERNAL_MANAGER_CREATE_ADMIN_ONLY = False` so does anyone.

New `patients/tests/test_extraction_matching.py`:

- **Order A** — extraction exists, then the VCF imports with `extraction` metadata → `Sample.extraction`
  set, status Matched.
- **Order B** — the VCF imports first → status Pending, no FK; the extraction is then created →
  `reconcile_pending_extractions` matches it.
- **Order C** — the link call arrives after the VCF imported → the reconcile task's
  sequencing-sample pass carries it down to `Sample`.
- A matched row is left alone by a later reconcile pass (`apply_extraction_match` returns False).
- A Pending row older than `PATIENT_EXTRACTION_MATCH_PENDING_DAYS` promotes to Needs attention.
- Ambiguous reference parks as Needs attention and is not promoted or re-matched.
- With `PATIENT_EXTRACTION_SAMPLE_NAME_REGEX` set, a VCF sample named
  `ExampleSample_DNA_2600000001C` and no posted reference derives `2600000001C` and matches; the
  parked reference records `derived`.
- The fallback is not consulted where a reference was posted, and creates no `Specimen` or
  `Extraction` when it resolves to nothing.
- The health check reports nothing when every row is matched, and counts Needs attention when one is
  not.

New `seqauto/tests/test_extraction_link.py`:

- Link call with an unknown sequencing sample → 400; with an unknown extraction → 202 and a parked row.
- One link call reaches all three of the DNA arm's VCF samples through
  `link_samples_and_vcfs_to_sequencing`.
- A re-sent sample sheet carries the link across to the new `SequencingSample` rows.

Extend `upload/tests/test_api.py` / `upload/tests/vcf/test_vcf_processors.py`:

- `extraction` accepted as a bare string and as an object; `sample_extractions` as a map; a malformed
  one is a 400 at upload; an unknown *extraction* is **not**, since that is the ordering race.
- Sending both `extraction` and `sample_extractions` is a 400 at upload.
- A VCF sample named `code` or `genome_build` under `sample_extractions` resolves to that sample and
  nothing else — sample names share a namespace with nothing.
- `sample_extractions` naming a sample the VCF doesn't have fails the import.
- A file whose metadata and seqauto link name different extractions fails the import.

New `patients/tests/test_specimen_measures.py`:

- Bulk post of TMB/MSI/GIS against a specimen reference.
- Re-posting a specimen's MSI with a new value replaces it — one row, the new value, and `modified`
  moved on.
- A measure naming an unknown specimen is a 400.
- `source_payload` round-trips.

`patients/tests/test_urls.py` gains the new API names so `URLS_NAME_REGISTER` coverage stays honest.

---

## §10 — Done when

The whole test run lands against one specimen with both arms attached:

- Patient, specimen and both extractions (`2600000001C` DNA, `2600000001B` RNA) created over the API.
- The DNA arm's three VCFs (`hard-filtered`, `cnv`, `_DragenExonCNV`) reach `Sample.extraction`
  through **one** seqauto link call.
- A hand-uploaded file reaches `Sample.extraction` through upload metadata alone, with no seqauto
  records involved.
- A VCF uploaded *before* its extraction exists parks, and attaches itself once the extraction is
  created.
- TMB/MSI/GIS posted against that specimen show on `view_specimen`.

Then raise the client-side issue in `SACGF/variantgrid_api`.

---

## Settled during review

Nothing outstanding. Recorded so the reasoning survives into implementation:

| Question | Decision |
|---|---|
| Can a bare identifier reach an `ExternalPK`? | No — a bare string is the local reference. An `ExternalPK` takes `code` **and** `external_type`, together or neither (§1) |
| One key or two for per-sample extractions? | Two — `extraction` and `sample_extractions`, so each key has one type and sample names share a namespace with nothing (§4) |
| `SpecimenMeasure` uniqueness | `(specimen, measure_type)` — the report wants one current MSI, so a resend replaces (§6) |
| Unknown `external_manager` on the API | 400 unless the poster is a superuser, under `PATIENTS_API_EXTERNAL_MANAGER_CREATE_ADMIN_ONLY = True`; a public server flips it (§2) |
| Permission scoping for `SequencingSample` | Unscoped — seqauto only runs on intranet deployments where records are global (§1) |
| Pending window | 3 days, the window Mocha settled on (§5) |
| Surfacing Needs attention | A health check on the overall stats signal, silent when everything matches (§7) |
| Where the sample-name fallback runs | At import, guarded on nothing having been posted — so it cannot override a client, and late arrivals reconcile for free (§5) |
