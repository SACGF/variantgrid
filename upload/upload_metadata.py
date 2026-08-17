"""
Per-upload metadata - facts about a file that its own contents do not reliably carry, supplied by
whatever submits it.

Stored as a blob on FileUpload.metadata rather than typed columns because the keys are per file type -
a VCF takes 'genome_build' and 'source', other types take neither. Each ImportTaskFactory declares the
keys it accepts via get_metadata_keys().

Validation happens at upload time, while the client is still connected, so a mistyped key is a 400
rather than a silently ignored value and a REQUIRES_USER_INPUT several pipeline stages later.
"""
import json
from typing import Optional

from patients.external_references import ExternalReference, ExternalReferenceError
from snpdb.models.models_genome import GenomeBuild
from upload.import_task_factories.import_task_factory import get_import_task_factories

GENOME_BUILD = "genome_build"
SOURCE = "source"
EXTRACTION = "extraction"                  # one extraction for every sample in the file
SAMPLE_EXTRACTIONS = "sample_extractions"  # VCF sample name -> extraction, for a multi-sample VCF

VCF_METADATA_KEYS = frozenset({GENOME_BUILD, SOURCE, EXTRACTION, SAMPLE_EXTRACTIONS})


class UploadMetadataError(ValueError):
    """ Client-supplied metadata we can reject before an import starts """


def _validate_genome_build(value) -> str:
    """ Clients should send a build's own name ('GRCh37') rather than an alias ('hg19'): an alias can
        collide with another build's name, and only the build itself says which patch you meant.
        Aliases still resolve - we store the canonical name, which is GenomeBuild's primary key. """
    try:
        return GenomeBuild.get_name_or_alias(str(value)).name
    except GenomeBuild.DoesNotExist as e:
        available = GenomeBuild.available_names_or_aliases()
        raise UploadMetadataError(f"Unknown '{GENOME_BUILD}': '{value}' - must be one of: {available}") from e
    except GenomeBuild.MultipleObjectsReturned as e:
        # Only where a deployment has both builds enabled, eg 'hg19' is a build name and GRCh37's alias
        raise UploadMetadataError(f"Ambiguous '{GENOME_BUILD}': '{value}' matches more than one enabled "
                                  f"build - use the build's own name") from e


def _validate_source(value) -> str:
    source = str(value).strip()
    if not source:
        raise UploadMetadataError(f"'{SOURCE}' must not be empty")
    return source


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


_VALIDATORS = {
    GENOME_BUILD: _validate_genome_build,
    SOURCE: _validate_source,
    EXTRACTION: _validate_extraction,
    SAMPLE_EXTRACTIONS: _validate_sample_extractions,
}


def get_metadata_keys_for_file_type(file_type) -> frozenset[str]:
    """ Keys accepted for a file type - union across the factories claiming it """
    keys = set()
    if file_type:
        for itf in get_import_task_factories():
            if itf.get_uploaded_file_type() == file_type:
                keys.update(itf.get_metadata_keys())
    return frozenset(keys)


def validate_upload_metadata(metadata: Optional[dict], allowed_keys) -> dict:
    """ Returns the cleaned metadata to store, raising UploadMetadataError on anything we can reject now """
    if not metadata:
        return {}
    if not isinstance(metadata, dict):
        raise UploadMetadataError(f"Upload metadata must be an object, got {type(metadata).__name__}")

    allowed_keys = frozenset(allowed_keys)
    if unknown_keys := set(metadata) - allowed_keys:
        unknown = ", ".join(sorted(unknown_keys))
        accepted = ", ".join(sorted(allowed_keys)) if allowed_keys else "(none for this file type)"
        raise UploadMetadataError(f"Unknown upload metadata key(s): {unknown}. Accepted keys: {accepted}")

    if EXTRACTION in metadata and SAMPLE_EXTRACTIONS in metadata:
        # One says the whole file is an extraction and the other says it is per-sample - guessing
        # which wins is the kind of silent decision upload-time validation exists to avoid
        raise UploadMetadataError(f"Send '{EXTRACTION}' or '{SAMPLE_EXTRACTIONS}', not both")

    return {key: _VALIDATORS[key](value) for key, value in metadata.items()}


def get_metadata_genome_build(file_upload) -> Optional[GenomeBuild]:
    """ The build the submitter declared, or None. Validated at upload so this will resolve """
    if file_upload and (build_name := (file_upload.metadata or {}).get(GENOME_BUILD)):
        # Validation stored the canonical name, which is the PK - don't re-run ambiguous alias resolution
        return GenomeBuild.objects.get(pk=build_name)
    return None


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


def get_metadata_source(file_upload) -> Optional[str]:
    """ The caller/software the submitter declared, matched against VCFSourceSettings.source_regex """
    if file_upload:
        return (file_upload.metadata or {}).get(SOURCE)
    return None
