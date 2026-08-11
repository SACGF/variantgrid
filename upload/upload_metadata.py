"""
Per-upload metadata - facts about a file that its own contents do not reliably carry, supplied by
whatever submits it.

Stored as a blob on FileUpload.metadata rather than typed columns because the keys are per file type -
a VCF takes 'genome_build' and 'source', other types take neither. Each ImportTaskFactory declares the
keys it accepts via get_metadata_keys().

Validation happens at upload time, while the client is still connected, so a mistyped key is a 400
rather than a silently ignored value and a REQUIRES_USER_INPUT several pipeline stages later.
"""
from typing import Optional

from snpdb.models.models_genome import GenomeBuild
from upload.import_task_factories.import_task_factory import get_import_task_factories

GENOME_BUILD = "genome_build"
SOURCE = "source"

VCF_METADATA_KEYS = frozenset({GENOME_BUILD, SOURCE})


class UploadMetadataError(ValueError):
    """ Client-supplied metadata we can reject before an import starts """


def _validate_genome_build(value) -> str:
    """ Stored as the build name - GenomeBuild's primary key is its name, so this is nearly an FK """
    try:
        return GenomeBuild.get_name_or_alias(str(value)).name
    except GenomeBuild.DoesNotExist as e:
        available = GenomeBuild.available_names_or_aliases()
        raise UploadMetadataError(f"Unknown '{GENOME_BUILD}': '{value}' - must be one of: {available}") from e
    except GenomeBuild.MultipleObjectsReturned as e:
        # e.g. 'hg19' is both a build name and GRCh37's alias - the submitter has to say which
        raise UploadMetadataError(f"Ambiguous '{GENOME_BUILD}': '{value}' matches more than one build - "
                                  f"use the build's own name") from e


def _validate_source(value) -> str:
    source = str(value).strip()
    if not source:
        raise UploadMetadataError(f"'{SOURCE}' must not be empty")
    return source


_VALIDATORS = {
    GENOME_BUILD: _validate_genome_build,
    SOURCE: _validate_source,
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

    return {key: _VALIDATORS[key](value) for key, value in metadata.items()}


def get_metadata_genome_build(file_upload) -> Optional[GenomeBuild]:
    """ The build the submitter declared, or None. Validated at upload so this will resolve """
    if file_upload and (build_name := (file_upload.metadata or {}).get(GENOME_BUILD)):
        # Validation stored the canonical name, which is the PK - don't re-run ambiguous alias resolution
        return GenomeBuild.objects.get(pk=build_name)
    return None


def get_metadata_source(file_upload) -> Optional[str]:
    """ The caller/software the submitter declared, matched against VCFSourceSettings.source_regex """
    if file_upload:
        return (file_upload.metadata or {}).get(SOURCE)
    return None
