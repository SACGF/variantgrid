"""
Annotation data files the resolved settings no longer refer to: a superseded dbNSFP or COSMIC release, an older
gnomAD VCF, an `old/` directory from a previous upgrade. They keep tens of GB per build until the disk fills.

Walks `annotation_data/<build>` for each enabled build and `annotation_data/all_builds` under
ANNOTATION_VEP_BASE_DIR (so never the VEP cache or code), and keeps every file that the enabled builds'
ANNOTATION entries or VCF_IMPORT_COMMON_FILTERS name - with its index files and sidecars (`<file>.tbi`,
`<file>.csi`) - and anything under a named directory. Working from the resolved settings means a deployment pinned
to an older columns version keeps that version's files. It only reports: a file may be kept on purpose.

Entry point: `unused_annotation_data`, reported by `snpdb/checks.py:check_unused_annotation_data` under
`manage.py check --deploy`.
"""
import os
from dataclasses import dataclass

from django.conf import settings

ALL_BUILDS_DIR = "all_builds"
# Index files and checksums named after the data file they belong to
SIDECAR_SUFFIXES = (".tbi", ".csi", ".gzi", ".fai", ".md5")


@dataclass(frozen=True)
class UnusedAnnotationData:
    path: str
    size: int  # bytes, the total for a directory


def _setting_paths(value) -> set[str]:
    """ Every string in a nested settings value, as an absolute path - a relative one is under
        ANNOTATION_VEP_BASE_DIR, as VEPConfig resolves it. Strings that aren't paths match no file """
    if isinstance(value, str):
        return {os.path.normpath(os.path.join(settings.ANNOTATION_VEP_BASE_DIR, value))}
    if isinstance(value, dict):
        value = list(value.values())
    if isinstance(value, (list, tuple, set)):
        return set().union(*(_setting_paths(v) for v in value))
    return set()


def _enabled_build_names() -> list[str]:
    return [name for name, config in settings.ANNOTATION.items() if config.get("enabled")]


def used_annotation_paths() -> set[str]:
    used = set()
    for build_name in _enabled_build_names():
        used |= _setting_paths(settings.ANNOTATION[build_name])
        used |= _setting_paths(settings.VCF_IMPORT_COMMON_FILTERS.get(build_name, {}))
    return used


def _sidecar_of(path: str) -> str | None:
    """ The data file an index / checksum belongs to """
    for suffix in SIDECAR_SUFFIXES:
        if path.endswith(suffix):
            return path.removesuffix(suffix)
    return None


def _is_used(path: str, used: set[str]) -> bool:
    """ The file itself, its index / sidecar, or something under a used directory """
    return path in used or _sidecar_of(path) in used or any(path.startswith(u + os.sep) for u in used)


def _contains_used(directory: str, used: set[str]) -> bool:
    return any(u.startswith(directory + os.sep) for u in used)


def _size(path: str) -> int:
    if not os.path.isdir(path):
        return os.path.getsize(path)
    return sum(os.path.getsize(os.path.join(root, f)) for root, _dirs, files in os.walk(path) for f in files)


def _unused_under(directory: str, used: set[str]) -> list[UnusedAnnotationData]:
    """ A directory holding nothing used is reported whole, rather than file by file, and a file's index and
        sidecars are counted in its own entry """
    unused_sizes: dict[str, int] = {}
    nested = []
    for entry in sorted(os.scandir(directory), key=lambda e: e.name):
        path = os.path.normpath(entry.path)
        if _is_used(path, used):
            continue
        if entry.is_dir(follow_symlinks=False) and _contains_used(path, used):
            nested.extend(_unused_under(path, used))
            continue
        # Sorted, so "x.vcf.gz" is seen before "x.vcf.gz.tbi"
        owner = data_file if (data_file := _sidecar_of(path)) in unused_sizes else path
        unused_sizes[owner] = unused_sizes.get(owner, 0) + _size(path)
    return [UnusedAnnotationData(path, size) for path, size in unused_sizes.items()] + nested


def unused_annotation_data() -> dict[str, list[UnusedAnnotationData]]:
    """ {annotation data directory: what in it nothing refers to}, only for directories with something to report """
    used = used_annotation_paths()
    data_dir = os.path.join(settings.ANNOTATION_VEP_BASE_DIR, "annotation_data")
    unused = {}
    for name in [*_enabled_build_names(), ALL_BUILDS_DIR]:
        directory = os.path.join(data_dir, name)
        if os.path.isdir(directory) and (found := _unused_under(directory, used)):
            unused[directory] = found
    return unused
