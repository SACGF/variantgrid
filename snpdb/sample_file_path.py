"""
Sample files (BAM/CRAM etc, SampleFilePath) beyond the per-sample formset: mapping stored paths to the user's
local view of them (UserDataPrefix) for IGV links, and bulk creation for every sample in a VCF from a
%-style pattern (#647).

Entry points: get_paths_and_user_data_paths, get_example_replacements, validate_sample_file_path_pattern,
resolve_sample_file_paths, create_sample_file_paths
"""
from collections import OrderedDict
from dataclasses import dataclass
from typing import Optional

from django.db.models import QuerySet

from library.common_dir import get_common_prefix_dirs
from snpdb.models.models_enums import SampleFileType
from snpdb.models.models_user_settings import UserDataPrefix
from snpdb.models.models_vcf import VCF, Sample, SampleFilePath

# Keys of Sample._get_sample_formatter_params - the same placeholders as the grid sample label template
PATTERN_KEYS = ("sample_id", "sample", "vcf_sample_name", "patient_id", "patient_code", "patient",
                "specimen_id", "specimen")


def get_paths_and_user_data_paths(user, file_paths):
    replace_dict = UserDataPrefix.get_replace_dict(user)

    paths_and_user_data_paths = OrderedDict()
    for from_path in file_paths:
        to_path = from_path
        for prefix, replacement in replace_dict.items():
            if to_path.startswith(prefix):
                to_path = to_path.replace(prefix, replacement)
                break

        paths_and_user_data_paths[from_path] = to_path

    return paths_and_user_data_paths


def get_example_replacements(user):
    example_replacements = {}
    sample_qs = Sample.filter_for_user(user)
    file_paths = SampleFilePath.objects.filter(sample__in=sample_qs).values_list("file_path", flat=True)
    if file_paths.exists():
        fp_set = set(file_paths)
        prefix_dirs = get_common_prefix_dirs(fp_set)
        example_replacements = get_paths_and_user_data_paths(user, prefix_dirs)

    return example_replacements


@dataclass
class SampleFilePathResolution:
    sample: Sample
    file_path: Optional[str]  # None when the pattern could not be resolved for this sample
    error: Optional[str]  # eg "No patient_code" - shown in the preview, sample skipped on save
    already_exists: bool  # a SampleFilePath with this sample/file_type/file_path is already stored

    @property
    def is_new(self) -> bool:
        return self.file_path is not None and not self.already_exists


class _KeyRecordingDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.used_keys = set()

    def __getitem__(self, key):
        self.used_keys.add(key)
        return super().__getitem__(key)


def validate_sample_file_path_pattern(pattern: str):
    """ Raises ValueError if the pattern can't give each sample its own path """
    trial_params = _KeyRecordingDict({k: 1 if k.endswith("_id") else k for k in PATTERN_KEYS})
    try:
        pattern % trial_params
    except KeyError as e:
        raise ValueError(f"Unknown placeholder {e}, valid placeholders are: {', '.join(PATTERN_KEYS)}") from e
    except (ValueError, TypeError) as e:
        raise ValueError(f"Invalid pattern ({e}) - placeholders are written %(name)s and a literal '%' as %%") from e
    if not trial_params.used_keys:
        raise ValueError("Pattern needs at least one placeholder, eg %(vcf_sample_name)s, "
                         "otherwise every sample would get the same file")


def get_vcf_samples_with_files(vcf: VCF) -> QuerySet[Sample]:
    """ In VCF column order, with what the pattern placeholders and existing files table read """
    return vcf.sample_set.select_related(
        "patient__external_pk", "extraction__specimen",
    ).prefetch_related("samplefilepath_set").order_by("pk")


def _resolve_sample_file_path(sample: Sample, pattern: str) -> tuple[Optional[str], Optional[str]]:
    """ returns (file_path, error) """
    params = {k: v for k, v in sample._get_sample_formatter_params().items() if v not in (None, "")}
    try:
        return pattern % params, None
    except KeyError as e:
        return None, f"No {e.args[0]}"
    except (ValueError, TypeError) as e:
        return None, str(e)


def resolve_sample_file_paths(vcf: VCF, pattern: str, file_type: SampleFileType) -> list[SampleFilePathResolution]:
    resolutions = []
    for sample in get_vcf_samples_with_files(vcf):
        file_path, error = _resolve_sample_file_path(sample, pattern)
        already_exists = any(sfp.file_type == file_type and sfp.file_path == file_path
                             for sfp in sample.samplefilepath_set.all())
        resolutions.append(SampleFilePathResolution(sample=sample, file_path=file_path, error=error,
                                                    already_exists=already_exists))
    return resolutions


def create_sample_file_paths(resolutions: list[SampleFilePathResolution], file_type: SampleFileType,
                             label: Optional[str]) -> list[SampleFilePath]:
    """ Additive only - errored and already present samples are skipped, existing rows are left alone.
        SampleFilePath has no unique constraint (databases already hold duplicates) so already_exists is the dedupe """
    new_sample_file_paths = [
        SampleFilePath(sample=r.sample, file_type=file_type, label=label or None, file_path=r.file_path)
        for r in resolutions if r.is_new
    ]
    return SampleFilePath.objects.bulk_create(new_sample_file_paths)
