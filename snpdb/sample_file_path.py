from collections import OrderedDict

from library.common_dir import get_common_prefix_dirs
from snpdb.models import Sample, UserDataPrefix, SampleFilePath


def get_bam_paths_and_user_data_paths(user, bam_file_paths):
    replace_dict = UserDataPrefix.get_replace_dict(user)

    bam_paths_and_user_data_paths = OrderedDict()
    for from_bam in bam_file_paths:
        to_bam = from_bam
        for prefix, replacement in replace_dict.items():
            if to_bam.startswith(prefix):
                to_bam = to_bam.replace(prefix, replacement)
                break

        bam_paths_and_user_data_paths[from_bam] = to_bam

    return bam_paths_and_user_data_paths


def get_example_replacements(user):
    example_replacements = {}
    sample_qs = Sample.filter_for_user(user)
    file_paths = SampleFilePath.objects.filter(sample__in=sample_qs).values_list("file_path", flat=True)
    if file_paths.exists():
        fp_set = set(file_paths)
        prefix_dirs = get_common_prefix_dirs(fp_set)
        example_replacements = get_bam_paths_and_user_data_paths(user, prefix_dirs)

    return example_replacements
