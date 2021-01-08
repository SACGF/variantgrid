from collections import OrderedDict

from library.common_dir import get_common_prefix_dirs
from snpdb.models import Sample, UserDataPrefix


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

    bam_file_paths = Sample.filter_for_user(user).filter(bam_file_path__isnull=False).values_list("bam_file_path", flat=True)
    if bam_file_paths.exists():
        bfp_set = set(bam_file_paths)
        prefix_dirs = get_common_prefix_dirs(bfp_set)
        example_replacements = get_bam_paths_and_user_data_paths(user, prefix_dirs)

    return example_replacements
