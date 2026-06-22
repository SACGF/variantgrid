import os
import re
from copy import deepcopy
from typing import Optional

from django.conf import settings

from snpdb.models.models_genome import GenomeBuild


def parse_gnomad_version_from_filename(path: str) -> Optional[str]:
    basename = os.path.basename(path)
    if m := re.match(r"^gnomad(.*?)_(GRCh37|GRCh38|hg19|hg38|T2T-CHM13v2.0)",
                     basename, flags=re.IGNORECASE):
        return m.group(1)
    return None


class VEPConfig:

    def __init__(self, genome_build: GenomeBuild):
        self.annotation_consortium = genome_build.annotation_consortium
        self.genome_build = genome_build
        self.vep_version = int(settings.ANNOTATION_VEP_VERSION)

        # We'll strip out any config - anything left is files/data
        vep_config = deepcopy(genome_build.settings["vep_config"])
        self.use_sift = vep_config.pop("sift", False)
        self.cache_version = vep_config.pop("cache_version", self.vep_version)

        self.vep_data = vep_config
        self.columns_version = genome_build.settings["columns_version"]

    def __getitem__(self, key):
        """ Throws KeyError if missing """
        value = self.vep_data[key]  # All callers need to catch KeyError
        if value is None:
            raise KeyError(key)
        return os.path.join(settings.ANNOTATION_VEP_BASE_DIR, value)

    @property
    def gnomad4_minor_version(self) -> str:
        try:
            gnomad4_path = self["gnomad4"]
        except KeyError:
            return "4.0"
        version = parse_gnomad_version_from_filename(gnomad4_path)
        if version is None:
            return "4.0"
        return version
