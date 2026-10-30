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


def parse_cosmic_version_from_filename(path: str) -> Optional[int]:
    """ COSMIC release from the distributed VCF name, e.g.
        CosmicCodingMuts_v95_20211101_grch37.normal.vcf.gz -> 95
        Cosmic_GenomeScreensMutant_Normal_v101_GRCh37.vcf.gz -> 101 """
    basename = os.path.basename(path)
    if m := re.match(r"^Cosmic.*_v(\d{2,})_.*.vcf.gz", basename):
        return int(m.group(1))
    return None


# vep_config data files whose name carries no version - we pin the basename, so swapping the installed
# file registers as a change, as does removing it (which drops that data's columns via has_data_files)
_BASENAME_PIN_FIELDS = {
    "alphagenome_splicing": "alphagenome_splicing",
    "avi": "avi",
    "eve": "eve",
    "maxentscan": "maxentscan",
    "promoter_ai": "promoter_ai",
    "protvar": "protvar",
    "repeatmasker": "repeat_masker",
    "splicevault": "splicevault",
    "transcript_blocklist": "transcript_blocklist",
    "utr_annotator": "utr_annotator",
}

# Data files with a version in the name that don't already have a VariantAnnotationVersion field.
# Anything that doesn't match falls back to the basename rather than failing the pipeline.
_VERSION_PIN_FIELDS = {
    "dbscsnv": ("dbscsnv", re.compile(r"^dbscSNV(?P<version>[\d.]+?)_", flags=re.IGNORECASE)),
    "gnomad_sv": ("gnomad_sv", re.compile(r"^gnomad[._]v?(?P<version>[\d.]+)[._]sv", flags=re.IGNORECASE)),
    "mastermind": ("mastermind", re.compile(r"^mastermind_cited_variants_reference-(?P<version>[\d.]+)-")),
    "topmed": ("topmed", re.compile(r"^TOPMED_(?:GRCh37|GRCh38)_(?P<version>\d{8})", flags=re.IGNORECASE)),
    "uk10k": ("uk10k", re.compile(r"^UK10K_COHORT\.(?P<version>\d{8})")),
}

# phastCons/phyloP bigwigs are installed and swapped as a set, and which ones apply is a function of
# the build, so they share a single pin rather than a field each
_CONSERVATION_KEYS = ("phastcons100way", "phastcons30way", "phastcons46way",
                      "phylop100way", "phylop30way", "phylop46way")


# vep_config data licensed for non-commercial use only - never passed to VEP (so never pinned or shown) unless
# settings.ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED
NON_COMMERCIAL_VEP_CONFIG_KEYS = ("alphagenome_splicing",)


def allowed_vep_config_data(vep_config_data: dict) -> dict:
    """ vep_config with the non-commercial data cleared when the deployment may not use it """
    if settings.ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED:
        return vep_config_data
    return {key: None if key in NON_COMMERCIAL_VEP_CONFIG_KEYS else value for key, value in vep_config_data.items()}


def disallowed_non_commercial_keys(vep_config_data: dict) -> list[str]:
    """ Non-commercial data configured, but left out because the deployment may not use it """
    if settings.ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED:
        return []
    return [key for key in NON_COMMERCIAL_VEP_CONFIG_KEYS if vep_config_data.get(key)]


def withheld_non_commercial_keys(variant_annotation_version) -> list[str]:
    """ Non-commercial data a VariantAnnotationVersion can't offer: it wasn't annotated with it (its pin is null),
        or the deployment may not use it - a version annotated before the setting was turned off still holds it """
    if not settings.ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED:
        return list(NON_COMMERCIAL_VEP_CONFIG_KEYS)
    return [key for key in NON_COMMERCIAL_VEP_CONFIG_KEYS
            if not getattr(variant_annotation_version, _BASENAME_PIN_FIELDS[key])]


def _pin_basename(path: Optional[str]) -> Optional[str]:
    if not path:
        return None
    return os.path.basename(path)


def vep_component_version_kwargs(build_settings: dict) -> dict:
    """ VariantAnnotationVersion pins for the VEP command line components the ##VEP= header doesn't
        report - data files without a dedicated version field, and the settings that alter the command
        line. Takes the raw settings.ANNOTATION[build] dict (@see #462).

        Values come from the configured paths alone, with no filesystem access, so this returns the same
        answer on a web node as on an annotation node. """

    vep_config_data = allowed_vep_config_data(build_settings.get("vep_config", {}))
    kwargs = {field: _pin_basename(vep_config_data.get(key))
              for key, field in _BASENAME_PIN_FIELDS.items()}

    for key, (field, pattern) in _VERSION_PIN_FIELDS.items():
        basename = _pin_basename(vep_config_data.get(key))
        if basename and (m := pattern.match(basename)):
            kwargs[field] = m.group("version")
        else:
            kwargs[field] = basename

    conservation = [b for k in _CONSERVATION_KEYS if (b := _pin_basename(vep_config_data.get(k)))]
    kwargs["conservation"] = ",".join(conservation) or None

    # Mirrors _get_vep_fasta - the vep_config override only exists on deployments pinned below VEP 112
    kwargs["fasta"] = _pin_basename(vep_config_data.get("fasta") or build_settings.get("reference_fasta"))
    kwargs["sift_enabled"] = bool(vep_config_data.get("sift", False))
    kwargs["pick_order"] = settings.ANNOTATION_VEP_PICK_ORDER
    kwargs["vep_args"] = " ".join(settings.ANNOTATION_VEP_ARGS) or None
    kwargs["sv_max_size"] = settings.ANNOTATION_VEP_SV_MAX_SIZE
    kwargs["sv_overlap_min_fraction"] = settings.ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION
    return kwargs


class VEPConfig:

    def __init__(self, genome_build: GenomeBuild):
        self.annotation_consortium = genome_build.annotation_consortium
        self.genome_build = genome_build
        self.vep_version = int(settings.ANNOTATION_VEP_VERSION)

        # We'll strip out any config - anything left is files/data
        vep_config = deepcopy(allowed_vep_config_data(genome_build.settings["vep_config"]))
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

    @property
    def cosmic_version(self) -> Optional[int]:
        """ COSMIC release of the configured VCF - gates which INFO field carries the sample count
            (see the cosmic_count VEPColumnDefs). None where COSMIC isn't configured (T2T) or the
            release can't be read off the filename, which leaves the count columns ungated. """
        try:
            cosmic_path = self["cosmic"]
        except KeyError:
            return None
        return parse_cosmic_version_from_filename(cosmic_path)
