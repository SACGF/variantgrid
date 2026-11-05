import os
import re
from collections.abc import Callable
from copy import deepcopy
from dataclasses import dataclass
from enum import Enum
from typing import Any, Optional

from django.conf import settings

from genes.enums import AnnotationConsortium
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
    "spliceai_indel": "spliceai_indel",  # the spliceai pin reads its version off the SNV file alone
    "splicevault": "splicevault",
    "transcript_blocklist": "transcript_blocklist",
    "utr_annotator": "utr_annotator",
}

# Data files with a version in the name that don't already have a VariantAnnotationVersion field.
# Anything that doesn't match falls back to the basename rather than failing the pipeline.
_VERSION_PIN_FIELDS = {
    "dbscsnv": ("dbscsnv", re.compile(r"^dbscSNV(?P<version>[\d.]+?)_", flags=re.IGNORECASE)),
    "gnomad_sv": ("gnomad_sv", re.compile(r"^gnomad[._]v?(?P<version>[\d.]+)[._]sv", flags=re.IGNORECASE)),
    # The second --custom over gnomAD SV - usually the same file as gnomad_sv, but configured separately
    "gnomad_sv_name": ("gnomad_sv_name", re.compile(r"^gnomad[._]v?(?P<version>[\d.]+)[._]sv", flags=re.IGNORECASE)),
    "mastermind": ("mastermind", re.compile(r"^mastermind_cited_variants_reference-(?P<version>[\d.]+)-")),
    "topmed": ("topmed", re.compile(r"^TOPMED_(?:GRCh37|GRCh38)_(?P<version>\d{8})", flags=re.IGNORECASE)),
    "uk10k": ("uk10k", re.compile(r"^UK10K_COHORT\.(?P<version>\d{8})")),
}

# phastCons/phyloP bigwigs are installed and swapped as a set, and which ones apply is a function of
# the build, so they share a single pin rather than a field each
_CONSERVATION_KEYS = ("phastcons100way", "phastcons30way", "phastcons46way",
                      "phylop100way", "phylop30way", "phylop46way")


# vep_config data whose licence rules out some uses, by the settings that declare a deployment makes that use. Where
# the deployment makes one, the data is never passed to VEP (so never pinned or shown)
LICENCE_RESTRICTED_VEP_CONFIG_KEYS = {
    # AlphaGenome merged splicing: non-commercial use only, and not for diagnostic procedures
    "alphagenome_splicing": ("DEPLOYMENT_COMMERCIAL_USE", "DEPLOYMENT_DIAGNOSTIC_USE"),
}


def licence_excluded_keys() -> list[str]:
    """ Licence-restricted data whose licence rules out a use this deployment makes """
    return [key for key, ruled_out_uses in LICENCE_RESTRICTED_VEP_CONFIG_KEYS.items()
            if any(getattr(settings, use) for use in ruled_out_uses)]


def licence_excluded_pin_fields() -> list[str]:
    """ VariantAnnotationVersion pins of the data licence_excluded_keys leaves out """
    return [_BASENAME_PIN_FIELDS[key] for key in licence_excluded_keys()]


def allowed_vep_config_data(vep_config_data: dict) -> dict:
    """ vep_config with the data cleared that the deployment may not use """
    if excluded_keys := licence_excluded_keys():
        return {key: None if key in excluded_keys else value for key, value in vep_config_data.items()}
    return vep_config_data


def configured_licence_excluded_keys(vep_config_data: dict) -> list[str]:
    """ Licence-restricted data configured, but left out because the deployment may not use it """
    return [key for key in licence_excluded_keys() if vep_config_data.get(key)]


def withheld_licence_restricted_keys(variant_annotation_version) -> list[str]:
    """ Licence-restricted data a VariantAnnotationVersion can't offer: the deployment may not use it, or the version
        wasn't annotated with it (its pin is null) - a version annotated before the deployment's use changed still
        holds it """
    excluded_keys = licence_excluded_keys()
    return [key for key in LICENCE_RESTRICTED_VEP_CONFIG_KEYS
            if key in excluded_keys or not getattr(variant_annotation_version, _BASENAME_PIN_FIELDS[key])]


def _pin_basename(path: Optional[str]) -> Optional[str]:
    if not path:
        return None
    return os.path.basename(path)


def vep_component_version_kwargs(build_settings: dict) -> dict:
    """ VariantAnnotationVersion pins for the vep_config data the ##VEP= header doesn't report - files without
        a dedicated version field, and whether SIFT runs. Takes the raw settings.ANNOTATION[build] dict (@see #462);
        the scalar settings are VAV_SETTINGS.

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
    return kwargs


class SettingRole(Enum):
    """ What happens to the VariantAnnotationVersions in use when the setting changes """
    # In the lookup and checked against the running VEP / the annotated file's header: the version halts with
    # VEPVersionMismatchError until a new one is made
    PIN = "pin"
    # In the lookup, but read off the version wherever it's used (the VEP command line or the import) and left out of
    # the checks: the version keeps going and create_new_variant_annotation_version mints a new one
    SNAPSHOT = "snapshot"
    # Stamped on creation only, never looked up: versions made before the field existed (null) still match, and a
    # changed setting takes effect once something else mints a new version
    DEFAULT = "default"


@dataclass(frozen=True)
class VAVSettingContext:
    columns_version: int
    annotation_consortium: str

    @classmethod
    def for_vep_config(cls, vep_config: "VEPConfig") -> "VAVSettingContext":
        return cls(vep_config.columns_version, vep_config.annotation_consortium)


@dataclass(frozen=True)
class VAVSetting:
    field: str  # VariantAnnotationVersion field
    setting: str  # the Django setting it records
    role: SettingRole
    value: Callable[[VAVSettingContext], Any]  # what a version created now records


def _min_pathogenic_rankscore(context: VAVSettingContext) -> Optional[float]:
    """ Rankscores count towards predictions_num_pathogenic at columns_version 2/3 only (1 uses the predictions,
        4+ the raw scores), so elsewhere the threshold isn't recorded and changing it doesn't make a new version """
    if context.columns_version in (2, 3):
        return settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE
    return None


def _gencode_subset(context: VAVSettingContext) -> Optional[str]:
    if context.annotation_consortium == AnnotationConsortium.ENSEMBL:
        return getattr(settings, "ANNOTATION_VEP_ENSEMBL_GENCODE", None)
    return None


def _distance(_context: VAVSettingContext) -> int:
    distance = getattr(settings, "ANNOTATION_VEP_DISTANCE", None)
    return 5000 if distance is None else distance


# Every Django setting a VariantAnnotationVersion records, with what a change to it does to the versions in use.
# The data files are pinned separately above (vep_component_version_kwargs); this is the scalar settings.
VAV_SETTINGS: tuple[VAVSetting, ...] = (
    VAVSetting("pick_order", "ANNOTATION_VEP_PICK_ORDER", SettingRole.PIN,
               lambda c: settings.ANNOTATION_VEP_PICK_ORDER),
    VAVSetting("vep_args", "ANNOTATION_VEP_ARGS", SettingRole.PIN,
               lambda c: " ".join(settings.ANNOTATION_VEP_ARGS) or None),
    VAVSetting("sv_max_size", "ANNOTATION_VEP_SV_MAX_SIZE", SettingRole.PIN,
               lambda c: settings.ANNOTATION_VEP_SV_MAX_SIZE),
    # gnomAD SV overlap_cutoff, applied by VEP to the gnomAD SV and by SVOverlapProcessor to ours (#2589)
    VAVSetting("sv_overlap_min_fraction", "ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION", SettingRole.PIN,
               lambda c: settings.ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION),
    # A symbolic del/dup/inv shorter than this goes to the STANDARD pipeline (annotation/annotation_pipeline_routing.py)
    VAVSetting("structural_variant_min_size", "ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE", SettingRole.PIN,
               lambda c: settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE),
    VAVSetting("distance", "ANNOTATION_VEP_DISTANCE", SettingRole.SNAPSHOT, _distance),  # VEP --distance
    VAVSetting("gencode_subset", "ANNOTATION_VEP_ENSEMBL_GENCODE", SettingRole.SNAPSHOT, _gencode_subset),
    VAVSetting("min_pathogenic_rankscore", "ANNOTATION_MIN_PATHOGENIC_RANKSCORE", SettingRole.SNAPSHOT,
               _min_pathogenic_rankscore),
    VAVSetting("sv_overlap_same_type", "ANNOTATION_VEP_SV_OVERLAP_SAME_TYPE", SettingRole.SNAPSHOT,
               lambda c: settings.ANNOTATION_VEP_SV_OVERLAP_SAME_TYPE),
    # Which overlapping gnomAD SV supplies the gnomAD columns (SVOverlapProcessor)
    VAVSetting("sv_overlap_single_value_method", "ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD",
               SettingRole.SNAPSHOT, lambda c: settings.ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD),
    # An explicit (sequence) variant whose ref or alt is at least this long goes to the STRUCTURAL_VARIANT pipeline;
    # the same setting as structural_variant_min_size, whose pin makes the new version when it changes (#1708)
    VAVSetting("explicit_structural_variant_min_size", "ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE", SettingRole.DEFAULT,
               lambda c: settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE),
)

SETTINGS_SNAPSHOT_FIELDS = frozenset(s.field for s in VAV_SETTINGS if s.role is SettingRole.SNAPSHOT)


def vav_settings_kwargs(context: VAVSettingContext, *roles: SettingRole) -> dict:
    """ field -> what a VariantAnnotationVersion created now records, for the settings with these roles """
    return {s.field: s.value(context) for s in VAV_SETTINGS if s.role in roles}


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
