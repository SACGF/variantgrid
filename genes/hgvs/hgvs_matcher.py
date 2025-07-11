import logging
import re
import sys
from dataclasses import dataclass
from typing import Optional, Union, Tuple

from django.conf import settings
from django.core.cache import cache
from django.db.models import Max, Min

from genes.hgvs import HGVSVariant, CHGVS, HGVSImplementationException, HGVSNomenclatureException
from genes.hgvs.biocommons_hgvs.hgvs_converter_biocommons import BioCommonsHGVSConverter
from genes.hgvs.hgvs_converter import HGVSConverterType, HgvsMatchRefAllele, HGVSConverter, HgvsOriginallyNormalized
from genes.hgvs.hgvs_converter_combo import ComboCheckerHGVSConverter
from genes.hgvs.pyhgvs.hgvs_converter_pyhgvs import PyHGVSConverter
from genes.models import TranscriptVersion, Transcript, LRGRefSeqGene, BadTranscript, \
    NoTranscript, TranscriptParts
from genes.models_enums import HGVSKind
from genes.transcripts_utils import transcript_is_lrg, looks_like_transcript, looks_like_hgvs_prefix
from library.constants import WEEK_SECS
from library.log_utils import report_exc_info
from library.utils import clean_string, FormerTuple
from snpdb.clingen_allele import get_clingen_allele_from_hgvs, \
    ClinGenAlleleServerException, ClinGenAlleleAPIException, get_clingen_allele_for_variant_coordinate, \
    clingen_check_variant_length
from snpdb.models import Variant, ClinGenAllele
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import VariantCoordinate


@dataclass
class FakeTranscriptVersion:
    transcript_id: str
    version: int
    hgvs_ok: bool = False
    gene_symbol = None

    @property
    def accession(self) -> str:
        return f"{self.transcript_id}.{self.version}"


@dataclass(frozen=True)
class VariantCoordinateAndDetails(FormerTuple):
    variant_coordinate: VariantCoordinate
    transcript_accession: str
    kind: str
    used_converter_type: Optional[HGVSConverterType]
    method: str  # human readable contains a trail of what was tried
    matches_reference: Union[bool, HgvsMatchRefAllele]
    originally_normalized: Union[bool, HgvsOriginallyNormalized]

    @property
    def as_tuple(self) -> tuple:
        return self.variant_coordinate, self.transcript_accession, self.kind, self.used_converter_type, self.method, self.matches_reference, self.originally_normalized


class ClinGenHGVSConverter(BioCommonsHGVSConverter):
    @classmethod
    def build_supported(cls, genome_build) -> bool:
        return genome_build in cls.supported_builds

    @classmethod
    @property
    def supported_builds(cls) -> list[GenomeBuild]:
        return [GenomeBuild.grch37(), GenomeBuild.grch38()]

    def __init__(self, genome_build, local_resolution=False, clingen_resolution=True):
        if not self.build_supported(genome_build):
            supported_builds = ", ".join([str(gb) for gb in self.supported_builds])
            raise ValueError(f"ClinGen: unsupported {genome_build=}, {supported_builds=}")

        super().__init__(genome_build,
                         local_resolution=local_resolution, clingen_resolution=clingen_resolution)

    def get_hgvs_converter_type(self) -> HGVSConverterType:
        return HGVSConverterType.CLINGEN_ALLELE_REGISTRY

    def get_version(self) -> str:
        return f"From (cached) API, uses {super().get_hgvs_converter_type().name} models)"


class HGVSConverterFactory:
    """
    For choice of PyHGVS vs BioCommons @see https://github.com/SACGF/variantgrid/issues/839
    """

    @staticmethod
    def factory(genome_build: GenomeBuild, hgvs_converter_type: Optional[HGVSConverterType] = None,
                local_resolution=True, clingen_resolution=True) -> HGVSConverter:
        if hgvs_converter_type is None:
            # if settings.DEBUG:  # TODO: Disabled
            #     hgvs_converter_type = HGVSConverterType.COMBO
            # else:
            hgvs_converter_type = settings.HGVS_DEFAULT_METHOD

        if clingen_resolution and not ClinGenHGVSConverter.build_supported(genome_build):
            logging.warning("ClinGen does not support build %s", genome_build)
            clingen_resolution = False

        if isinstance(hgvs_converter_type, str):
            hgvs_converter_type = HGVSConverterType[hgvs_converter_type.upper()]

        if hgvs_converter_type == HGVSConverterType.BIOCOMMONS_HGVS:
            return BioCommonsHGVSConverter(genome_build, local_resolution=local_resolution, clingen_resolution=clingen_resolution)
        elif hgvs_converter_type == HGVSConverterType.PYHGVS:
            return PyHGVSConverter(genome_build, local_resolution=local_resolution, clingen_resolution=clingen_resolution)
        elif hgvs_converter_type == HGVSConverterType.COMBO:
            converters = [
                BioCommonsHGVSConverter(genome_build, local_resolution=local_resolution, clingen_resolution=clingen_resolution),
                PyHGVSConverter(genome_build, local_resolution=local_resolution, clingen_resolution=clingen_resolution),
            ]
            return ComboCheckerHGVSConverter(genome_build, converters, die_on_error=False)
        elif hgvs_converter_type == HGVSConverterType.CLINGEN_ALLELE_REGISTRY:
            return ClinGenHGVSConverter(genome_build)

        raise HGVSImplementationException(f"HGVSConverter type {hgvs_converter_type} not supported")


class VariantResolvingError(ValueError):

    def __init__(self, message: str, technical_message: Optional[str] = None):
        self.technical_message = technical_message
        super().__init__(message)


class HGVSMatcher:
    # "NR", "NM", "NC", "ENST", "LRG_", "XR"}

    TRANSCRIPT_PREFIX = re.compile(r"^(NR|NM|NC|ENST|LRG|XR)", re.IGNORECASE)
    TRANSCRIPT_NO_UNDERSCORE = re.compile(r"^(NM|NC)(\d+)")
    TRANSCRIPT_UNDERSCORE_REPLACE = r"\g<1>_\g<2>"
    # noinspection RegExpSingleCharAlternation
    HGVS_CLEAN_PATTERN = re.compile(
        r"^(?P<transcript_prefix>:NR|NM|NC|ENST|LRG|XR)(?P<transcript_value>[0-9_]*)(\.(?P<transcript_version>[0-9]+))?([(]?(?P<gene_symbol>[^):\.]{2,8})[)]?)?:?(?P<letter>c|g|n|p)[.]?(?P<nomen>[^:\.]*)$",
        re.IGNORECASE
    )
    # Only try to fix things with the sloppy pattern if the clean pattern doesn't match
    # Replace is now done with lambda instead of regex so we can lowercase c,g,n,p
    HGVS_TRANSCRIPT_NO_CDOT = re.compile(r"^(NM_|ENST)\d+.*:\d+")
    HGVS_CONTIG_NO_GDOT = re.compile(r"^NC_\d+.*:\d+")

    C_DOT_REF_ALT_NUC = re.compile("(?P<ref>[gatc]+)>(?P<alt>[gatc=]+)$", re.IGNORECASE)

    C_DOT_REF_DEL_INS_DUP_NUC = re.compile(r"(del(?P<del>[gatc]+))?(?P<op>ins|dup|del)(?P<ins>[gatc]*)$")
    # captures things in teh form of delG, insG, dupG & delGinsC
    # the del pattern at the start is only for delins, as otherwise the op del is captured

    P_HGVS_REMOVAL = re.compile(r"^(?P<main>.*?) (p.*|[(]p.*[)])$")

    HGVS_METHOD_INTERNAL_LIBRARY = "Internally converted using library"
    # External calls to ClinGen
    HGVS_METHOD_CLINGEN_ALLELE_REGISTRY = "ClinGen Allele Registry"

    class TranscriptContigMismatchError(ValueError):
        pass

    def __init__(self, genome_build: GenomeBuild, hgvs_converter_type=None,
                 local_resolution=True, clingen_resolution=True, allow_alternative_transcript_version=True):
        self.genome_build = genome_build
        self.allow_alternative_transcript_version = allow_alternative_transcript_version
        self.attempt_clingen = True  # Set to False if server issues
        self.hgvs_converter = HGVSConverterFactory.factory(genome_build, hgvs_converter_type,
                                                           local_resolution=local_resolution,
                                                           clingen_resolution=clingen_resolution)

    def _clingen_get_variant_coordinate_matches_reference_and_normalized(self, hgvs_string: str, match_ref_allele=None) -> tuple[VariantCoordinate, bool, Union[HgvsOriginallyNormalized, bool]]:
        hgvs_name = self.create_hgvs_variant(hgvs_string)
        if hgvs_name.kind == 'g':
            clingen_check_variant_length(hgvs_string, hgvs_name.length, is_dup=hgvs_name.mutation_type == 'dup')
        cleaned_hgvs = self.hgvs_converter.c_hgvs_remove_gene_symbol(hgvs_string)

        try:
            ca = get_clingen_allele_from_hgvs(cleaned_hgvs, require_allele_id=False)
            variant_coord = ca.get_variant_coordinate(self.genome_build)
            # Was converted to internal, need to return raw strings so standard base validation is OK
            if variant_coord.alt == Variant.REFERENCE_ALT:
                variant_coord = VariantCoordinate(chrom=variant_coord.chrom,
                                                  start=variant_coord.start,
                                                  ref=variant_coord.ref,
                                                  alt=variant_coord.ref)  # ref == alt
            if match_ref_allele is None:
                match_ref_allele = True

            transcript_accession = self.get_transcript_accession(cleaned_hgvs)
            normalized_hgvs = ca.get_c_hgvs_variant(self.hgvs_converter, transcript_accession)
            no_gene_generated = self.hgvs_converter.c_hgvs_remove_gene_symbol(str(normalized_hgvs))
            originally_normalized = HgvsOriginallyNormalized(original_hgvs=self.create_hgvs_variant(cleaned_hgvs),
                                                             normalized_hgvs=self.create_hgvs_variant(no_gene_generated))
            return variant_coord, match_ref_allele, originally_normalized
        except ClinGenAlleleAPIException:
            self.attempt_clingen = False
            raise
        except ClinGenAlleleServerException as cga_se:
            # If it's unknown reference we can just retry with another version, other errors are fatal
            if cga_se.is_unknown_reference():
                transcript_accession = self.hgvs_converter.get_transcript_accession(hgvs_string)
                self._set_clingen_allele_registry_missing_transcript(transcript_accession)
            else:
                if settings.CLINGEN_ALLELE_REGISTRY_REATTEMPT_WITH_ACTUAL_REF and match_ref_allele is None:
                    # Don't do if already swapped (stop infinite recursion)
                    rjson = cga_se.response_json
                    if rjson["errorType"] == 'IncorrectReferenceAllele':
                        actual_allele = rjson['actualAllele']
                        given_allele = rjson['givenAllele']
                        transcript_reference_sequence = rjson["referenceSequence"]
                        hgvs_variant = self.create_hgvs_variant(cleaned_hgvs)
                        if hgvs_variant.ref_allele == given_allele:
                            hgvs_variant.ref_allele = actual_allele
                            hgvs_swapped_ref = hgvs_variant.format()
                            match_ref_allele = HgvsMatchRefAllele(provided_ref=given_allele,
                                                                  calculated_ref=actual_allele,
                                                                  ref_type=f"transcript {transcript_reference_sequence}",
                                                                  ref_source="ClinGen Allele Registry")
                            return self._clingen_get_variant_coordinate_matches_reference_and_normalized(hgvs_swapped_ref, match_ref_allele=match_ref_allele)

                self.attempt_clingen = False
            raise

    @staticmethod
    def _get_renamed_lrg_transcript_hgvs_variant_and_transcript_version(genome_build: GenomeBuild, hgvs_variant: HGVSVariant):
        lrg_identifier = hgvs_variant.transcript

        if transcript_version := LRGRefSeqGene.get_transcript_version(genome_build, lrg_identifier):
            if transcript_version.hgvs_ok:
                # Replace LRG transcript with local RefSeq
                hgvs_variant.transcript = transcript_version.accession
                return hgvs_variant, transcript_version
        return None, None

    def _lrg_get_variant_coordinate_used_transcript_method_and_matches_reference(self, hgvs_variant: HGVSVariant) -> tuple[VariantCoordinate, str, HGVSConverterType, str, Union[bool, HgvsMatchRefAllele], Union[bool, HgvsOriginallyNormalized]]:
        lrg_transcript_accession = hgvs_variant.transcript
        new_hgvs_variant, transcript_version = self._get_renamed_lrg_transcript_hgvs_variant_and_transcript_version(self.genome_build, hgvs_variant)
        if new_hgvs_variant:
            new_hgvs_string = new_hgvs_variant.format()
            method = f"{self.hgvs_converter.description(describe_fallback=False)} as '{new_hgvs_string}' (from LRG_RefSeqGene)"
            variant_coordinate, matches_reference, originally_normalized = self.hgvs_converter.hgvs_to_variant_coordinate_reference_match_and_normalized(new_hgvs_string, transcript_version)
            internal_converter_type = self.hgvs_converter.get_hgvs_converter_type()
            return variant_coordinate, lrg_transcript_accession, internal_converter_type, method, matches_reference, originally_normalized

        hgvs_string = hgvs_variant.format()
        try:
            method = HGVSConverterType.CLINGEN_ALLELE_REGISTRY.name
            variant_coordinate, matches_reference, originally_normalized = self._clingen_get_variant_coordinate_matches_reference_and_normalized(hgvs_string)
            return variant_coordinate, lrg_transcript_accession, HGVSConverterType.CLINGEN_ALLELE_REGISTRY, method, matches_reference, originally_normalized
        except ClinGenAllele.ClinGenAlleleRegistryException as cga_re:
            raise HGVSImplementationException(f"Could not retrieve {hgvs_string} from ClinGen Allele Registry") from cga_re

    @staticmethod
    def _get_clingen_allele_registry_key(transcript_accession: str) -> str:
        return f"clingen_allele_registry_unknown_reference_{transcript_accession}"

    def _clingen_allele_registry_ok(self, transcript_accession: str) -> bool:
        """ So we don't keep hammering their server for the same transcript they don't have, we store in Redis
            cache that they don't have that transcript - this expires over time so we'll check again in case
            they updated their references and now have it """

        if not self.attempt_clingen:
            return False  # Had non-recoverable errors before

        key = HGVSMatcher._get_clingen_allele_registry_key(transcript_accession)
        return not cache.get(key)

    @staticmethod
    def _set_clingen_allele_registry_missing_transcript(transcript_accession: str):
        key = HGVSMatcher._get_clingen_allele_registry_key(transcript_accession)
        cache.set(key, True, timeout=WEEK_SECS)

    def get_variant_coordinate(self, hgvs_string: str) -> VariantCoordinate:
        vcd = self.get_variant_coordinate_and_details(hgvs_string)
        return vcd.variant_coordinate

    @staticmethod
    def _get_sort_key_transcript_version_and_methods(version, prefer_local=True, closest=False):
        def get_sort_key(item):
            tv, method = item

            if version:
                # Ask for 3, have [1, 2, 3, 4, 5, 6]
                # Closest:           4, 5, 3, 6, 2, 1
                # Up then down:      4, 5, 6, 3, 2, 1
                version_distance = abs(version-tv.version)
                prefer_later = tv.version < version
                if closest:
                    sort_keys = [version_distance, prefer_later]
                else:
                    sort_keys = [prefer_later, version_distance]
            else:
                # Latest to earliest
                sort_keys = [-tv.version]

            if method == HGVSMatcher.HGVS_METHOD_INTERNAL_LIBRARY:
                method_sort = 1
            else:
                method_sort = 2

            if prefer_local:
                sort_keys.insert(0, method_sort)
            else:
                sort_keys.append(method_sort)

            return tuple(sort_keys)

        return get_sort_key

    def create_hgvs_variant(self, hgvs_string) -> HGVSVariant:
        return self.hgvs_converter.create_hgvs_variant(hgvs_string)

    def _normalized_check(self, hgvs_variant: HGVSVariant) -> HgvsOriginallyNormalized:
        normalized_hgvs = self.hgvs_converter.normalize(hgvs_variant)
        return HgvsOriginallyNormalized(original_hgvs=hgvs_variant, normalized_hgvs=normalized_hgvs)

    def filter_best_transcripts_and_converter_type_by_accession(self, transcript_accession, prefer_local=True, closest=False) -> list[tuple[TranscriptVersion, HGVSConverterType]]:
        """ Get the best transcripts you'd want to match a HGVS against - assuming you will try multiple in order """

        transcript_id: str
        version: int
        try:
            transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
        except ValueError:
            raise HGVSNomenclatureException(f"Error parsing transcript version from \"{transcript_accession}\"")

        tv_qs = TranscriptVersion.objects.filter(genome_build=self.genome_build, transcript_id=transcript_id)
        tv_by_version = {tv.version: tv for tv in tv_qs}
        if not tv_by_version:
            # If we don't have any in DB - we should check that it's actually real
            try:
                # The only thing we care about is BadTranscript - otherwise can carry on
                TranscriptVersion.raise_bad_or_missing_transcript(transcript_accession)
            except BadTranscript:
                raise  # RefSeq/Ensembl def don't have this transcript
            except NoTranscript:
                pass  # ok

        # When looking at the range of versions to check, we'll use the lowest/highest we've seen in any build
        data = Transcript.objects.filter(pk=transcript_id).aggregate(min_tv=Min("transcriptversion__version"),
                                                                     max_tv=Max("transcriptversion__version"),
                                                                     min_tvsi=Min("transcriptversionsequenceinfo__version"),
                                                                     max_tvsi=Max("transcriptversionsequenceinfo__version"))
        # If we have no local transcript versions we'll just try 1
        version_if_no_local = version or 1
        min_versions = [v for v in [version_if_no_local, data.get("min_tv"), data.get("min_tvsi")] if v is not None]
        max_versions = [v for v in [version_if_no_local, data.get("max_tv"), data.get("max_tvsi")] if v is not None]

        if self.allow_alternative_transcript_version:
            min_version = min(min_versions)
            max_version = max(max_versions)
        else:
            if not version:
                raise HGVSNomenclatureException(f"{transcript_accession=} must have version when allow_alternative_transcript_version=False")
            min_version = version
            max_version = version

        local_converter_type = self.hgvs_converter.get_hgvs_converter_type()
        tv_and_converter_type = []
        for v in range(min_version, max_version+1):
            transcript_version = tv_by_version.get(v)
            if not transcript_version:
                transcript_version = FakeTranscriptVersion(transcript_id=transcript_id, version=v)
            if transcript_version.hgvs_ok and self.hgvs_converter.local_resolution:
                tv_and_converter_type.append((transcript_version, local_converter_type))
            if self.hgvs_converter.clingen_resolution:
                tv_and_converter_type.append((transcript_version, HGVSConverterType.CLINGEN_ALLELE_REGISTRY))

        # TODO: Maybe we should filter transcript versions that have the same length
        sort_key = self._get_sort_key_transcript_version_and_methods(version, prefer_local=prefer_local, closest=closest)
        return sorted(tv_and_converter_type, key=sort_key)

    def get_variant_coordinate_and_details(self, hgvs_string: str) -> VariantCoordinateAndDetails:
        """ Returns variant_coordinate and method for HGVS resolution = """

        transcript_accession = self.hgvs_converter.get_transcript_accession(hgvs_string)
        used_transcript_accession = None
        used_converter_type = None
        method = None
        matches_reference = None
        originally_normalized = None
        hgvs_variant = self.create_hgvs_variant(hgvs_string)
        kind = hgvs_variant.kind
        error_messages: list[str] = []
        combined_error_message = None

        if transcript_is_lrg(transcript_accession):
            variant_coordinate, used_transcript_accession, hgvs_converter_type, method, matches_reference, originally_normalized = self._lrg_get_variant_coordinate_used_transcript_method_and_matches_reference(hgvs_variant)
            used_converter_type = hgvs_converter_type
        elif hgvs_variant.kind in ('c', 'n'):
            if not transcript_accession:
                msg = f"Could not parse \"{hgvs_string}\" c.HGVS requires a transcript or LRG."
                if hgvs_variant.gene:
                    # We should only be here if gene based HGVS searching is disabled
                    msg += f"\nGene appears to be \"{hgvs_variant.gene}\""
                raise ValueError(msg)

            variant_coordinate = None
            hgvs_methods = []
            for tv, potential_converter_type in self.filter_best_transcripts_and_converter_type_by_accession(transcript_accession):
                used_transcript_accession = tv.accession
                hgvs_variant.transcript = tv.accession
                # Ensure that ref is always present so we can give warning about provided reference
                hgvs_string_for_version = hgvs_variant.format(max_ref_length=sys.maxsize)
                if potential_converter_type.is_internal_type():
                    method = self.hgvs_converter.description(describe_fallback=False)
                    variant_coordinate, matches_reference, originally_normalized = self.hgvs_converter.hgvs_to_variant_coordinate_reference_match_and_normalized(hgvs_string_for_version, tv)
                elif potential_converter_type == HGVSConverterType.CLINGEN_ALLELE_REGISTRY:
                    method = HGVSConverterType.CLINGEN_ALLELE_REGISTRY.name
                    if self._clingen_allele_registry_ok(tv.accession):
                        error_message = f"Could not convert \"{hgvs_string}\" using ClinGenAllele Registry"
                        try:
                            variant_coordinate, matches_reference, originally_normalized = self._clingen_get_variant_coordinate_matches_reference_and_normalized(hgvs_string_for_version)
                        except ClinGenAlleleServerException as cga_se:
                            # If it's unknown reference we can just retry with another version, other errors are fatal
                            if cga_se.is_unknown_reference():
                                self._set_clingen_allele_registry_missing_transcript(tv.accession)
                            else:
                                msg = f"{error_message} : {cga_se}"
                                logging.error(msg)
                                error_messages.append(msg)
                        except ClinGenAllele.ClinGenAlleleRegistryException as cgare:
                            # API or other recoverable error - try again w/another transcript
                            msg = f"{error_message} : {cgare}"
                            logging.error(msg)
                            error_messages.append(msg)

                if method:
                    if hgvs_string != hgvs_string_for_version:
                        method += f" as \"{hgvs_string_for_version}\""
                    hgvs_methods.append(method)

                if variant_coordinate:
                    used_converter_type = potential_converter_type
                    break

            if variant_coordinate is None:
                if error_messages:
                    combined_error_message = "\n".join(sorted(set(error_messages)))

                if hgvs_methods:
                    attempts = "\n".join(hgvs_methods)
                    raise VariantResolvingError(f"Could not convert \"{hgvs_string}\" - tried:\n{attempts}", technical_message=combined_error_message)
                else:
                    raise ValueError(f"\"{transcript_accession}\": No transcripts found")
        else:
            # g. HGVS
            method = self.hgvs_converter.description(describe_fallback=False)
            variant_coordinate, matches_reference, originally_normalized = self.hgvs_converter.hgvs_to_variant_coordinate_reference_match_and_normalized(hgvs_string, None)
            used_converter_type = self.hgvs_converter.get_hgvs_converter_type()

        return VariantCoordinateAndDetails(
            variant_coordinate=variant_coordinate.as_internal_symbolic(self.genome_build),
            transcript_accession=used_transcript_accession,
            kind=kind,
            used_converter_type=used_converter_type,
            method=method,
            matches_reference=matches_reference,
            originally_normalized=originally_normalized
        )

    def get_transcript_accession(self, hgvs_string: str) -> str:
        return self.hgvs_converter.get_transcript_accession(hgvs_string)

    def get_transcript_parts(self, hgvs_string: str) -> TranscriptParts:
        accession = self.hgvs_converter.get_transcript_accession(hgvs_string)
        return TranscriptVersion.get_transcript_id_and_version(accession)

    def _lrg_variant_coordinate_to_hgvs(self, variant_coordinate: VariantCoordinate, lrg_identifier: str = None) -> tuple[HGVSVariant, HGVSConverterType, str]:
        if transcript_version := LRGRefSeqGene.get_transcript_version(self.genome_build, lrg_identifier):
            if transcript_version.hgvs_ok:
                hgvs_variant, hgvs_converter_type, hgvs_method = self.variant_coordinate_to_hgvs_used_converter_type_and_method(variant_coordinate, transcript_version.accession)
                if hgvs_variant.transcript != transcript_version.accession:
                    msg = f"Error creating HGVS for {variant_coordinate}, LRG '{lrg_identifier}' asked for HGVS " \
                          f"'{transcript_version.accession}' but got '{hgvs_variant.transcript}'"
                    raise ValueError(msg)
                # Replace with our LRG
                hgvs_variant.transcript = lrg_identifier
                hgvs_variant.gene = str(transcript_version.gene_symbol)
                return hgvs_variant, hgvs_converter_type, hgvs_method

        problems = ["No transcript via LRGRefSeqGene"]

        # Use ClinGen - will raise exception if it can't get it
        if ca := get_clingen_allele_for_variant_coordinate(self.genome_build, variant_coordinate, self, require_allele_id=False):
            if hgvs_variant := ca.get_c_hgvs_variant(self.hgvs_converter, lrg_identifier):
                return hgvs_variant, HGVSConverterType.CLINGEN_ALLELE_REGISTRY, self.HGVS_METHOD_CLINGEN_ALLELE_REGISTRY
            else:
                problems.append(f"{ca} didn't contain HGVS for '{lrg_identifier}'")

        problem_str = ", ".join(problems)
        raise HGVSImplementationException(f"Could not convert {variant_coordinate} to HGVS using '{lrg_identifier}': {problem_str}")

    def variant_coordinate_to_hgvs_used_converter_type_and_method(self, variant_coordinate: VariantCoordinate,
                                                                  transcript_accession: str = None) -> tuple[HGVSVariant, HGVSConverterType, str]:
        """
            returns (hgvs, method) - hgvs is c.HGVS is transcript provided, g.HGVS if not

            This handles using different methods (eg local/ClinGen allele registry)
            We always generate the HGVS with full-length reference bases etc, as we adjust that in HGVSExtra.format()
        """

        hgvs_variant = None
        hgvs_converter_type = None
        hgvs_method = None
        if transcript_accession:
            if transcript_is_lrg(transcript_accession):
                return self._lrg_variant_coordinate_to_hgvs(variant_coordinate, transcript_accession)

            hgvs_methods = {}
            for transcript_version, potential_converter_type in self.filter_best_transcripts_and_converter_type_by_accession(transcript_accession):
                hgvs_method = f"{potential_converter_type}: {transcript_version}"
                hgvs_methods[hgvs_method] = None
                if potential_converter_type.is_internal_type():
                    # Sanity Check - make sure contig is the same
                    vc_contig = self.genome_build.chrom_contig_mappings[variant_coordinate.chrom]
                    if vc_contig != transcript_version.contig:
                        contig_msg = f"Variant contig={vc_contig} while Transcript " \
                                     f"{transcript_version.accession} contig={transcript_version.contig}"
                        raise self.TranscriptContigMismatchError(contig_msg)

                    hgvs_variant = self.hgvs_converter.variant_coordinate_to_c_hgvs(variant_coordinate, transcript_version)
                elif potential_converter_type == HGVSConverterType.CLINGEN_ALLELE_REGISTRY:
                    if self._clingen_allele_registry_ok(transcript_version.accession):
                        error_message = f"Could not convert '{variant_coordinate}' ({transcript_version}) using {potential_converter_type}: %s"
                        # TODO: We could also use VEP then add reference bases on our HGVSs
                        try:
                            if ca := get_clingen_allele_for_variant_coordinate(self.genome_build, variant_coordinate, self, require_allele_id=False):
                                if hgvs_variant := ca.get_c_hgvs_variant(self.hgvs_converter, transcript_version.accession):
                                    # Use our latest symbol as ClinGen can be out of date, and this keeps it consistent
                                    # regardless of whether we use PyHGVS or ClinGen to resolve
                                    if gene_symbol := transcript_version.gene_symbol:
                                        hgvs_variant.gene = str(gene_symbol)
                                    hgvs_method = potential_converter_type
                        except ClinGenAlleleServerException as cga_se:
                            # If it's unknown reference we can just retry with another version, other errors are fatal
                            if cga_se.is_unknown_reference():
                                self._set_clingen_allele_registry_missing_transcript(transcript_version.accession)
                            else:
                                logging.error(error_message, cga_se)
                                hgvs_methods[hgvs_method] = str(cga_se)
                        except ClinGenAllele.ClinGenAlleleRegistryException as cgare:
                            # API or other recoverable error - try again w/another transcript
                            logging.error(error_message, cgare)
                            hgvs_methods[hgvs_method] = str(cgare)

                if hgvs_variant:
                    hgvs_converter_type = potential_converter_type
                    break

            if hgvs_methods:
                if hgvs_variant is None:
                    method_and_errors = []
                    for method, errors in hgvs_methods.items():
                        if errors:
                            method += ": " + errors
                        method_and_errors.append(method)
                    attempts = ", ".join(method_and_errors)
                    raise HGVSImplementationException(f"Could not convert {variant_coordinate} to HGVS - tried: {attempts}")
            else:
                # No methods tried, mustn't have had any transcripts
                TranscriptVersion.raise_bad_or_missing_transcript(transcript_accession)
        else:
            # No transcript = convert to Genomic HGVS
            # We need to return a HGVSVariant - self.variant_coordinate_to_g_hgvs may go down fast path returns a str
            hgvs_variant = self.hgvs_converter.variant_coordinate_to_g_hgvs(variant_coordinate)
            hgvs_converter_type = self.hgvs_converter.get_hgvs_converter_type()
            hgvs_method = self.hgvs_converter.description(describe_fallback=False)

        return hgvs_variant, hgvs_converter_type, hgvs_method

    def variant_to_hgvs_variant_used_converter_type_and_method(self, variant: Variant, transcript_name=None) -> Tuple[HGVSVariant, HGVSConverterType, str]:
        return self.variant_coordinate_to_hgvs_used_converter_type_and_method(variant.coordinate, transcript_name)

    def variant_to_hgvs_variant(self, variant: Variant, transcript_name=None) -> HGVSVariant:
        """ returns c.HGVS is transcript provided, g.HGVS if no transcript"""
        return self.variant_coordinate_to_hgvs_variant(variant.coordinate, transcript_name=transcript_name)

    def variant_coordinate_to_hgvs_variant(self, variant_coordinate: VariantCoordinate, transcript_name=None) -> HGVSVariant:
        variant_coordinate = variant_coordinate.as_external_explicit(self.genome_build)
        return self.variant_coordinate_to_hgvs_used_converter_type_and_method(variant_coordinate, transcript_name)[0]

    def _fast_variant_coordinate_to_g_hgvs(self, refseq_accession, offset, ref, alt) -> str:
        """ This only works for SNPs (ie not indels etc) """
        if ref == alt:
            hgvs_allele = f"{ref}="
        else:
            hgvs_allele = f"{ref}>{alt}"

        if refseq_accession == self.genome_build.mitochondria_accession:
            kind = 'm'
        else:
            kind = 'g'
        return f"{refseq_accession}:{kind}.{offset}{hgvs_allele}"

    def variant_to_g_hgvs(self, variant: Variant) -> str:
        return self.variant_coordinate_to_g_hgvs(variant.coordinate)

    def variant_coordinate_to_g_hgvs(self, variant_coordinate: VariantCoordinate) -> str:
        variant_coordinate = variant_coordinate.as_external_explicit(self.genome_build)
        (chrom, position, ref, alt, _svlen) = variant_coordinate
        if len(alt) == 1 and len(ref) == 1:
            contig = self.genome_build.chrom_contig_mappings[chrom]
            hgvs_str = self._fast_variant_coordinate_to_g_hgvs(contig.refseq_accession, position, ref, alt)
        else:
            hgvs_variant = self.hgvs_converter.variant_coordinate_to_g_hgvs(variant_coordinate)
            hgvs_str = str(hgvs_variant)
        return hgvs_str

    def variant_to_c_hgvs_parts(self, variant: Variant, transcript: Optional[str], throw_on_issue: bool = False) -> Optional[CHGVS]:
        try:
            hgvs_variant = self.variant_to_hgvs_variant(variant, transcript)
            if hgvs_variant:
                return CHGVS(hgvs_variant.format(), transcript)
        except:
            if throw_on_issue:
                raise
            report_exc_info()
        return None

    def clean_hgvs(self, hgvs_string: str) -> tuple[str, list[str]]:
        search_messages = []
        if p_hgvs_remove_match := HGVSMatcher.P_HGVS_REMOVAL.match(hgvs_string):
            hgvs_string = p_hgvs_remove_match.group("main")

        cleaned_hgvs = clean_string(hgvs_string)  # remove non-printable characters
        cleaned_hgvs = cleaned_hgvs.replace(" ", "")  # No whitespace in HGVS
        cleaned_hgvs = cleaned_hgvs.replace("::", ":")  # Fix double colon
        cleaned_hgvs = cleaned_hgvs.replace("__", "_")  # fix double underscore

        # Fix double kind (eg c.c.)
        for kind in HGVSKind:
            cleaned_hgvs = cleaned_hgvs.replace(f"{kind}.{kind}.", f"{kind}.")

        if cleaned_hgvs[0:2].upper() in ("M_", "C_", "R_"):
            cleaned_hgvs = "N" + cleaned_hgvs
        # Lowercase mutation types, e.g. NM_032638:c.1126_1133DUP - won't matter if also changes gene name as that's
        # case-insensitive
        MUTATION_TYPES = ["ins", "del", "dup", "inv"]  # Will also handle delins and del...ins
        for mt in MUTATION_TYPES:
            cleaned_hgvs = cleaned_hgvs.replace(mt.upper(), mt)

        # Handle unbalanced brackets or >1 of each type
        open_bracket = cleaned_hgvs.count("(")
        close_bracket = cleaned_hgvs.count(")")
        if open_bracket - close_bracket or open_bracket > 1 or close_bracket > 1:
            # Best bet is to just strip all of them
            cleaned_hgvs = cleaned_hgvs.replace("(", "").replace(")", "")

        # if transcript_prefix_match := self.TRANSCRIPT_PREFIX.search(cleaned_hgvs):
        #     transcript_prefix = transcript_prefix_match.group(0)
        #     if transcript_prefix != transcript_prefix.upper():
        #         cleaned_hgvs = self.TRANSCRIPT_PREFIX.sub(transcript_prefix.upper(), cleaned_hgvs)

        cleaned_hgvs = self.TRANSCRIPT_NO_UNDERSCORE.sub(self.TRANSCRIPT_UNDERSCORE_REPLACE, cleaned_hgvs)
        # r"\g<1>:\g<2>.\g<3>"

        if match := self.HGVS_CLEAN_PATTERN.match(cleaned_hgvs):
            cleaned_text = match.group("transcript_prefix").upper() + match.group("transcript_value")
            if transcript_version := match.group("transcript_version"):
                cleaned_text += "." + transcript_version
            if gene_symbol := match.group("gene_symbol"):
                cleaned_text += "(" + gene_symbol + ")"
            cleaned_text += ":" + match.group("letter").lower() + "." + match.group("nomen")
            cleaned_hgvs = cleaned_text

        def fix_ref_alt(m):
            return m.group('ref').upper() + '>' + m.group('alt').upper()

        def fix_del_ins(m):
            parts = []
            if del_nucs := m.group('del'):
                parts.append(f"del{del_nucs.upper()}")
            parts.append(m.group('op'))
            parts.append(m.group('ins').upper())
            return "".join(parts)

        cleaned_hgvs = self.C_DOT_REF_ALT_NUC.sub(fix_ref_alt, cleaned_hgvs)
        cleaned_hgvs = self.C_DOT_REF_DEL_INS_DUP_NUC.sub(fix_del_ins, cleaned_hgvs)

        # If it contains a transcript and a colon, but no "c." then add it
        if self.HGVS_TRANSCRIPT_NO_CDOT.match(cleaned_hgvs):
            cleaned_hgvs = cleaned_hgvs.replace(":", ":c.")
        elif self.HGVS_CONTIG_NO_GDOT.match(cleaned_hgvs):
            cleaned_hgvs = cleaned_hgvs.replace(":", ":g.")

        if hgvs_string != cleaned_hgvs:
            # WARNING, THIS GETS IGNORED IN SEARCH, calling code just checks itself if there's been any difference
            search_messages.append(f'Cleaned "{hgvs_string}" =>"{cleaned_hgvs}"')

        fixed_hgvs, fixed_messages = self.fix_gene_transcript(cleaned_hgvs)
        if fixed_hgvs != cleaned_hgvs:
            search_messages.extend(fixed_messages)
            cleaned_hgvs = fixed_hgvs
        return cleaned_hgvs, search_messages

    @staticmethod
    def fix_gene_transcript(hgvs_string: str) -> tuple[str, list[str]]:
        """ Fix common case of 'GATA2(NM_032638.5):c.1082G>C' and lower case transcript IDs """

        fixed_messages = []

        try:
            prefix, allele = hgvs_string.split(":")
        except ValueError:
            return hgvs_string, []  # Can't do anything here

        if m := re.match(r"(.+)\((.*)\)", prefix):  # Both provided
            transcript_accession, gene_symbol = m.groups()
        else:
            transcript_accession = prefix
            gene_symbol = None

        transcript_ok = looks_like_transcript(transcript_accession)
        if not transcript_ok and gene_symbol:
            # fix gene/transcript swap and lower case separately to get separate warnings.
            uc_gene = gene_symbol.upper()
            if looks_like_transcript(uc_gene):  # Need to upper here
                old_transcript = transcript_accession
                transcript_accession = gene_symbol
                gene_symbol = old_transcript
                if gene_symbol:
                    fixed_messages.append("Swapped gene/transcript")
                transcript_ok = looks_like_transcript(transcript_accession)
            elif looks_like_hgvs_prefix(uc_gene):
                gene_symbol = uc_gene
                fixed_messages.append("Upper cased HGVS prefix")

        if not transcript_ok:
            if transcript_accession:
                uc_transcript = transcript_accession.upper()
                if looks_like_transcript(uc_transcript):
                    transcript_accession = uc_transcript
                    fixed_messages.append("Upper cased transcript")

        if gene_symbol:
            prefix = f"{transcript_accession}({gene_symbol})"
        else:
            prefix = transcript_accession
        fixed_hgvs_string = f"{prefix}:{allele}"
        return fixed_hgvs_string, fixed_messages

    def get_gene_symbol_if_no_transcript(self, hgvs_string: str) -> Optional[str]:
        """ If HGVS uses gene symbol instead of transcript, return symbol """
        hgvs_variant = self.create_hgvs_variant(hgvs_string)
        return hgvs_variant.get_gene_symbol_if_no_transcript()


def get_hgvs_variant_coordinate(hgvs_string: str, genome_build: GenomeBuild) -> VariantCoordinate:
    """ Convenience method for 1 off HGVS - for batches use HGVSMatcher """
    matcher = HGVSMatcher(genome_build)
    return matcher.get_variant_coordinate(hgvs_string)


def get_hgvs_variant(hgvs_name: str, genome_build: GenomeBuild) -> Optional[Variant]:
    """ Convenience method for 1 off HGVS - for batches use HGVSMatcher """
    vc = get_hgvs_variant_coordinate(hgvs_name, genome_build)
    try:
        variant = Variant.get_from_variant_coordinate(vc, genome_build)
    except Variant.DoesNotExist:
        variant = None
    return variant
