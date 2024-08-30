import re
from dataclasses import dataclass
from functools import cached_property
from typing import Dict, Any, Mapping, Optional, Union, List, TypedDict

from django.conf import settings

from annotation.models import CitationFetchRequest
from annotation.models.models_citations import CitationFetchResponse
from classification.criteria_strengths import CriteriaStrength, CriteriaStrengths
from classification.enums import SpecialEKeys
from genes.hgvs import CHGVS, PHGVS
from library.log_utils import report_message
from library.utils import empty_to_none
from snpdb.models import GenomeBuild, GenomeBuildPatchVersion


class VCDbRefDict(TypedDict, total=False):
    id: str
    db: str
    idx: Union[str, int]
    url: str
    summary: str
    internal_id: int


class VCValidation(TypedDict, total=True):
    severity: str
    code: str
    message: str
    options: Optional[Any]


class VCBlobDict(TypedDict, total=False):
    value: Any
    note: str
    explain: str
    immutable: str
    db_refs: List[VCDbRefDict]
    validation: List[VCValidation]


class SomaticValueDict(TypedDict):
    somatic_clinical_significance: str
    amp_level: Optional[str]


@dataclass(frozen=True)
class SomaticClinicalSignificanceValue:
    tier_level: str
    amp_level: Optional[str] = None

    @property
    def without_amp_level(self) -> 'SopmaticClinicalSignificanceValue':
        return SomaticClinicalSignificanceValue(tier_level=self.tier_level)

    @property
    def sort_value(self) -> Optional[int]:
        if sort_value := _SOMATIC_CLINICAL_SIGNIFICANCE_SORT_VALUES.get(self):
            return sort_value
        elif self.amp_level:
            if sort_value := _SOMATIC_CLINICAL_SIGNIFICANCE_SORT_VALUES.get(self.without_amp_level):
                return sort_value
        return None

    def __lt__(self, other):
        return self.sort_value or 0 < other.sort_value or 0

    def as_json(self):
        return {
            "somatic_clinical_significance": self.tier_level,
            "amp_level": self.level
        }

    @property
    def as_str(self):
        if self.amp_level:
            return f"{self.tier_level}|{self.amp_level}"
        else:
            return self.tier_level

    @staticmethod
    def from_str(value: str):
        parts = value.split("|")
        if len(parts > 1):
            return SomaticClinicalSignificanceValue(parts[0], parts[1])
        else:
            return SomaticClinicalSignificanceValue(parts[0])



_SOMATIC_CLINICAL_SIGNIFICANCE_SORT_VALUES = {
    SomaticClinicalSignificanceValue("tier_1", "A"): 9,
    SomaticClinicalSignificanceValue("tier_1", "B"): 8,
    SomaticClinicalSignificanceValue("tier_1"): 7,
    SomaticClinicalSignificanceValue("tier_1_or_2"): 6,
    SomaticClinicalSignificanceValue("tier_2", "C"): 5,
    SomaticClinicalSignificanceValue("tier_2", "D"): 4,
    SomaticClinicalSignificanceValue("tier_2"): 3,
    SomaticClinicalSignificanceValue("tier_3"): 2,
    SomaticClinicalSignificanceValue("tier_4"): 1
}


VCStoreValue = VCBlobDict
VCPatchValue = Union[None, VCStoreValue]
VCStore = Dict[str, VCStoreValue]
VCPatch = Dict[str, VCPatchValue]


class EvidenceMixin:
    """
    For methods common between Classification and ClassificationModification
    Specifically anything that simply needs there to be a dictionary of evidence
    """

    @property
    def _evidence(self) -> VCStore:
        raise NotImplementedError('EvidenceMixin must implement evidence_dict')

    @staticmethod
    def get_optional_value_from(evidence: dict, key: str):
        if evidence is None:
            return None

        blob = evidence.get(key)
        if blob is None:
            return None

        if not isinstance(blob, Mapping):
            return blob

        return blob.get('value')

    def get(self, key: str, default=None):
        value = self.get_optional_value_from(self._evidence or {}, key)
        if value is None:
            return default
        return value

    def get_value_list(self, key: str) -> list:
        value = self.get_optional_value_from(self._evidence or {}, key)
        if value is None:
            return []
        elif isinstance(value, list):
            return value
        else:
            return [value]

    def __getitem__(self, key: str):
        value = self.get(key)
        if value is None:
            raise KeyError(f'No value for {key}')
        return value

    def get_genome_build(self) -> GenomeBuild:
        build_name: str
        try:
            build_name = self[SpecialEKeys.GENOME_BUILD]
        except KeyError:
            raise ValueError("Classification does not have a value for genome build")

        try:
            return GenomeBuild.get_name_or_alias(build_name)
        except GenomeBuild.DoesNotExist:
            raise ValueError(f"Unsupported GenomeBuild {build_name}")

    def get_genome_build_opt(self) -> Optional[GenomeBuild]:
        try:
            return self.get_genome_build()
        except ValueError:
            return None

    def get_genome_build_patch_version(self) -> GenomeBuildPatchVersion:
        if build_name := self.get(SpecialEKeys.GENOME_BUILD):
            try:
                return GenomeBuildPatchVersion.get_or_create(build_name)
            except ValueError:
                return GenomeBuildPatchVersion.get_unspecified_patch_version_for(self.get_genome_build())
        else:
            raise ValueError("Classification does not have a value for genome build")

    @cached_property
    def db_refs(self) -> List[VCDbRefDict]:
        all_db_refs = []
        for blob in self._evidence.values():
            db_refs = blob.get('db_refs')
            if db_refs:
                all_db_refs.extend(db_refs)
        return all_db_refs

    def loaded_citations(self) -> CitationFetchResponse:
        return CitationFetchRequest.fetch_all_now(self.db_refs)

    @property
    def is_likely_acmg(self) -> bool:
        if non_standard_res := settings.CLASSIFICATION_NON_ACMG_ASSERTION_METHOD:
            for assertion_method in self.get_value_list(SpecialEKeys.ASSERTION_METHOD):
                for non_standard_re in non_standard_res:
                    if non_standard_re.match(assertion_method):
                        return False
        return True

    def criteria_strengths(self, e_keys: Optional['EvidenceKeyMap'] = None) -> CriteriaStrengths:
        from classification.models import EvidenceKeyMap
        if not e_keys:
            e_keys = EvidenceKeyMap.instance()

        criteria: List[CriteriaStrength] = []
        for ek in e_keys.criteria():
            if strength := self.get(ek.key):
                criteria.append(CriteriaStrength(ek, strength))

        return CriteriaStrengths(strengths=criteria, is_acmg_standard=self.is_likely_acmg)

    def criteria_strength_summary(self, ekeys: Optional['EvidenceKeyMap'] = None, only_acmg: bool = False) -> str:
        strengths = self.criteria_strengths(e_keys=ekeys)
        return strengths.summary_string(acmg_only=only_acmg)

    @cached_property
    def amp_level(self) -> Optional[str]:
        for key, level in SpecialEKeys.AMP_LEVELS_TO_LEVEL.items():
            if self.get(key):
                return level
        return None

    @property
    def somatic_clinical_significance_value(self) -> SomaticClinicalSignificanceValue:
        if somatic_clin_sig := self.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE):
            return SomaticClinicalSignificanceValue(
                tier_level=somatic_clin_sig,
                amp_level=self.amp_level
            )
        return None

    @cached_property
    def c_parts(self) -> CHGVS:
        return CHGVS(self.get(SpecialEKeys.C_HGVS) or "")

    @cached_property
    def p_parts(self) -> PHGVS:
        return PHGVS.parse(self.get(SpecialEKeys.P_HGVS), override_is_confirmed_to=False)

    @property
    def transcript(self) -> Optional[str]:
        """
        Returns the transcript exactly how it was imported (don't bump versions even if it's correct to do so)
        :return:
        """
        if c_hgvs := self.get(SpecialEKeys.C_HGVS):
            parts = CHGVS(c_hgvs)
            if parts.transcript:
                return parts.transcript

        for transcript_key in settings.VARIANT_ANNOTATION_TRANSCRIPT_PREFERENCES:
            transcript_key = {
                'lrg_identifier': SpecialEKeys.LRG_ID,
                'refseq_transcript_accession': SpecialEKeys.REFSEQ_TRANSCRIPT_ID,
                'ensembl_transcript_accession': SpecialEKeys.ENSEMBL_TRANSCRIPT_ID
            }.get(transcript_key)
            if transcript := self.get(transcript_key):
                return transcript

        # note old code would try to fix the transcript version of the record
        # but only if refseq_transcript_id and ensemble were blank
        """
        matcher = HGVSMatcher(self.get_genome_build())
        transcript_id = matcher.get_transcript_id(c_hgvs, transcript_version=True)
        if transcript_id:
            return transcript_id
        """

        return None

    def calc_clinical_significance_choice(self) -> Optional[str]:
        cs = self.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        if not cs:
            return None
        from classification.models import EvidenceKeyMap
        options = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).matched_options(cs)
        if options:
            return options[0].get('vg')

    @staticmethod
    def _clean_key(key):
        key = str(key).lower()
        # Remove all non-word characters (everything except numbers and letters)
        key = re.sub(r"[^\w\s:]", ' ', key).strip()

        # Replace all runs of whitespace with a single underscore
        key = re.sub(r"\s+", '_', key)
        return key

    @staticmethod
    def to_patch(raw: Dict[str, Any]) -> VCPatch:
        """
        Cleans up a dictionary significantly ready for processing.
        Converts keys to numbers and letters, and all whitespace to underscores.
        For each entry in the source dict -
        * Converts Nones to dicts with None for all entries
        * Converts values to dicts with an entry of value (and None entries for note and explain)
        * Tidies Dicts to have all empty to none
        :param raw: A dictionary, presumably from JSON
        :return: A VCStore where all keys are clean str keys and all values are dicts
        """

        from classification.models import EvidenceKeyMap
        keys = EvidenceKeyMap.instance()

        clean: VCPatch = {}
        for key, value_obj in raw.items():
            key = keys.with_namespace_if_required(EvidenceMixin._clean_key(key))
            if key in clean:
                report_message(message=f'Multiple keys have been normalised to {key}',
                               extra_data={'raw_keys': list(raw.keys())},
                               level='warning')

            if value_obj is not None:
                if not isinstance(value_obj, Mapping):
                    value_obj = {"value": empty_to_none(value_obj), "note": None, "explain": None}

                else:
                    for attr_key, attr_value in value_obj.items():
                        value_obj[attr_key] = empty_to_none(attr_value)

            clean[key] = value_obj

        return clean

    @staticmethod
    def patch_with(target: dict, patch: dict, tidy_nones=False):
        """
        Update the evidence with normalised patch values.
        """
        for key, value in patch.items():
            # providing None means to wipe the whole object
            if value is None:
                target.pop(key, None)

            elif key not in target:
                target[key] = value

            else:
                existing = target[key]
                for sub_key, sub_value in value.items():
                    if sub_value is None:
                        existing.pop(sub_key, None)
                    else:
                        existing[sub_key] = sub_value

        if tidy_nones:
            for key in patch.keys():
                if (blob := target.get(key)) and isinstance(blob, dict):
                    for value in blob.values():
                        if value is not None:
                            break
                    else:
                        target.pop(key, None)
