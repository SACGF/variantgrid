import re
from functools import total_ordering
from typing import Dict, Any, Mapping, Optional, Union, List, TypedDict

from django.conf import settings
from lazy import lazy

from annotation.models import Citation, CitationSource
from classification.enums import SpecialEKeys, CriteriaEvaluation
from genes.hgvs import CHGVS, PHGVS
from library.log_utils import report_message
from library.utils import empty_to_none
from snpdb.models import GenomeBuild


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


VCStoreValue = VCBlobDict
VCPatchValue = Union[None, VCStoreValue]
VCStore = Dict[str, VCStoreValue]
VCPatch = Dict[str, VCPatchValue]


@total_ordering
class CriteriaStrength:

    def __init__(self, ekey: 'EvidenceKey', strength: Optional[str]):
        self.ekey = ekey
        self.strength = strength or ekey.default_crit_evaluation

    @property
    def strength_value(self) -> int:
        try:
            return CriteriaEvaluation.ALL_STRENGTHS.index(self.strength)
        except:
            return 0

    def __str__(self) -> str:

        # just need special handling of X
        def strength_suffix(strength: str):
            if strength == "X":
                return "unspecified"
            elif strength.endswith("X"):
                return strength[0] + "_unspecified"
            return strength

        # Make sure criteria are in camel case so removing spaces still leaves it readable
        pretty_label = self.ekey.pretty_label.replace(" ", "")
        suffix = self.strength

        matches_direction = False
        if not self.ekey.namespace:
            if self.ekey.default_crit_evaluation == self.strength:
                return pretty_label
            criteria_first_letter = self.ekey.key[0].upper()
            matches_direction = criteria_first_letter == suffix[0]

        if matches_direction:
            suffix = suffix[1:]

        return f"{pretty_label}_{strength_suffix(suffix)}"

    def __eq__(self, other) -> bool:
        if not isinstance(other, CriteriaStrength):
            return False
        return self.ekey == other.ekey and self.strength == other.strength

    def __hash__(self):
        return hash(self.ekey) + hash(self.strength)

    def __lt__(self, other) -> bool:
        if not isinstance(other, CriteriaStrength):
            raise ValueError(f'Cannot sort CriteriaStrength and {other}')
        if self.strength_value == other.strength_value:
            return self.ekey.pretty_label < other.ekey.pretty_label
        return self.strength_value < other.strength_value


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

    @lazy
    def db_refs(self) -> List[Dict]:
        all_db_refs = []
        for blob in self._evidence.values():
            db_refs = blob.get('db_refs')
            if db_refs:
                all_db_refs.extend(db_refs)
        return all_db_refs

    @property
    def citations(self) -> List[Citation]:
        """
        Returns the entire list of citations through the evidence
        :return: A list of Citations
        """
        citations = []
        for db_ref in self.db_refs:
            source = CitationSource.CODES.get(db_ref.get('db'))
            if source:
                citation, _ = Citation.objects.get_or_create(citation_source=source, citation_id=db_ref.get('idx'))
                citations.append(citation)
        return citations

    def criteria_strength_summary(self, ekeys: Optional['EvidenceKeyMap'] = None, only_acmg: bool = False) -> str:
        if ekeys is None:
            from classification.models import EvidenceKeyMap
            ekeys = EvidenceKeyMap.instance()

        criteria: List[CriteriaStrength] = []
        for ek in ekeys.criteria():
            if only_acmg and ek.namespace:
                continue
            strength = self.get(ek.key)
            if CriteriaEvaluation.is_met(strength):  # exclude neutral, not met, not applicable
                criteria.append(CriteriaStrength(ek, strength))

        criteria.sort()
        return ", ".join([str(c) for c in criteria])

    @lazy
    def c_parts(self) -> CHGVS:
        return CHGVS(self.get(SpecialEKeys.C_HGVS) or "")

    @lazy
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
        clean: VCPatch = {}
        for key, value_obj in raw.items():
            key = EvidenceMixin._clean_key(key)

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
