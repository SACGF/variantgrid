# used for validating multiple keys when one changes
import operator
from dataclasses import dataclass
from enum import Enum
from functools import reduce
from typing import List, Set, Optional, Union, Any, Dict, Iterable

from django.contrib.auth.models import User
from django.db.models import Q

from classification.enums import SpecialEKeys
from classification.models.evidence_mixin import VCPatch, VCStore
from flags.models import FlagCollection
from genes.models import GeneSymbol, Gene, Transcript
from library.utils import VarsDict
from snpdb.models import VariantCoordinate, Allele, Lab


@dataclass
class PatchMessage(VarsDict):
    code: str
    message: str
    key: Optional[str]

    def to_json(self):
        json_blob = {"code": self.code, "message": self.message}
        if self.key:
            json_blob["key"] = self.key
        return json_blob


class ClassificationPatchStatus(str, Enum):
    # doesn't capture super subtle actions like publishing when there were outstanding changes
    # but will capture everything we want to know about for a multi-record import
    UNKNOWN = "unknown"
    NEW = "new"
    UPDATE = "update"
    NO_CHANGE = "no-change"
    WITHDRAWN = "withdrawn"
    DELETED = "deleted"
    UN_WITHDRAWN = "un-withdrawn"
    ALREADY_WITHDRAWN = "already-withdrawn"


class ClassificationPatchResponse(VarsDict):

    def __init__(self):
        self.warnings: List[PatchMessage] = list()
        self.modified_keys: Set[str] = set()
        self.classification_json: Optional[Dict] = None
        self.internal_error: Optional[Any] = None
        self.withdrawn = None
        self.deleted = None
        self.published = False
        self.saved = False
        self.status = ClassificationPatchStatus.UNKNOWN

    def __iadd__(self, other: 'ClassificationPatchResponse'):
        self.warnings += other.warnings
        self.modified_keys |= other.modified_keys
        self.internal_error = self.internal_error or other.internal_error
        self.withdrawn = self.withdrawn or other.withdrawn
        self.deleted = self.deleted or other.deleted
        self.published = self.published or other.published
        self.saved = self.saved or other.saved
        # classification json is a hard field to merge, just take the 2nd?
        self.classification_json = other.classification_json or self.classification_json
        # status is a hard field to merge
        return self

    def append_warning(self, code: str, message: str, key: Optional[str] = None):
        self.warnings.append(PatchMessage(code=code, message=message, key=key))

    def to_json(self):
        # we have two modes, just from the patch, and from the API where we send a bunch of data
        # a lot of this is legacy, bringing it all into ClassificationPatchResponse will hopefully let
        # us consolidate
        if self.internal_error:
            return {"internal_error": str(self.internal_error)}

        response_json = {}
        if self.classification_json:
            response_json = self.classification_json.copy()

        if self.warnings:
            response_json["patch_messages"] = self.warnings
        if self.deleted:
            response_json["deleted"] = True
        if self.withdrawn:
            response_json["withdrawn"] = True
        return response_json


class ValidationMerger:
    """
    Used by multi-value validation to say what errors can be removed (as they've been re-evaluated)
    and what ones need to be added
    """

    @staticmethod
    def union_send_responses(tuples: List) -> 'ValidationMerger':
        validations = [v[1] for v in tuples]
        merged = ValidationMerger()
        for vm in validations:
            if vm:
                merged.remove_codes.update(vm.remove_codes)
                merged.flat_messages = merged.flat_messages + vm.flat_messages
        return merged

    def __init__(self):
        self.remove_codes = {}
        self.flat_messages = []

    def tested(self, keys: Iterable[str], codes: Iterable[str]):
        if isinstance(codes, str):
            for key in keys:
                self.remove_codes[key] = (self.remove_codes.get(key, [])) + [codes]
        else:
            for code in codes:
                self.tested(keys, code)

    #  FIXME fix terminology of message/validation
    def add_message(self, key: str, code: str, severity: str, message: str, options=None, link=None):
        valid_dict = {
            "key": key,
            "severity": severity,
            "code": code,
            "message": message
        }
        if options:
            valid_dict['options'] = options
        if link:
            valid_dict['link'] = link

        self.flat_messages.append(valid_dict)

    @property
    def messages_dict(self):
        message_d = {}
        for message in self.flat_messages:
            copied = message.copy()
            key = copied.pop('key')
            message_d[key] = message_d.get(key, []) + [message]
        return message_d

    def apply(self, patch: dict, evidence: dict):
        for key, codes in self.remove_codes.items():
            entry = None
            if key in patch:
                entry = patch.get(key)
                if entry is None:
                    # we're just removing values, don't need to make entry to remove it
                    continue
            elif key in evidence:
                existing = evidence[key]
                if existing and any(v.get('code') in codes for v in existing.get('validation', [])):
                    entry = existing
                    patch[key] = entry

            if entry and 'validation' in entry:
                # filter out errors that we've tested for
                entry['validation'] = [v for v in entry['validation'] if v.get('code') not in codes]
                if len(entry['validation']) == 0:
                    entry.pop('validation')

        for key, messages in self.messages_dict.items():
            if key in patch:
                entry = patch[key]
                if entry is None:
                    entry = {}
                    patch[key] = entry
            else:
                entry = evidence.get(key, {})
                patch[key] = entry

            entry['validation'] = entry.get('validation', []) + messages


class ClassificationJsonParams:
    """
    When converting Classifications to JSON use these parameters to indicate how.
    """
    def __init__(self,
                 current_user: User,
                 include_data=False,
                 version: Optional[Union['ClassificationModification', float]] = None,
                 flatten=False,
                 include_lab_config=False,
                 include_messages=True,
                 strip_complicated=False,
                 api_version=1,
                 hardcode_extra_data: Dict = None,
                 fix_data_types=False):
        """
        :param current_user: The user who will be consuming this data
        :param include_data: Include all the evidence for this classification (typically True for a GET and False for a POST)
        :param version: Which version of the data to show, if omitted will be latest version
        :param flatten: If True, instead of "c_hgvs": {"value":"x","note":"Oops"} we'd get "c_hgvs.value":"x","c_hgvs.note":"Oops"
        :param include_lab_config: If True, return the lab's EvidenceKey overrides
        :param include_messages: If True, return any validation messages
        :param strip_complicated: If True, note, explain, db_refs will be ignored - used to get to the raw data
        :param api_version: 1 (typically for our own forms that haven't upgraded) 2 where the format of the data more closely mimics what you upload in a VC import
        :param hardcode_extra_data: if provided, will inject itself into the resulting JSON
        :param fix_data_types: if an evidence key has changed to/from a multiselect
        """
        self.current_user = current_user
        self.include_data = include_data
        self.version = version
        self.flatten = flatten
        self.include_lab_config = include_lab_config
        self.include_messages = include_messages
        self.strip_complicated = strip_complicated
        self.api_version = api_version
        self.hardcode_extra_data = hardcode_extra_data
        self.fix_data_types = fix_data_types


class VariantCoordinateFromEvidence:
    """
    Used to extract the variant coordinates that we can then use for variant matching
    Also keeps track of which values we used or which values we couldn't for diagnostic purposes
    """
    def __init__(self, classification: 'Classification'):
        """
        :param classification: The classification we're extracting data from
        """
        from classification.models import classification_flag_types
        self.variant_coordinate = None
        self.messages = []
        self.matching_flag = classification.flag_collection_safe.get_flag_of_type(
            flag_type=classification_flag_types.matching_variant_flag, open_only=True)
        try:
            genome_build = classification.get_genome_build()
            self.genome_build_str = genome_build.name
        except:
            self.genome_build_str = 'No genome build'

    def record(self, value: str, variant_coordinate: Optional[VariantCoordinate] = None,
               message: Optional[str] = None, error: Optional[str] = None) -> None:
        """
        Record what value we're using (or attempted to use) to get the variant coordinates
        :param value: The source value we tried to extract VariantCoordinate from e.g. a c.hgvs, g.hgvs or variant coordinate string
        :param variant_coordinate: The variant coordinate we processed from the value (or None if we couldn't)
        :param message: If present, add message
        :param error: If present, indicates the reason why we couldn't extract a variant_coordinate from the value
        """
        if variant_coordinate:
            self.variant_coordinate = variant_coordinate

        if error:
            self.messages.append(f'Attempted to match {self.genome_build_str} {value}, error = {error}')
        elif variant_coordinate:
            ref = variant_coordinate.ref
            if len(ref) > 50:
                ref = ref[:50] + '...'
            self.messages.append(f'Matching on {self.genome_build_str} {value} resolved to {variant_coordinate.chrom}:{variant_coordinate.pos} {ref}->{variant_coordinate.alt}')
        else:
            self.messages.append(f'Attempted to match {self.genome_build_str} {value}, could not derive coordinate')

        if message:
            self.messages.append(message)

    def report(self) -> None:
        """
        Record all information about the VariantMatching evidence on the VariantMatching flag of the Classification
        """
        if self.matching_flag:
            if not self.messages:
                self.messages.append('No variant related values to match on')
            message = '\n\n'.join(self.messages)
            self.matching_flag.flag_action(comment=message)


class PatchMeta:

    def __init__(self, patch: VCPatch, existing: VCStore, revalidate_all: bool = False):
        """
        :param patch: The patch of new values that's about to overwrite existing
        :param existing: The complete current state of data
        :param revalidate_all: If true, methods that see if things have changed will always return True
        """
        self.patch = patch
        self.existing = existing
        self.revalidate_all = revalidate_all

        if patch is None:
            raise ValueError('patch must have a value')
        if existing is None:
            raise ValueError('existing must have a value')

        self.modified_keys = set(self.patch.keys())

    def patch_value(self, key: str, value: Any):
        if key in self.patch:
            if self.patch[key] is None:
                if value is None:
                    return
                self.patch[key] = {"value": value}

            self.patch[key]['value'] = value
        else:
            self.patch[key] = {'value': value}
        self.modified_keys.add(key)

    def remove_patch_value(self, key: str):
        if key in self.patch or key in self.existing:
            self.patch[key] = {'value': None, 'explain': None, 'note': None, 'immutable': None}

    def get(self, key: str, fallback_existing=True):
        """
        :param key: The str evidence key
        :param fallback_existing: If true, will get the value from existing if it's not in the patch
        :return: The value from patch for the key
        """
        # do or with {} as some dicts might have None as values so providing a default value will still return None
        if key in self.patch:
            return (self.patch.get(key) or {}).get('value')
        if fallback_existing:
            return (self.existing.get(key, {}) or {}).get('value')
        return None

    def is_modified(self, key: str) -> bool:
        """
        Returns if the key should be revalidated
        :param key: the str evidence key
        :return: True if the key has changed in the patch or if we're revalidating_all
        """
        return self.revalidate_all or key in self.modified_keys

    def intersection_modified(self, key_set: Set[str]) -> Set[str]:
        """
        Performs is_modified but over a set
        :param key_set: The set of str evidence keys we're asking about
        :return: The subset of keys that returned True for is_modified
        """
        if self.revalidate_all:
            return key_set
        return key_set.intersection(self.modified_keys)


class UserClassificationStats:
    def __init__(self, user: User):
        self.user = user

    @property
    def issue_count(self) -> int:
        from classification.models import Classification

        return FlagCollection.filter_for_open_flags(
            Classification.filter_for_user(user=self.user)
        ).filter(lab__in=Lab.valid_labs_qs(self.user, admin_check=True)).order_by('-created').exclude(withdrawn=True).count()


def classification_gene_symbol_filter(gene_symbol: Union[str, GeneSymbol]) -> Optional[Q]:

    # We want to filter using the genes set via variant annotation
    genes: Optional[Iterable[Gene]] = None
    symbols: Optional[Iterable[str]] = None
    if gene_symbol := GeneSymbol.cast(gene_symbol):
        genes = gene_symbol.alias_meta.genes
        symbols = gene_symbol.alias_meta.alias_symbol_strs

        if genes:
            allele_qs = Allele.objects.filter(variantallele__variant__variantannotation__gene__in=genes)
            match_gene = Q(classification__variant__variantallele__allele__in=allele_qs)
            evidence_q_list = list()

            for symbol in symbols:
                evidence_q_list.append(Q(published_evidence__gene_symbol__value__iexact=symbol))

            t_qs = Transcript.objects.filter(transcriptversion__gene_version__gene__in=genes).distinct()
            evidence_q_list.append(Q(classification__transcript_version_grch37__transcript__in=t_qs))
            evidence_q_list.append(Q(classification__transcript_version_grch38__transcript__in=t_qs))
            """
            for transcript_id, annotation_consortium in t_qs.values_list("identifier", "annotation_consortium"):
                e_key = SpecialEKeys.ANNOTATION_CONSORTIUM_KEYS[annotation_consortium]
                # We want to match transcript versions (ie X should match X, X.1, X.2 etc)
                # This is quick but wrong - startswith=NM_001099 wrongly matches NM_00109911
                q_transcript_startswith = Q(**{f"published_evidence__{e_key}__value__startswith": transcript_id})
                # This regex correctly handles transcripts but is very slow (10x previous query)
                regex = rf"^{transcript_id}(\.\d+)?$"
                q_regex = Q(**{f"published_evidence__{e_key}__value__regex": regex})
                # Run both quick startwith and regex queries (only slightly slower than original startswith)
                evidence_q_list.append(q_transcript_startswith & q_regex)
            """
            match_evidence = reduce(operator.or_, evidence_q_list)
            return match_gene | match_evidence

    return None
