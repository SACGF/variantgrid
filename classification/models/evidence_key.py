import math
import re
from collections import defaultdict
from copy import deepcopy
from dataclasses import field
from enum import Enum
from functools import cached_property
from typing import Any, List, Optional, Dict, Iterable, Mapping, Union, Set, TypedDict, cast

import pydantic
from django.conf import settings
from django.db import models
from django.db.models.deletion import SET_NULL
from django_extensions.db.models import TimeStampedModel

from classification.enums import CriteriaEvaluation, SubmissionSource, SpecialEKeys
from classification.enums.classification_enums import EvidenceCategory, \
    EvidenceKeyValueType, ShareLevel
from classification.models.evidence_mixin import VCBlobDict, VCPatchValue, VCPatch, VCDbRefDict, EvidenceMixin
from library.cache import timed_cache
from library.utils import empty_to_none, strip_json, first
from snpdb.models import VariantGridColumn

CLASSIFICATION_VALUE_TOLERANCE = 0.00000001
"""
The amount a float can change, and for VG to not consider it a change (this is to stop churning changes through
caused by floating point rounding issues)
"""


class EvidenceKeyOption(TypedDict):

    key: str
    """
    The key for the option (as what will be sorted in the DB)
    """

    label: Optional[str]
    """
    The label to display for the option
    """

    default: Optional[bool]
    """
    Is this the default option - only applies to (ACMG) criteria (e.g. BA1's default is BA, PM2 is PM)
    """

    override: Optional[bool]
    """
    Is this considered an override strength - only applies to criteria (ACMG) criteria - not default and not "not met"
    """

    bucket: Optional[int]
    """
    Only used for clinical_significance, what discordant bucket does each value fall into
    """

    namespace: Optional[str]
    """
    If a namespace is provided, the option will only be enabled if that namespace is enabled
    """


class EvidenceKeyOverrides(pydantic.BaseModel):
    evidence_key_config: dict[str, dict[str, Any]] = pydantic.Field(default_factory=dict)
    """
    Entries should confirm to the EvidenceKey JSON but only the elements that require overriding
    """

    namespaces: set[str] = pydantic.Field(default_factory=set)

    def to_json(self):
        data = self.evidence_key_config.copy()
        data["namespaces"] = list(sorted(self.namespaces))
        return data

    @staticmethod
    def from_dict(config_dict: Optional[dict[str, Any]]) -> 'EvidenceKeyOverrides':
        """
        The dictionary is to have keys that match evidence keys.key
        and values of either true/false (for visibility) or a JSON representation
        of parts of an evidence key to override those attributes.

        If there's a key called "namespaces" it's expected to be a list of strings
        these strings should match the namespace of some evidence keys
        :param config_dict: A diction of evidence key overrides and namespaces
        :return: A well-structured EvidenceKeyOverrides
        """
        if not config_dict:
            return EvidenceKeyOverrides()

        evidence_key_config: dict[str, dict[str, Any]] = {}
        namespaces: set[str] = set()

        for key, value in config_dict.items():
            if key == "namespaces":
                namespaces = set(value)
                continue
            if isinstance(value, bool):
                value = {"hide": not value}
            elif not isinstance(value, dict):
                raise ValueError(f"Received illegal value for classification config: {key}: {value}")
            evidence_key_config[key] = value

        return EvidenceKeyOverrides(
            evidence_key_config=evidence_key_config,
            namespaces=namespaces
        )

    @staticmethod
    def merge(*args: 'EvidenceKeyOverrides') -> 'EvidenceKeyOverrides':
        """
        Merges multiple ClassificationConfigs into one. Provide them in increasing priority.
        Note the evidence key attributes will overwrite each other (only if the attributes
        match, different attributes will be merged).
        In addition. namespaces will be unioned together.
        """
        evidence_key_config: dict[str, dict[str, Any]] = {}
        namespaces: set[str] = set()

        for config in reversed(args):
            for key, value in config.evidence_key_config.items():
                existing = evidence_key_config.get(key)
                if existing:
                    for sub_key, sub_value in value.items():
                        existing[sub_key] = sub_value
                else:
                    evidence_key_config[key] = value.copy()
            namespaces |= config.namespaces

        return EvidenceKeyOverrides(
            evidence_key_config=evidence_key_config,
            namespaces=namespaces
        )


class EvidenceKey(TimeStampedModel):
    key = models.TextField(primary_key=True)
    mandatory = models.BooleanField(default=False)

    max_share_level = models.CharField(max_length=16, choices=ShareLevel.choices(), default='logged_in_users')
    """
    max_share_level restricts sharing on an evidence-key level e.g. Public allows evidence
    key level sharing to anon users or exporting to external systems, while setting to
    INSTITUTION restricts visibility of a particular field even if the classification is public
    """

    @property
    def max_share_level_enum(self) -> ShareLevel:
        return ShareLevel(self.max_share_level)

    order = models.IntegerField(default=0)
    """
    Order within the section (if tied on order within a section, sorted by label)
    """

    label = models.TextField(null=True, blank=True)
    sub_label = models.TextField(null=True, blank=True)
    description = models.TextField(null=True, blank=True)
    examples = models.JSONField(null=True, blank=True)
    options = models.JSONField(null=True, blank=True)
    see = models.TextField(null=True, blank=True)
    """
    Primary URL to describe EvidenceKey, might be worth deprecating in favour of URLs in description
    """
    evidence_category = models.CharField(max_length=3, choices=EvidenceCategory.CHOICES, null=False, blank=False)
    value_type = models.CharField(max_length=1, choices=EvidenceKeyValueType.CHOICES, null=False, blank=True, default=EvidenceKeyValueType.FREE_ENTRY)

    # TODO rename to crit_default_evaluation
    default_crit_evaluation = models.TextField(max_length=3, null=True, blank=True, choices=CriteriaEvaluation.CHOICES)
    crit_allows_override_strengths = models.BooleanField(default=False, null=False, blank=True)
    crit_uses_points = models.BooleanField(default=False, null=False, blank=True)

    allow_custom_values = models.BooleanField(default=False, null=False, blank=True)
    hide = models.BooleanField(default=False, null=False, blank=True)

    namespace_overrides = models.JSONField(null=True, blank=True)

    immutable = models.BooleanField(default=False, null=False, blank=True)

    copy_consensus = models.BooleanField(default=True, null=False, blank=True)

    variantgrid_column = models.ForeignKey(VariantGridColumn, blank=True, null=True, on_delete=SET_NULL)
    """
    If provided, column is auto-populated from annotation data
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.is_dummy = False
        self.default_value = None  # only set by org/lab overrides
        self.exclude_namespace = False  # if the key was in a namespace not used by the lab

    def redcap_key(self, suffix: int = None, is_note=False):
        """
        redcap warns if variables go over 26 characters for a key
        BUT we get conflicts if we shorten them... so leave it as it is for now
        (needed to shorten to 20 to deal with suffixes and prefixes)
        if len(lower_key) > 20:
            lower_key = lower_key[:20]
        """
        suffix_str = ''
        if suffix:
            if is_note:
                suffix_str = '_n%s' % str(suffix)
            else:
                suffix_str = '_%s' % str(suffix)

        lower_key = self.key.lower().replace(' ', '_')
        if re.match('^[0-9].*', lower_key):
            lower_key = 'x' + lower_key

        return lower_key + suffix_str

    def matched_options(self, normal_value_obj) -> List[Dict[str, str]]:
        """
        Given a value (or possibly list of values) generate or find the options that match.
        e.g., for the value ["maternal","xsdfdwerew"] for variant inheritance we'd get
        [{"key":"maternal", "value":"Maternal"}, {"key":"xsdfdwerew", "label":"xsdfdwerew"}]
        See in the 2nd example that we still try to provide an option for unmatched.
        :param normal_value_obj: A value or dict with a key value
        :return: A list of options, guaranteed to be the same length as normal_value_obj if that is a list, otherwise
        will return a length of 1.
        """
        value: Optional[Any]
        if isinstance(normal_value_obj, Mapping):
            value = normal_value_obj.get('value')
        else:
            value = normal_value_obj

        options = self.virtual_options
        if options is not None and len(options) > 0:
            matched_options = []
            if not isinstance(value, list):
                value = [value]
            for val in value:
                match = next((option for option in options if option.get('key') == val), None)
                if match:
                    matched_options.append(match)
                else:
                    matched_options.append({'key': val, 'label': val})

            return matched_options

        return [{'key': value, 'label': value}]

    @staticmethod
    def __special_up_options(options: List[Dict[str, Any]], default: str, allow_overrides=True) -> List[EvidenceKeyOption]:
        """
        Converts a CriteriaEvaluation options (e.g., CriteriaEvaluation.BENIGN_OPTIONS) and a default strength (e.g., BM)
        to a list of evidence keys with default and override populated appropriately
        """
        use_options = []
        for entry in options:
            if entry.get('key') == default:
                entry = entry.copy()
                entry['default'] = True
            elif entry.get('key') == 'NM' or entry.get('key') == 'NA':
                pass
            elif allow_overrides:
                entry = entry.copy()
                entry['override'] = True
            else:
                continue

            use_options.append(entry)
        return use_options

    def _virtual_options_for_criteria(self):
        default_strengths = None
        default_crit = self.default_crit_evaluation

        for criteria in [CriteriaEvaluation.BENIGN_OPTIONS,
                         CriteriaEvaluation.NEUTRAL_OPTIONS,
                         CriteriaEvaluation.PATHOGENIC_OPTIONS]:
            if default_crit in (o.get('key') for o in criteria):
                default_strengths = criteria
                break

        if not default_strengths:
            return []
        else:
            use_options = []
            for entry in criteria:
                entry = entry.copy()
                entry_key = entry.get('key')
                if entry_key == default_crit:
                    entry['default'] = True
                elif entry_key in ('NM', 'NA'):
                    pass
                elif self.crit_allows_override_strengths:
                    entry['override'] = True
                else:
                    continue

                if self.crit_uses_points:
                    if points := CriteriaEvaluation.POINTS.get(entry_key):
                        entry["label"] = f"Met: {points} Point{'s' if points > 1 else ''}"
                if entry_key == CriteriaEvaluation.PATHOGENIC_UNSPECIFIED and "horak" in self.key:
                    # TODO hardcoded support for horak saying ONCOGENIC instead of PATHOGENIC, might be a better way to do this
                    entry["label"] = "Oncogenic Unspecified Strength"

                use_options.append(entry)
            return use_options

    @property
    def virtual_options(self) -> Optional[List[EvidenceKeyOption]]:
        if self.options:
            return cast(List[EvidenceKeyOption], self.options)
        if self.value_type == EvidenceKeyValueType.CRITERIA:
            return self._virtual_options_for_criteria()

        return None

    @cached_property
    def option_indexes(self) -> Optional[Dict[str, int]]:
        index_map: Optional[Dict[str, int]] = None
        if options := self.virtual_options:
            index_map = {}
            for index, option in enumerate(options):
                index_map[option.get('key')] = index + 1
        return index_map

    def classification_sorter_value(self, val: Any) -> Union[int, Any]:
        if index_map := self.option_indexes:
            return index_map.get(val, 0)
        return val

    def sort_values(self, values: Iterable) -> list:
        sorter = self.classification_sorter_value
        sorted_list = sorted(values, key=lambda x: (sorter(x), x))
        return sorted_list

    def classification_sorter(self, evidence: Dict[str, Any]) -> Union[int, Any]:
        """
        Provide .classification_sorter as a Callable[dict aka ClassificationData] -> sortable
        """
        val = evidence.get(self.key)
        return self.classification_sorter_value(val)

    def validate(self):
        if options := self.options:
            for option in options:
                option_key = option.get('key')
                if option_key is not None:
                    if ' ' in option_key:
                        raise ValueError(f'{self.key} option with space in it "{option_key}"')

    @property
    def option_dictionary(self) -> Dict[str, str]:
        if options := self.virtual_options:
            return {x.get('key'): x.get('label') for x in options}
        else:
            return {}

    def option_dictionary_property(self, prop: str) -> Dict[str, Any]:
        if options := self.virtual_options:
            return {x.get('key'): x.get(prop) for x in options if prop in x}
        else:
            return {}

    @staticmethod
    def dummy_key(key: str) -> 'EvidenceKey':
        """
        For use when an EvidenceKey (that doesn't exist) is referenced, and you don't want None exceptions
        """
        dummy = EvidenceKey()
        dummy.key = key
        dummy.value_type = EvidenceKeyValueType.FREE_ENTRY
        dummy.is_dummy = True
        return dummy

    @cached_property
    def namespace(self) -> Optional[str]:
        try:
            return self.key[0:self.key.index(':')]
        except:
            return None

    @property
    def without_namespace(self) -> str:
        try:
            return self.key[self.key.index(':')+1:]
        except:
            return self.key

    @property
    def pretty_label(self) -> str:
        if empty_to_none(self.label):
            return self.label
        return EvidenceKey.pretty_label_from_string(self.key)

    @staticmethod
    def pretty_label_from_string(key: str) -> str:
        if re.search("^([A-Z|a-z])[_-]", key):
            key = key[:1] + '-' + key[2:]

        return key[:1].upper() + key[1:].replace('_', ' ')

    def pretty_value_from(self, evidence: EvidenceMixin, empty_value=None):
        raw_value = evidence.get(self.key)
        if raw_value is None:
            return empty_value
        return self.pretty_value(raw_value)

    def pretty_value(self, normal_value_obj: Any, dash_for_none: bool = False) -> Optional[str]:
        """
        :param normal_value_obj: The blob for the evidence key, e.g. {"value":x} or just x
        :param dash_for_none: If the value obj doesn't contain a value, return "-" if no value is selected
        :return
        """
        value: Any
        if isinstance(normal_value_obj, Mapping):
            value = normal_value_obj.get('value')
        else:
            value = normal_value_obj

        if options := self.virtual_options:
            if not isinstance(value, list):
                value = [value]
            str_values = []
            for val in value:
                matched_option: Dict
                if val == '' or val is None:
                    matched_option = next((option for option in options if option.get('key') is None or option.get('key') == ''), None)
                else:
                    matched_option = next((option for option in options if option.get('key') == val), None)

                if matched_option:
                    part_value = matched_option.get('label') or EvidenceKey.pretty_label_from_string(matched_option.get('key'))
                else:
                    part_value = val
                if part_value is not None:
                    str_values.append(part_value)

            if len(str_values) == 0:
                return '-' if dash_for_none else None

            return ', '.join(str_values)

        if (value == '' or value is None) and dash_for_none:
            return '-'

        return value

    @property
    def striped_options(self):
        options = self.virtual_options
        if not options:
            return None
        return [{
            'key': o.get('key'),
            'label': o.get('label'),
            'default': o.get('default'),
            'override': o.get('override'),
            'namespace': o.get('namespace')
        } for o in options]

    def to_json(self) -> dict:
        return strip_json({
            'key': self.key,
            'allow_custom_values': self.allow_custom_values,
            'namespace_overrides': self.namespace_overrides,
            'order': self.order,
            'label': self.label,
            'sub_label': self.sub_label,
            'description': self.description,
            'evidence_category': self.evidence_category,
            'value_type': self.value_type,
            'options': self.striped_options,
            'examples': self.examples,
            'hide': self.hide,
            'mandatory': self.mandatory,
            'see': self.see
        })

    @staticmethod
    def get_value(blob: Union[Mapping, Any]):
        """ Returns value from an Ekey blob (can be {'value': VAL} or just VAL) """
        if isinstance(blob, Mapping):
            return blob.get('value')
        return blob

    def __str__(self):
        if self.label and self.label != self.key:
            name = f"{self.key}: {self.label}"
        else:
            name = self.key
        return name

    @property
    def is_vital_key(self) -> bool:
        if settings.CLASSIFICATION_DOWNLOADABLE_FIELDS == "*":
            return True
        else:
            return self.key in settings.CLASSIFICATION_DOWNLOADABLE_FIELDS


class EvidenceKeyMap:

    @staticmethod
    def cached_key(key: str) -> EvidenceKey:
        return EvidenceKeyMap.instance().get(key)

    @staticmethod
    def pretty_value_for(item, key: str) -> str:
        return EvidenceKeyMap.cached_key(key).pretty_value(item.get(key))

    @staticmethod
    def __ordered_keys() -> List[EvidenceKey]:
        # sort in code (rather than sql) as pretty_label isn't available normally
        key_entries = list(EvidenceKey.objects.all())
        key_entries.sort(key=lambda k: (k.order, k.pretty_label.lower()))

        keys_by_category = {}
        for key_entry in key_entries:
            category_key_entries = keys_by_category.get(key_entry.evidence_category, [])
            category_key_entries.append(key_entry)
            keys_by_category[key_entry.evidence_category] = category_key_entries

        ordered_by_category_key_entries = []
        for categoryTuple in list(EvidenceCategory.CHOICES):
            category = categoryTuple[0]
            category_key_entries = keys_by_category.pop(category, [])
            ordered_by_category_key_entries.extend(category_key_entries)

        return ordered_by_category_key_entries

    @staticmethod
    @timed_cache(ttl=60)
    def instance() -> 'EvidenceKeyMap':
        return EvidenceKeyMap(
            ordered_keys=EvidenceKeyMap.__ordered_keys()
        )

    @staticmethod
    def cached() -> 'EvidenceKeyMap':
        return EvidenceKeyMap.instance()

    def __init__(self, ordered_keys: List[EvidenceKey], config: Optional[EvidenceKeyOverrides] = None):
        if not config:
            config = EvidenceKeyOverrides()
        self._ordered_keys = ordered_keys
        self.key_dict: dict[str, EvidenceKey] = {}
        self.all_keys: list[EvidenceKey] = []
        self.un_namespaced: dict[str, list[EvidenceKey]] = defaultdict(list)

        for key in ordered_keys:
            if key.namespace:
                self.un_namespaced[key.without_namespace].append(key)

            in_namespace = not key.namespace or key.namespace in config.namespaces
            override_config = config.evidence_key_config.get(key.key)
            if not in_namespace or override_config:
                # keep ordered_keys untouched, so we can copy them again with different configs
                key = deepcopy(key)
                if not in_namespace:
                    key.exclude_namespace = True
                if override_config:
                    for override_attr, override_val in override_config.items():
                        setattr(key, override_attr, override_val)
            self.key_dict[key.key] = key
            self.all_keys.append(key)

    def with_namespace_if_required(self, key: str):
        """
        This is primarily to prefix "acmg:" to evidence keys where the ACMG codes didn't use to have a prefix
        """
        if key not in self.key_dict and key in self.un_namespaced:
            key = first(self.un_namespaced.get(key)).key
        return key

    def without_namespace_if_required(self, key: str):
        """
        This is for the few instances where a prefix has been removed, e.g. "somatic:testing_context" going to "testing_con
        """
        if key not in self.key_dict and ":" in key and key.split(":")[1] in self.key_dict:
            key = key.split(":")[1]
        return key

    def with_overrides(self, config: EvidenceKeyOverrides) -> 'EvidenceKeyMap':
        return EvidenceKeyMap(
            ordered_keys=self._ordered_keys,
            config=config
        )

    def get(self, key: str) -> EvidenceKey:
        if key in self.key_dict:
            return self.key_dict[key]
        return EvidenceKey.dummy_key(key)

    def __getitem__(self, item) -> EvidenceKey:
        return self.get(item)

    def __contains__(self, item):
        return item in self.key_dict

    def immutable(self) -> List[EvidenceKey]:
        return [eKey for eKey in self.all_keys if eKey.immutable]

    def mandatory(self) -> List[EvidenceKey]:
        return [eKey for eKey in self.all_keys if eKey.mandatory]

    def share_level(self, sl: ShareLevel) -> List[EvidenceKey]:
        return [eKey for eKey in self.all_keys if eKey.max_share_level_enum == sl]

    def share_level_and_higher(self, sl: ShareLevel) -> List[EvidenceKey]:
        return [eKey for eKey in self.all_keys if eKey.max_share_level_enum in ShareLevel.same_and_higher(sl)]

    def criteria(self) -> List[EvidenceKey]:
        """
        :return: A list of ALL criteria EvidenceKeys, includes standard ACMG and custom ones with namespaces
        """
        return [eKey for eKey in self.all_keys if eKey.value_type == EvidenceKeyValueType.CRITERIA]

    def vital(self) -> List[EvidenceKey]:
        return [e_key for e_key in self.all_keys if e_key.is_vital_key]

    def acmg_criteria(self) -> List[EvidenceKey]:
        """
        :return: A list of STANDARD ACMG criteria EvidenceKeys
        """
        return self.criteria_for("acmg")

    def criteria_for(self, namespace) -> List[EvidenceKey]:
        crit = [eKey for eKey in self.criteria() if eKey.namespace == namespace]
        crit.sort(key=lambda k: k.pretty_label.lower())
        return crit

    @staticmethod
    @timed_cache(ttl=60)
    def clinical_significance_to_bucket() -> dict[str, int]:
        return EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_dictionary_property("bucket")


class WipeMode(Enum):
    """
    When wiping out a value from a patch... how do we wipe it
    """
    SET_NONE = 'set_none'  # None
    SET_EMPTY = 'set_empty'  # {}
    POP = 'pop'  # remove it from the patch dictionary all together
    ATTRIBUTES_TO_NONE = 'atts_to_none'  # {"value": None, "note": None, "explain": None}


class VCDataCell:
    """
    Representation of the value provided for a single evidence key
    Can exist even when there's no value at all in dictionary
    Changes here will automatically be applied back to the data
    """

    def __init__(self, data: VCPatch, e_key: EvidenceKey):
        self.data = data
        self.e_key = e_key
        self.validate = True

    def __str__(self):
        return f'"{self.e_key.key}": {str(self.raw)}'

    @property
    def _my_data(self) -> VCBlobDict:
        return self.data.get(self.e_key.key) or {}

    def _ensure_my_data(self) -> VCBlobDict:
        existing = self.data.get(self.e_key.key)
        if existing is None:
            existing = {}
            self.data[self.e_key.key] = existing
        return existing

    @property
    def raw(self) -> VCPatchValue:
        """
        What is being stored in the dictionary for this key,
        note both missing from the dictionary or being None in the dictionary will both return None
        """
        return self.data.get(self.e_key.key)

    @raw.setter
    def raw(self, value: VCPatchValue):
        self.data[self.e_key.key] = value

    @property
    def provided(self) -> bool:
        """
        Returns true if this key is in the parent dictionary (as a sub-dictionary or None)
        """
        return self.e_key.key in self.data

    @property
    def immutability(self) -> Optional[SubmissionSource]:
        val = self._my_data.get('immutable')
        if val:
            try:
                return SubmissionSource(val)
            except:
                return SubmissionSource.FORM
        else:
            return None

    @immutability.setter
    def immutability(self, value: SubmissionSource):
        str_value = None
        if value:
            str_value = value.value

        self._ensure_my_data()['immutable'] = str_value

    def pop_immutability(self):
        raw = self.raw
        if raw:
            raw.pop('immutable', None)

    @property
    def has_data(self) -> bool:
        """
        Returns true if there is at least one key in this sub-dictionary
        (even if the value for it is None) e.g. {"value": None} has_data
        """
        return bool(self.raw)

    def wipe(self, mode: WipeMode):
        """
        Clear out the contents of the key
        """
        if mode == WipeMode.SET_NONE:
            self.raw = None
        elif mode == WipeMode.SET_EMPTY:
            self.raw = {}
        elif mode == WipeMode.POP:
            self.data.pop(self.e_key.key, None)
        elif mode == WipeMode.ATTRIBUTES_TO_NONE:
            self.raw = dict({'value': None, 'explain': None, 'note': None})

    def clear_validation(self):
        self.raw.pop('validation', None)

    @property
    def value(self) -> Any:
        return self._my_data.get('value')

    @value.setter
    def value(self, value: Any):
        self._ensure_my_data()['value'] = value

    @property
    def explain(self) -> Optional[str]:
        return self._my_data.get('explain')

    @explain.setter
    def explain(self, value: Optional[str]):
        self._ensure_my_data()['explain'] = value

    @property
    def note(self) -> Optional[str]:
        return self._my_data.get('note')

    @note.setter
    def note(self, value: Optional[str]):
        self._ensure_my_data()['note'] = value

    @property
    def db_refs(self) -> Optional[List[VCDbRefDict]]:
        return self._my_data.get('db_refs')

    @db_refs.setter
    def db_refs(self, db_refs: Optional[List[VCDbRefDict]]):
        self._ensure_my_data()['db_refs'] = db_refs

    def strip_non_client_submission(self):
        """ Normalises the entry to only contain what a user can upload """
        stripped = {}
        raw = self.raw
        for key in ['value', 'note', 'explain']:
            if key in raw:
                stripped[key] = raw[key]
        self.raw = stripped

    def has_change(self, other: 'VCDataCell'):
        raw = self.raw
        other_raw = other.raw
        if raw == other_raw:
            return False
        if raw is None or other_raw is None:
            return True
        for attr in ['value', 'note', 'explain']:
            if raw.get(attr) != other_raw.get(attr):
                return True

        if self.immutability != other.immutability:
            return True

        return False

    def merge_from(self, other: 'VCDataCell'):
        """
        Copy attributes from another cell (presumably with the same key) if we don't
        have a value for that attribute yet
        """
        self_raw = self.raw
        other_raw = other.raw
        for attr in ['value', 'note', 'explain']:
            if attr not in self_raw and attr in other_raw:
                self_raw[attr] = other_raw[attr]

    def __eq__(self, other):
        return self.e_key == other.e_key and \
            self.raw == other.raw

    def __contains__(self, item):
        return item in self._my_data

    def diff(self, dest: Optional['VCDataCell'], ignore_if_omitted: Optional[Set[str]] = None) -> Dict[str, Any]:
        """
        Given two dictionaries, returns only the entries that have changed
        :param dest: Another dict
        :param ignore_if_omitted: If a key name is in here, and if a value is in dest, but not in self, don't make any mention of it,
        otherwise it will be set to None.
        :return: The differences
        """
        diff_dict = {}
        source = self._my_data
        dest = dest._my_data if dest else {}
        for key, value in source.items():
            # no diff in the keys here
            if value is None:
                # value is None in source but has a value other than None in dest
                if key in dest and dest[key] is not None:
                    diff_dict[key] = None
            elif key in dest:
                # key is in both source and dest but with diff values
                dest_value = dest[key]
                if dest_value != value:
                    if isinstance(dest_value, float) and isinstance(value, float):
                        if not math.isclose(dest_value, value, abs_tol=CLASSIFICATION_VALUE_TOLERANCE):
                            diff_dict[key] = value
                    else:
                        diff_dict[key] = value
            else:
                # value is not None and key is not in dest
                diff_dict[key] = value

        for key, value in dest.items():
            if key not in source:
                if ignore_if_omitted and key in ignore_if_omitted:
                    pass
                else:
                    diff_dict[key] = None

        return diff_dict

    def add_validation(self, code: str, severity: str, message: str, options=None):
        """
        Add a validation message to this cell
        """
        if self.validate:
            VCDataDict.add_validation(self._ensure_my_data(), code=code, severity=severity, message=message, options=options)

    def has_validation_code(self, code: str) -> bool:
        """
        Returns true if there's a validation message with the given code
        """
        validations = self._my_data.get('validation', [])
        for validation in validations:
            if validation.get('code') == code:
                return True
        return False


class VCDataDict:
    """
    Represents an entire patch or base set of data for EvidenceKeys
    """

    def __init__(self, data: Dict[str, Any], evidence_keys: EvidenceKeyMap):
        if not isinstance(data, dict):
            raise ValueError('Data must be of type dict')
        self.data = data
        self.e_keys = evidence_keys

    def __str__(self):
        return str(self.data)

    def __getitem__(self, item: Union[str, EvidenceKey]) -> VCDataCell:
        """
        Note this will always return a VCDataCell even if there's no entry for it
        (inspect the data cell to see if there's any corresponding entry in this dict)
        """
        if isinstance(item, str):
            item = self.e_keys.get(item)
        return VCDataCell(data=self.data, e_key=item)

    def __contains__(self, item: Union[str, EvidenceKey]):
        if isinstance(item, EvidenceKey):
            item = item.key
        return item in self.data

    def keys(self) -> Iterable[EvidenceKey]:
        return [self.e_keys.get(key) for key in self.data.keys()]

    def cells(self) -> Iterable[VCDataCell]:
        return (self[key] for key in set(self.keys()))

    def ignore(self, key: str):
        self.data.pop(key, None)

    @staticmethod
    def add_validation(value_obj: VCBlobDict, code: str, severity: str, message: str, options=None):
        validation_key = 'validation'
        validation_array = value_obj.get(validation_key, [])
        value_obj[validation_key] = validation_array

        valid_dict = {
            "severity": severity,
            "code": code,
            "message": message
        }
        if options:
            valid_dict['options'] = options

        validation_array.append(valid_dict)
