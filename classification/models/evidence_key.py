from enum import Enum

from django.db import models
from django.db.models.deletion import SET_NULL
from django_extensions.db.models import TimeStampedModel
from lazy import lazy
from typing import Any, List, Optional, Dict, Iterable, Mapping, Union, Set
import re

from library.cache import timed_cache
from library.utils import empty_to_none
from snpdb.models import VariantGridColumn, Lab
from classification.enums import CriteriaEvaluation, SubmissionSource
from classification.enums.classification_enums import EvidenceCategory, \
    EvidenceKeyValueType, ShareLevel
from classification.json_serialize import strip_json
from snpdb.views.datatable_view import RichColumn


class EvidenceKey(TimeStampedModel):
    """ This is global / visible to everyone """
    key = models.TextField(primary_key=True)
    mandatory = models.BooleanField(default=False)
    # max_share_level restricts sharing on an evidence-key level eg Public allows evidence
    # key level sharing to anon users or exporting to external systems, while setting to
    # INSTITUTION restricts visibility of a particular field even if the classification is public
    max_share_level = models.CharField(max_length=16, choices=ShareLevel.choices(), default='logged_in_users')

    @property
    def max_share_level_enum(self) -> ShareLevel:
        return ShareLevel(self.max_share_level)

    order = models.IntegerField(default=0)
    # TODO provide a list of choices
    label = models.TextField(null=True, blank=True)
    sub_label = models.TextField(null=True, blank=True)
    description = models.TextField(null=True, blank=True)
    examples = models.JSONField(null=True, blank=True)
    options = models.JSONField(null=True, blank=True)
    # Primary URL for key
    see = models.TextField(null=True, blank=True)  # see this (ie a link)
    evidence_category = models.CharField(max_length=3, choices=EvidenceCategory.CHOICES, null=False, blank=False)
    value_type = models.CharField(max_length=1, choices=EvidenceKeyValueType.CHOICES, null=False, blank=True, default=EvidenceKeyValueType.FREE_ENTRY)

    default_crit_evaluation = models.TextField(max_length=3, null=True, blank=True, choices=CriteriaEvaluation.CHOICES)

    allow_custom_values = models.BooleanField(default=False, null=False, blank=True)
    hide = models.BooleanField(default=False, null=False, blank=True)

    immutable = models.BooleanField(default=False, null=False, blank=True)

    copy_consensus = models.BooleanField(default=True, null=False, blank=True)

    # This column ise used to auto-populate the keys
    variantgrid_column = models.ForeignKey(VariantGridColumn, blank=True, null=True, on_delete=SET_NULL)

    def __init__(self, *args, **kwargs):
        super(EvidenceKey, self).__init__(*args, **kwargs)
        self.is_dummy = False
        self.default_value = None  # only set by org/lab overrides
        self.exclude_namespace = False  # if the key was in a namespace not used by the lab

    def redcap_key(self, suffix: int = None, is_note=False):
        suffix_str = ''
        if suffix:
            if is_note:
                suffix_str = '_n%s' % str(suffix)
            else:
                suffix_str = '_%s' % str(suffix)

        lower_key = self.key.lower().replace(' ', '_')
        if re.match('^[0-9].*', lower_key):
            lower_key = 'x' + lower_key

        # redcap warns if variables go over 26 characters for a key
        # BUT we get conflicts if we shorten them... so leave it as it is for now
        # (needed to shorten to 20 to deal with suffixes and prefixes)
        # if len(lower_key) > 20:
        #    lower_key = lower_key[:20]

        return lower_key + suffix_str

    def matched_options(self, normalValueObj) -> List[Dict[str, str]]:
        """
        Given a value (or possibly list of values) generate or find the options that match.
        e.g. for the value ["maternal","xsdfdwerew"] for variant inheritance we'd get
        [{"key":"maternal", "value":"Maternal"}, {"key":"xsdfdwerew", "label":"xsdfdwerew"}]
        With he fist one being the result of a match, and the second on being the result of not being matched
        :param normalValueObj: A value or dict with a key value
        :return: A list of options that matched
        """
        value = None
        if isinstance(normalValueObj, Mapping):
            value = normalValueObj.get('value')
        else:
            value = normalValueObj

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
    def __special_up_options(options, default):
        use_options = []
        for entry in options:
            if entry.get('key') == default:
                entry = entry.copy()
                entry['default'] = True
            elif entry.get('key') == 'NM' or entry.get('key') == 'NA':
                pass
            else:
                entry = entry.copy()
                entry['override'] = True
            use_options.append(entry)
        return use_options

    @property
    def virtual_options(self):
        if self.options:
            return self.options
        if self.value_type == EvidenceKeyValueType.CRITERIA:
            for criteria in [CriteriaEvaluation.BENIGN_OPTIONS,
                             CriteriaEvaluation.NEUTRAL_OPTIONS,
                             CriteriaEvaluation.PATHOGENIC_OPTIONS]:
                if self.default_crit_evaluation in (o.get('key') for o in criteria):
                    return EvidenceKey.__special_up_options(criteria, self.default_crit_evaluation)
            print(f"Warning, could not work out what list {self.default_crit_evaluation} fits into")
            return []

        return None

    def validate(self):
        if self.options:
            for option in self.options:
                option_key = option.get('key')
                if option_key is not None:
                    if ' ' in option_key:
                        raise ValueError(f'{self.key} option with space in it "{option_key}"')

    @property
    def option_dictionary(self):
        return {x.get('key'): x.get('label') for x in self.virtual_options}

    @staticmethod
    def dummy_key(key: str):
        dummy = EvidenceKey()
        dummy.key = key
        dummy.value_type = EvidenceKeyValueType.FREE_ENTRY
        dummy.is_dummy = True
        return dummy

    @lazy
    def namespace(self):
        try:
            return self.key[0:self.key.index(':')]
        except:
            return None

    @property
    def pretty_label(self) -> str:
        if empty_to_none(self.label):
            return self.label
        return EvidenceKey.pretty_label_from_string(self.key)

    @staticmethod
    def pretty_label_from_string(key: str):
        if re.search("^([A-Z|a-z])[_-]", key):
            key = key[:1] + '-' + key[2:]

        return key[:1].upper() + key[1:].replace('_', ' ')

    @staticmethod
    def merge_config(config1: dict, config2: dict) -> dict:
        if config1 is None:
            config1 = {}
        if config2 is None:
            config2 = {}

        merged = config1.copy()
        for key, value in config2.items():
            if key in merged:
                existing = merged.get(key)
                if isinstance(existing, bool) and isinstance(value, bool):
                    merged[key] = value
                else:
                    if isinstance(existing, bool):
                        existing = {'hide': not existing}
                    else:
                        existing = existing.copy()
                    if isinstance(value, bool):
                        value = {'hide': not value}
                    for key2, value2 in value.items():
                        existing[key2] = value2
                    merged[key] = existing
            else:
                merged[key] = value
        return merged

    def special_to_key(self, option_key: str, val):
        if val is None:
            return None
        for option in self.virtual_options:
            if option_key in option and option[option_key] == val:
                return option.get('key')
        #couldnt match
        return val

    def pretty_value(self, normalValueObj: Any, dash_for_none: bool = False) -> Optional[str]:
        value = None
        if isinstance(normalValueObj, Mapping):
            value = normalValueObj.get('value')
        else:
            value = normalValueObj

        if (value == '' or value is None) and dash_for_none:
            return '-'

        options = self.virtual_options
        if options is not None and len(options) > 0:
            if not isinstance(value, list):
                value = [value]
            str_values = []
            for val in value:
                match: Dict
                if val == '' or val is None:
                    match = next((option for option in options if option.get('key') is None or option.get('key') == ''), None)
                else:
                    match = next((option for option in options if option.get('key') == val), None)

                if match:
                    part_value = match.get('label') or EvidenceKey.pretty_label_from_string(match.get('key'))
                else:
                    part_value = val
                if part_value is not None:
                    str_values.append(part_value)

            if len(str_values) == 0:
                return None
            return ', '.join(str_values)
        return value

    @property
    def label_dict(self) -> Dict[str, str]:
        options = self.virtual_options
        if options:
            return {option.get('key'): option.get('label') or EvidenceKey.pretty_label_from_string(option.get('key')) for option in options}
        return {}

    @property
    def striped_options(self):
        options = self.virtual_options
        if not options:
            return None
        return [{
            'key': o.get('key'),
            'label': o.get('label'),
            'default': o.get('default'),
            'override': o.get('override')
        } for o in options]

    def to_json(self) -> dict:
        return strip_json({
            'key': self.key,
            'allow_custom_values': self.allow_custom_values,
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

    def __str__(self):
        if self.label and self.label != self.key:
            name = f"{self.key}: {self.label}"
        else:
            name = self.key
        return name

    def rich_column(self, prefix: str):
        return RichColumn(
            key=f"{prefix}__{self.key}__value",
            orderable=True,
            client_renderer=f"VCTable.evidence_key.bind(null, '{self.key}')",
            label=self.pretty_label
        )


class EvidenceKeyMap:

    @staticmethod
    @timed_cache(ttl=60)
    def cached_key(key: str) -> EvidenceKey:
        return EvidenceKey.objects.get(pk=key)

    @staticmethod
    @timed_cache(ttl=60)
    def cached() -> 'EvidenceKeyMap':
        return EvidenceKeyMap()

    @staticmethod
    @timed_cache(ttl=60)
    def _ordered_keys() -> List[EvidenceKey]:
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

    def __init__(self, lab: Lab = None):
        """
        @param lab If provided the keys will be overridden with org/lab specific config
        """
        self.key_dict = {}
        self.all_keys = EvidenceKeyMap._ordered_keys()
        self.lab = lab
        self.namespaces = {None}

        for evidence_key in self.all_keys:
            self.key_dict[evidence_key.key] = evidence_key

        #update keys to have overridden values
        if lab:
            merged = EvidenceKey.merge_config(lab.classification_config, lab.organization.classification_config)
            if 'namespaces' in merged:
                self.namespaces = {None}.union(set(merged.pop('namespaces', [])))

            if merged:
                for okey, ovalue in merged.items():
                    ekey = self.get(okey)
                    if ovalue is False:
                        ekey.hide = True
                    elif ovalue is True:
                        ekey.hide = False
                    else:
                        for param, paramv in ovalue.items():
                            setattr(ekey, param, paramv)

        for evidence_key in self.all_keys:
            if evidence_key.namespace not in self.namespaces:
                evidence_key.exclude_namespace = True

    def get(self, key: str) -> EvidenceKey:
        if key in self.key_dict:
            return self.key_dict[key]
        return EvidenceKey.dummy_key(key)

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

    def acmg_criteria(self) -> List[EvidenceKey]:
        """
        :return: A list of STANDARD ACMG criteria EvidenceKeys
        """
        acmg_crit = [eKey for eKey in self.criteria() if eKey.namespace is None]
        acmg_crit.sort(key=lambda k: k.pretty_label.lower())
        return acmg_crit


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
    Can existing even when there's no value at all in dictionary and changes
    here will automatically be applied back to the data
    """

    def __init__(self, data: Dict[str, Any], e_key: EvidenceKey):
        self.data = data
        self.e_key = e_key

    def __str__(self):
        return f'"{self.e_key.key}": {str(self.raw)}'

    @property
    def _my_data(self) -> Dict[str, Any]:
        return self.data.get(self.e_key.key) or {}

    def _ensure_my_data(self) -> Dict[str, Any]:
        existing = self.data.get(self.e_key.key)
        if existing is None:
            existing = {}
            self.data[self.e_key.key] = existing
        return existing

    @property
    def raw(self) -> Optional[Dict[str, Any]]:
        """
        What is being stored in the dictionary for this key,
        note both missing from the dictionary or being None in the dictioanry will both return None
        """
        return self.data.get(self.e_key.key)

    @raw.setter
    def raw(self, value):
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
    def db_refs(self) -> Optional[List]:
        return self._my_data.get('db_refs')

    @db_refs.setter
    def db_refs(self, db_refs: Optional[List]):
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

    def diff(self, dest: Optional['VCDataCell'], ignore_if_omitted: Optional[Set[str]] = None) -> Dict[str, Any]:
        """
        Given two dictionaries, returns only the entries that have changed
        :param dest: Another dict
        :param ignore_omitted: If a key name is in here, and if a value is in dest, but not in self, don't make any mention of it,
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
                if dest[key] != value:
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
        VCDataDict.add_validation(self._ensure_my_data(), code=code, severity=severity, message=message, options=options)

    def has_validiation_code(self, code: str) -> bool:
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

    def __init__(self, data: Dict[str, Any]):
        if not isinstance(data, dict):
            raise ValueError('Data must be of type dict')
        self.data = data
        self.e_keys = EvidenceKeyMap.cached()

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
    def add_validation(valueObj: dict, code: str, severity: str, message: str, options=None):
        validation_key = 'validation'
        validation_array = valueObj.get(validation_key, [])
        valueObj[validation_key] = validation_array

        valid_dict = {
            "severity": severity,
            "code": code,
            "message": message
        }
        if options:
            valid_dict['options'] = options

        validation_array.append(valid_dict)
