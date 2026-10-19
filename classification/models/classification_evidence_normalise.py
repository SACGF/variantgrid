"""
Evidence normalisation: turns a submitted evidence cell into the value we store and attaches validation to it.
normalise_evidence_cell is what Classification.patch_value runs over every patched cell - the shared tidying
(non-breaking spaces, escaped HTML, db_ref extraction, mandatory checks, single item lists), then the per value type
normaliser from VALUE_TYPE_NORMALISERS, then the key specific checks. Select values are matched to their options
through evidence_option_value_parts / match_evidence_option, and evidence_validation_messages collects the
validation messages stored in evidence.
"""
from __future__ import annotations

import re
from collections.abc import Callable, Mapping
from datetime import datetime
from typing import TYPE_CHECKING, Any, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models.query_utils import Q

from annotation.regexes import DbRegexes, db_ref_regexes
from classification.enums import CriteriaEvaluation, SpecialEKeys, SubmissionSource, ValidationCode
from classification.models.classification_patcher import patch_fuzzy_age
from classification.models.evidence_key import (
    EvidenceKey,
    EvidenceKeyValueType,
    VCDataCell,
    WipeMode,
)
from library.utils import cautious_attempt_html_to_text, empty_to_none
from snpdb.models import Sample

if TYPE_CHECKING:
    from classification.models.classification import Classification

_key_to_regex = {
    'db_rs_id': DbRegexes.SNP,
    'uniprot_id': DbRegexes.UNIPROTKB,
    'clinvar_variation_id': DbRegexes.CLINVAR,
    'cosmic_id': DbRegexes.COSMIC,
    'gtr_id': DbRegexes.GTR
}

def evidence_option_export_value(normal_value_obj: Any, key_data: EvidenceKey, export_key: str = 'key') -> Any:
    if options := key_data.virtual_options:
        if normal_value_obj and isinstance(normal_value_obj, Mapping):
            value = normal_value_obj.get('value')
        else:
            value = normal_value_obj

        for option in options:
            if export_key in option:
                export_value = option.get(export_key)
                if option.get('key') == value:
                    return export_value
    return None

def match_evidence_option(options, check_value) -> Optional[dict[str, str]]:
    for option in options:
        option_value = option.get('key')

        if check_value is True:
            if option.get('default'):
                return option
            continue
        if not isinstance(check_value, str):
            check_value = str(check_value)

        aliases = []
        if option_value:
            aliases.append(option_value)
        if option.get('label'):
            aliases.append(option.get('label'))
        # Originally had index as a value you could submit, but this starts to get dangerous for fields
        # that allow custom values, can't distinguish between index 3 and custom value 3
        # if option.get('index') is not None:
        #    aliases.append(str(option.get('index')))
        elif option_value:
            aliases.append(EvidenceKey.pretty_label_from_string(option_value))

        aliases.extend(option.get('aliases', []))

        match = next((alias for alias in aliases if alias.lower() == check_value.lower()), None)
        if match is not None:
            return option
        check_value = check_value.replace('-', '_').replace(' ', '_')
        match = next((alias for alias in aliases if alias.lower() == check_value.lower()), None)
        if match is not None:
            return option
    return None

def evidence_option_value_parts(e_key: EvidenceKey, value: Any) -> list[Any]:
    """ The separate options in a select value as submitted: a string is split on commas (or '|_' for a
        multi-select that uses it), a list is taken as is. Blanks are dropped """
    parts: list[Any]
    if isinstance(value, str):
        if e_key.value_type == EvidenceKeyValueType.MULTISELECT and '|_' in value:
            parts = value.split('|_')
        else:
            parts = value.split(',')
        parts = [p.strip() for p in parts]
    elif isinstance(value, list):
        parts = value
    else:
        parts = [value]
    return [part for part in parts if part is not None and part != ""]

def _normalise_option_values(cell: VCDataCell, values: list[Any]) -> Optional[list[str]]:
    e_key = cell.e_key
    options = e_key.virtual_options or []
    # Do a case-insensitive check for each value against the key and any aliases
    # if there's a match to any of those, normalise back to the key (with the case of the key)
    results: list[str] = []
    # remove duplicates
    values = list(set(values))
    for check_value in values:
        check_value = empty_to_none(check_value)
        if check_value is None:
            continue
        if check_value is True and e_key.default_crit_evaluation:
            results.append(e_key.default_crit_evaluation)
            continue
        if check_value is False and e_key.value_type == EvidenceKeyValueType.CRITERIA:
            results.append(CriteriaEvaluation.NOT_MET)
            continue

        check_value = str(check_value)
        matched_option = match_evidence_option(options, check_value)
        if matched_option:
            results.append(matched_option.get('key'))
            if settings.CLASSIFICATION_REQUIRE_OVERWRITE_NOTE and \
                    matched_option.get('override') and cell.note is None and cell.explain is None:
                cell.add_validation(
                    code=ValidationCode.REQUIRES_NOTE,
                    severity='warning',
                    message='Override value has been selected, but no note provided'
                )
        else:
            # add the value we couldn't match
            results.append(str(check_value))

            if not e_key.allow_custom_values:
                cell.add_validation(
                    code=ValidationCode.INVALID_VALUE,
                    severity='error',
                    message='Illegal value (' + str(check_value) + ')'
                )

    if results:
        results.sort()  # TODO, sort by presentation index
        return results
    return None


# Each normaliser takes the cell and its tidied value, sets cell.value (and any validation) and returns the value
# for the key specific checks in normalise_evidence_cell to work on

def _normalise_user(cell: VCDataCell, value: Any) -> Any:
    # convert users (username or email) to username if possible
    value = str(value)
    if user := User.objects.filter(Q(username=value) | Q(email=value)).first():
        value = user.username
    cell.value = value
    return value


def _normalise_boolean(cell: VCDataCell, value: Any) -> Any:
    # normalise booleans if it's unambiguous what the value is meant to be,
    # otherwise leave the value as is to be caught by validation
    if isinstance(value, list) and len(value) == 1:
        value = value[0]
    if not isinstance(value, bool):
        if value == 1:
            value = True
        elif value == 0:
            value = False
        elif isinstance(value, str):
            value = value.upper()
            if value in ('TRUE', 'Y', 'YES', 'T', '1'):
                value = True
            elif value in ('FALSE', 'N', 'NO', 'F', '0'):
                value = False
            else:
                cell.add_validation(code=ValidationCode.INVALID_VALUE, severity='error',
                                    message="Invalid value for yes/no field (" + str(value) + ")")
    cell.value = value
    return value


def _normalise_options(cell: VCDataCell, value: Any) -> Any:
    # handle selects with multiple values (comma sep string, array etc.)
    # sets the value of cell to the valid options
    e_key = cell.e_key
    value = _normalise_option_values(cell, evidence_option_value_parts(e_key, value))
    if value is not None:
        if e_key.value_type != EvidenceKeyValueType.MULTISELECT:
            value_len = len(value)
            if value_len > 1:
                cell.add_validation(code=ValidationCode.TOO_MANY_VALUES, severity='error',
                                    message="Can only have a single value for a select field")
            elif value_len == 1:
                value = value[0]
    cell.value = value
    return value


def _normalise_text(cell: VCDataCell, value: Any) -> Any:
    # base text tidy
    value = str(value)
    if cell.e_key.value_type == EvidenceKeyValueType.FREE_ENTRY:
        # strip out HTML from single row files to keep our data simple
        # allow HTML in text areas
        if not value.startswith("http"):  # makes beautiful soap angry thinking we want to go to the URL
            value = cautious_attempt_html_to_text(value)

    cell.value = value
    return value


def _normalise_unit(cell: VCDataCell, value: Any) -> Any:
    valid = False
    try:
        value = float(value)
        cell.value = value
        valid = 0 <= value <= 1
    except Exception:
        pass
    if not valid:
        cell.add_validation(code=ValidationCode.INVALID_UNIT, severity='warning',
                            message='Value should be a number between 0 and 1')
    return value


def _normalise_float(cell: VCDataCell, value: Any) -> Any:
    valid = False
    try:
        value = float(value)
        cell.value = value
        valid = True
    except Exception:
        pass
    if not valid:
        cell.add_validation(code=ValidationCode.INVALID_FLOAT, severity='warning',
                            message='Value should be a number')
    return value


def _normalise_integer(cell: VCDataCell, value: Any) -> Any:
    valid = False
    try:
        value = int(value)
        cell.value = value
        valid = True
    except Exception:
        pass
    if not valid:
        cell.add_validation(code=ValidationCode.INVALID_INTEGER, severity='warning',
                            message='Value should be a whole number')
    return value


def _validate_date(cell: VCDataCell, value: Any) -> Any:
    try:
        if value:
            datetime.strptime(value, '%Y-%m-%d')
    except ValueError:
        message = "Invalid date (expected yyyy-mm-dd)"
        cell.add_validation(code=ValidationCode.INVALID_DATE, severity='warning',
                            message=message)
    return value


VALUE_TYPE_NORMALISERS: dict[str, Callable[[VCDataCell, Any], Any]] = {
    EvidenceKeyValueType.USER: _normalise_user,
    EvidenceKeyValueType.BOOLEAN: _normalise_boolean,
    EvidenceKeyValueType.MULTISELECT: _normalise_options,
    EvidenceKeyValueType.SELECT: _normalise_options,
    EvidenceKeyValueType.CRITERIA: _normalise_options,
    EvidenceKeyValueType.FREE_ENTRY: _normalise_text,
    EvidenceKeyValueType.TEXT_AREA: _normalise_text,
    EvidenceKeyValueType.UNIT: _normalise_unit,
    EvidenceKeyValueType.FLOAT: _normalise_float,
    EvidenceKeyValueType.INTEGER: _normalise_integer,
    EvidenceKeyValueType.DATE: _validate_date,
}
""" Per EvidenceKey.value_type, the normaliser normalise_evidence_cell runs on a non-None value. AGE has none """


def normalise_evidence_cell(classification: Classification, cell: VCDataCell, source: str):
    """
    Normalises the value and adds any validation
    """
    # parse all notes looking for references
    value = cell.value
    e_key = cell.e_key
    note = cell.note

    if isinstance(value, str):
        if '\u00a0' in value:
            value = value.replace('\u00a0', ' ')
    if note and '\u00a0' in note:
        note = note.replace('\u00a0', ' ')
        cell.note = note

    if classification.lab.external:
        cell.validate = False

    if value and e_key.value_type in (EvidenceKeyValueType.FREE_ENTRY, EvidenceKeyValueType.TEXT_AREA):
        scan_value = str(value)
        if results := db_ref_regexes.search(scan_value, default_regex=_key_to_regex.get(e_key.key)):
            cell.db_refs = [r.to_json() for r in results]

    if note:
        if results := db_ref_regexes.search(note):
            cell.db_refs = (cell.db_refs or []) + [r.to_json() for r in results]

    if source == SubmissionSource.API and isinstance(value, str):
        # we often get escaped > sent to us as part of
        value = value.replace('&gt;', '>').replace('&lt;', '<')

    if value is None and e_key.mandatory and not e_key.exclude_namespace:
        if e_key.key == SpecialEKeys.CONDITION and classification.condition_resolution_obj:
            # Gene condition is not mandatory if we have a condition resolution
            pass
        else:
            # only provide mandatory validation for internal labs
            # don't want to dirty up imported read only classifications with errors
            cell.add_validation(code=ValidationCode.MANDATORY, severity='error', message='Missing mandatory value')

    # if we got an array of values
    if isinstance(value, list):
        if not value:
            # if the array is empty, treat that as a null value
            value = None
        elif e_key.value_type != EvidenceKeyValueType.MULTISELECT and len(value) == 1:
            # if we're not a multi-select key, and we have an array with 1 value in it.
            # use that 1 value
            value = value[0]

    # types below without a normaliser (age, date, an unmatched user) keep the cleanup above
    if value != cell.value:
        cell.value = value

    if value is not None:
        # remove support for existing_variant_id
        # variant matching is best done through allele_info stuff now
        if e_key.key == 'existing_sample_id':
            try:
                classification.sample = Sample.objects.get(pk=int(value))
                # self.requires_auto_population = True
                return None  # clear out the value
            except Exception:
                cell.add_validation(code=ValidationCode.MATCHING_ERROR, severity='error',
                                    message="Couldn't resolve Sample " + str(value) + ")")

        if normaliser := VALUE_TYPE_NORMALISERS.get(e_key.value_type):
            value = normaliser(cell, value)

        # only complain about unknown keys if we have a value
        if e_key.is_dummy:
            cell.add_validation(code=ValidationCode.UNKNOWN_KEY, severity='warning',
                                message='This key ' + str(e_key.key) + ' has not been registered')

    # SPECIAL HANDLING OF REFSEQ, basically looking for transcript version
    if e_key.key == SpecialEKeys.REFSEQ_TRANSCRIPT_ID:
        if isinstance(value, str):
            version_regex = re.compile('.*?[.][0-9]+')
            if not version_regex.match(value):
                cell.add_validation(code=ValidationCode.UNKNOWN_TRANSCRIPT, severity='warning',
                                    message='Transcript should include version e.g. NM_001256799.1')

    elif e_key.key == SpecialEKeys.AGE:
        if settings.CLASSIFICATION_AUTOFUZZ_AGE:
            cell.value = patch_fuzzy_age(value)

    # ensure we have one non None value before returning the structure
    if not cell.has_data:
        cell.wipe(WipeMode.SET_NONE)


def evidence_validation_messages(evidence: dict) -> list[dict]:
    messages = []

    for key, blob in evidence.items():
        if blob is not None:  # just required to not blow up on legacy, shouldn't have Nones in here
            validations = blob.get('validation')
            if validations is not None:
                if not isinstance(validations, list):
                    validations = [validations]  # just required to not blow up on legacy
                for validation in validations:
                    try:
                        valid_copy = validation.copy()
                        valid_copy['key'] = key
                        messages.append(valid_copy)
                    except Exception:
                        pass
    # make order of messages predictable
    messages.sort(key=lambda m: m.get('key') + ':' + m.get('code'))
    return messages
