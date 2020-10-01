from collections import Mapping
from django.contrib.auth.models import User
import copy

from library.utils import is_not_none
from classification.models.evidence_key import EvidenceKey, EvidenceKeyValueType
from classification.models.variant_classification import VariantClassification


class HistoricalEKeyConverter:
    """
    We rename keys or change types, but not all environments will be upgraded (ie SA Path validated releases)
    This converts historical classification keys to the latest ones in Shariant (using latest EKeys)
    """

    HISTORICAL_AND_SHARIANT = [
        # historical is (key_name, actual_historical_type) as the types were wrong...
        (("cadd", None), ("cadd", EvidenceKeyValueType.FLOAT)),
        (("1000_genomes", EvidenceKeyValueType.UNIT), ("1000_genomes_af", EvidenceKeyValueType.UNIT)),
        (("esp", EvidenceKeyValueType.PERCENT), ("esp_af", EvidenceKeyValueType.UNIT)),
        (("exac", EvidenceKeyValueType.PERCENT), ("exac_af", EvidenceKeyValueType.UNIT)),
        (("flossies", EvidenceKeyValueType.PERCENT), ("flossies_af", EvidenceKeyValueType.UNIT)),
        (("gnomad", EvidenceKeyValueType.PERCENT), ("gnomad_af", EvidenceKeyValueType.UNIT)),
        (("gnomad_hp", None), ("gnomad_popmax", EvidenceKeyValueType.SELECT)),
        (("gnomad_hp_maf", None), ("gnomad_popmax_af", EvidenceKeyValueType.UNIT)),
        (("grantham", None), ("grantham", EvidenceKeyValueType.INTEGER)),
        (("uk10k", EvidenceKeyValueType.UNIT), ("uk10k_af", EvidenceKeyValueType.UNIT)),
        (("sample", None), ("sample_type", EvidenceKeyValueType.SELECT)),  # renamed 2019-05-22
        (("variant_type", None), (None, None)),  # No easy mappings - don't send
    ]

    # Everything is stored as string - so arg = string
    TYPE_CONVERTERS = {
        EvidenceKeyValueType.FLOAT: {EvidenceKeyValueType.FREE_ENTRY: str,
                                     EvidenceKeyValueType.INTEGER: lambda s: int(float(s))},
        EvidenceKeyValueType.FREE_ENTRY: {EvidenceKeyValueType.FLOAT: float,
                                          EvidenceKeyValueType.INTEGER: int},
        EvidenceKeyValueType.INTEGER: {EvidenceKeyValueType.FREE_ENTRY: str,
                                       EvidenceKeyValueType.FLOAT: float},
        EvidenceKeyValueType.PERCENT: {EvidenceKeyValueType.UNIT: lambda p: float(p) / 100},
        EvidenceKeyValueType.UNIT: {EvidenceKeyValueType.PERCENT: lambda u: 100 * float(u)},
    }

    def __init__(self):
        historical_keys = filter(is_not_none, [s[0][0] for s in self.HISTORICAL_AND_SHARIANT])
        ekeys = EvidenceKey.objects.filter(key__in=historical_keys)
        historical_types = {ekey.key: ekey.value_type for ekey in ekeys}

        self.legacy_import_user = User.objects.filter(username='legacy_import').first()
        self.historical_to_shariant = {}

        for (historical_key, actual_historical_type), (shariant_key, shariant_type) in self.HISTORICAL_AND_SHARIANT:
            historical_type = historical_types.get(historical_key)
            if historical_type:  # Make sure it's in DB (ie not using latest keys)
                if actual_historical_type:
                    historical_type = actual_historical_type

                #print(f"({historical_key}, {actual_historical_type}), ({shariant_key}, {shariant_type})")
                if historical_key != shariant_key or historical_type != shariant_type:
                    to_shariant = None
                    if historical_key and shariant_key and historical_type != shariant_type:
                        try:
                            to_shariant = self.TYPE_CONVERTERS[historical_type][shariant_type]
                        except:
                            msg = f"{historical_key} -> {shariant_key}. Don't know how to convert from {historical_type} => {shariant_type}"
                            raise NotImplementedError(msg)

                    self.historical_to_shariant[historical_key] = (shariant_key, to_shariant)

    @staticmethod
    def _convert_keys(data, key_mappings, add_note=False):
        data = copy.deepcopy(data)
        converted_data = {}
        for from_key in data.keys():
            to_key_and_converter = key_mappings.get(from_key)
            if to_key_and_converter:
                (to_key, converter) = to_key_and_converter
                if to_key:
                    from_value = VariantClassification.get_optional_value_from(data, from_key)
                    if from_value is None:
                        # no need to convert Nones, just move the key over
                        converted_data[to_key] = data.get(from_key)
                    else:
                        if converter:
                            to_value = converter(from_value)
                        else:
                            to_value = from_value
                        to_value = str(to_value)

                        valueObj = data[from_key]
                        if not isinstance(valueObj, Mapping):
                            valueObj = {}

                        valueObj["value"] = to_value
                        if add_note and from_value != to_value:
                            notes = []
                            existing_note = valueObj.get("note")
                            if existing_note:
                                notes.append(existing_note)
                            notes.append(f"Converted from: '{from_key}'='{from_value}'")
                            valueObj["note"] = ". ".join(notes)

                        converted_data[to_key] = valueObj
            else:
                converted_data[from_key] = data[from_key]  # Straight copy

        return converted_data

    def to_shariant(self, vcm, data, add_note=True):
        if self.legacy_import_user and vcm.user == self.legacy_import_user:
            shariant_data = data  # God knows how this data was formatted....
        else:
            shariant_data = self._convert_keys(data, self.historical_to_shariant, add_note=add_note)

        return shariant_data

    # There's also a from_shariant() in source history if we ever want it again...
