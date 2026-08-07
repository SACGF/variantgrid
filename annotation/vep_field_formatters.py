""" Per-column value cleaning functions for VEP annotation fields.

    These are pure ``value -> value`` transforms attached to individual
    ``VEPColumnDef`` rows (via their ``formatter`` field) in
    ``annotation.vep_columns``. ``BulkVEPVCFAnnotationInserter`` applies them
    generically after copying raw VEP output into the destination columns.

    Cross-field / derived processing (ALoFT transcript picking, dbNSFP
    per-transcript picking, gnomAD SV overlap, calculated columns) is NOT here -
    it lives in the inserter because it cannot be expressed as a single-column
    transform.
"""
import operator

from annotation.models.damage_enums import SIFTPrediction, AlphaMissensePrediction
from library.utils import invert_dict

VEP_SEPARATOR = '&'
EMPTY_VALUES = {'', '.'}


def gnomad_filtered_func(raw_value):
    """ We use FILTER in Gnomad3 (GRCh38 only) - need to convert back to bool
        In the combined exomes/genomes (gnomad2, gnomad4) we use gnomad_filtered=1
        So don't need to format this etc
    """
    return raw_value not in (None, "PASS")


def format_hgnc_id(raw_value):
    """ VEP GRCh37 returns 55 while GRCh38 returns "HGNC:55" """
    return int(raw_value.replace("HGNC:", ""))


def format_vep_sift_to_choice(vep_sift):
    """ we ignore the low_confidence calls to make it simpler """
    if vep_sift.startswith("deleterious"):
        return SIFTPrediction.DAMAGING
    elif vep_sift.startswith("tolerated"):
        return SIFTPrediction.TOLERATED
    raise ValueError(f"Unknown SIFT value: '{vep_sift}'")


def get_format_alphamissense_class_func():
    """ GRCh37 has 'benign' while GRCh38 has 'likely_benign'
        @see https://github.com/Ensembl/VEP_plugins/issues/668
    """
    cff = get_choice_formatter_func(AlphaMissensePrediction.choices)

    def _format_alphamissense_class(alphamissense_class):
        if alphamissense_class == "benign":
            alphamissense_class = "likely_benign"
        return cff(alphamissense_class)
    return _format_alphamissense_class


def get_extract_existing_variation(prefix):
    def format_vep_existing_variation(vep_existing_variation):
        ev_list = vep_existing_variation.split(VEP_SEPARATOR)
        cosmic_ids = [ev for ev in ev_list if ev.startswith(prefix)]
        return VEP_SEPARATOR.join(sorted(set(cosmic_ids)))

    return format_vep_existing_variation


def format_vep_somatic(raw_value):
    return "1" in raw_value


def get_choice_formatter_func(choices, empty_values=None):
    lookup = invert_dict(dict(choices))

    def format_choice(raw_value):
        if empty_values is not None:
            if raw_value in empty_values:
                return None
        return lookup[raw_value]

    return format_choice


def get_clean_and_pick_single_value_func(pick_single_value_func, cast_func=None, empty_values=None):
    """ Returns a function to clean and pick single value.
        casting is performed before calling pick_single_value_func so you can call min/max """

    if empty_values is None:
        empty_values = EMPTY_VALUES

    def _clean_and_pick_single_value_func(raw_value):
        it = (tm for tm in raw_value.split(VEP_SEPARATOR) if tm != '')
        # Handle '.'
        if cast_func:
            values = [cast_func(v) for v in it if v not in empty_values]
        else:
            values = [v for v in it if v not in empty_values]
        value = None
        if values:
            value = pick_single_value_func(values)
        return value

    return _clean_and_pick_single_value_func


def join_uniq(multiple_values):
    return VEP_SEPARATOR.join(set(multiple_values))


def get_most_damaging_func(klass):
    def get_most_damaging(multiple_values):
        prediction_list = multiple_values.split(VEP_SEPARATOR)
        return klass.get_most_damaging(prediction_list)

    return get_most_damaging


def get_format_empty_as_none(empty_values: set):
    def format_empty_as_none(val):
        if val in empty_values:
            val = None
        return val
    return format_empty_as_none


def format_nmd_escaping_variant(value) -> bool:
    return value == "NMD_escaping_variant"


def format_aloft_high_confidence(value) -> bool:
    high_confidence = None
    if value == "High":
        high_confidence = True
    elif value == "Low":
        high_confidence = False
    return high_confidence


def format_canonical(value) -> bool:
    return value == "YES"


# dbNSFP 5.3.1a (columns_version >= 4) emits AlphaMissense_pred as 'B'/'LB'/'A'/'LP'/'P',
# optionally as &-separated arrays parallel to Ensembl_transcriptid. Map to the
# AlphaMissensePrediction code and pick the most-pathogenic value across the array.
_ALPHAMISSENSE_PRED_MAP = {
    "B": AlphaMissensePrediction.BENIGN,
    "LB": AlphaMissensePrediction.LIKELY_BENIGN,
    "A": AlphaMissensePrediction.AMBIGUOUS,
    "LP": AlphaMissensePrediction.LIKELY_PATHOGENIC,
    "P": AlphaMissensePrediction.PATHOGENIC,
}
_ALPHAMISSENSE_PRED_ORDER = [
    AlphaMissensePrediction.PATHOGENIC,
    AlphaMissensePrediction.LIKELY_PATHOGENIC,
    AlphaMissensePrediction.AMBIGUOUS,
    AlphaMissensePrediction.LIKELY_BENIGN,
    AlphaMissensePrediction.BENIGN,
]


def format_alphamissense_pred(raw_value):
    seen = {_ALPHAMISSENSE_PRED_MAP[v] for v in raw_value.split(VEP_SEPARATOR) if v in _ALPHAMISSENSE_PRED_MAP}
    for code in _ALPHAMISSENSE_PRED_ORDER:
        if code in seen:
            return code
    return None


# --- Shared factory instances -------------------------------------------------
# Pre-built once so VEPColumnDefs can reference them directly. These mirror the
# locals previously built inside BulkVEPVCFAnnotationInserter._add_vep_field_handlers.

# TOPMED and 1k genomes can return multiple values - take highest
format_pick_highest_float = get_clean_and_pick_single_value_func(max, float)
format_pick_highest_int = get_clean_and_pick_single_value_func(max, int)
format_sum_int = get_clean_and_pick_single_value_func(sum, int)
remove_empty_multiples = get_clean_and_pick_single_value_func(join_uniq)

# MaveDB emits "NA" for missing - treat "NA"/empty as null and take the lowest score
_EMPTY_MAVE_FLOAT_VALUES = EMPTY_VALUES | {"NA"}
format_pick_lowest_float = get_clean_and_pick_single_value_func(min, float,
                                                                empty_values=_EMPTY_MAVE_FLOAT_VALUES)
# OpenTargets emits "NA" for missing and &-joins multiple GWAS entries (e.g. "NA&NA") - treat
# "NA"/empty as null and take the highest (most confident) locus-to-gene score.
format_pick_highest_float_na = get_clean_and_pick_single_value_func(max, float,
                                                                    empty_values=_EMPTY_MAVE_FLOAT_VALUES)

format_empty_as_none = get_format_empty_as_none(empty_values=EMPTY_VALUES)

# COSMIC v90 (5/9/2019) switched to COSV (build independent identifiers)
extract_cosmic = get_extract_existing_variation("COSV")
extract_dbsnp = get_extract_existing_variation("rs")

# Mastermind_counts is "&"-joined MMCNT1&MMCNT2&MMCNT3 - pick the positional value
format_mastermind_count_1 = get_clean_and_pick_single_value_func(operator.itemgetter(0), int)
format_mastermind_count_2 = get_clean_and_pick_single_value_func(operator.itemgetter(1), int)
format_mastermind_count_3 = get_clean_and_pick_single_value_func(operator.itemgetter(2), int)
