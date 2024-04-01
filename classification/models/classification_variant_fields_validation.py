import re
from re import RegexFlag
from typing import Optional

from django.conf import settings
from django.dispatch.dispatcher import receiver

from classification.enums import ValidationCode, SpecialEKeys, AlleleOriginBucket
from classification.models import EvidenceKeyMap, PatchMeta
from classification.models.classification import Classification, \
    classification_validation_signal
from classification.models.classification_utils import ValidationMerger
from genes.hgvs import HGVSMatcher
from genes.models import NoTranscript
from snpdb.models import GenomeBuild, Variant, VariantCoordinate


__AT_LEAST_ONE_SET = {SpecialEKeys.CLINICAL_SIGNIFICANCE, SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE}
__CHECK_IF_CHANGED = __AT_LEAST_ONE_SET | {SpecialEKeys.ALLELE_ORIGIN}

@receiver(classification_validation_signal, sender=Classification)
def validate_variant_classification_significance(sender, patch_meta: PatchMeta, key_map: EvidenceKeyMap, **kwargs) -> Optional[ValidationMerger]:
    """
    Validates that at l
    east one of clinical significance or somatic clinical significance has a value
    """
    if patch_meta.intersection_modified(__CHECK_IF_CHANGED):
        vm = ValidationMerger()
        vm.tested(
            keys=[SpecialEKeys.CLINICAL_SIGNIFICANCE, SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE],
            codes=[
                ValidationCode.MISSING_CLINICAL_SIGNIFICANCE,
                ValidationCode.INVALID_FIELD_FOR_SOMATIC,
                ValidationCode.INVALID_FIELD_FOR_GERMLINE
            ]
        )

        allele_origin = patch_meta.get(SpecialEKeys.ALLELE_ORIGIN, fallback_existing=True)
        bucket = AlleleOriginBucket.bucket_for_allele_origin(allele_origin)

        if bucket == AlleleOriginBucket.GERMLINE:
            if not patch_meta.get(SpecialEKeys.CLINICAL_SIGNIFICANCE):
                vm.add_message(
                    SpecialEKeys.CLINICAL_SIGNIFICANCE,
                    code=ValidationCode.MISSING_CLINICAL_SIGNIFICANCE,
                    severity='error',
                    message=f"{key_map.get(SpecialEKeys.CLINICAL_SIGNIFICANCE).label} requires a value"
                )
            if patch_meta.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE):

                allele_origin_pretty_value = key_map.get(SpecialEKeys.ALLELE_ORIGIN).pretty_value(allele_origin)
                if not allele_origin_pretty_value:
                    allele_origin_pretty_value = "'blank'"

                vm.add_message(
                    SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE,
                    code=ValidationCode.INVALID_FIELD_FOR_GERMLINE,
                    severity='warning',
                    message=f"This field is not valid for the given {key_map.get(SpecialEKeys.ALLELE_ORIGIN).pretty_label} {allele_origin_pretty_value}"
                )
            return vm

        has_value = False
        for code in __AT_LEAST_ONE_SET:
            if patch_meta.get(code, fallback_existing=True):
                has_value = True
                break

        if not has_value:
            # TODO this test is pretty basic compared to what
            likely_somatic = "somatic" in (patch_meta.get(SpecialEKeys.ALLELE_ORIGIN, fallback_existing=True) or "").lower()
            message: str
            if bucket == AlleleOriginBucket.GERMLINE:
                message = f"{key_map.get(SpecialEKeys.CLINICAL_SIGNIFICANCE).label} or {key_map.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).label} requires a value"
            else:
                message = f"{key_map.get(SpecialEKeys.CLINICAL_SIGNIFICANCE).label} requires a value"

            vm.add_message(
                SpecialEKeys.CLINICAL_SIGNIFICANCE,
                code=ValidationCode.MISSING_CLINICAL_SIGNIFICANCE,
                severity='error',
                message=message
            )
        return vm


__LEVELS_AND_TIER = set(list(SpecialEKeys.AMP_LEVELS_TO_LEVEL.keys()) + [SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE])
__LEVELS_TO_TIER = {
    "A": "1",
    "B": "1",
    "C": "2",
    "D": "2"
}
__TIER_TO_ROMAN = {
    "1": "I",
    "2": "II",
    "3": "III",
    "4": "IV"
}


@receiver(classification_validation_signal, sender=Classification)
def validate_letter_to_tier(sender, patch_meta: PatchMeta, key_map: EvidenceKeyMap, **kwargs) -> Optional[ValidationMerger]:
    """
    Validates that at l
    east one of clinical significance or somatic clinical significance has a value
    """
    max_level: Optional[str] = None
    if patch_meta.intersection_modified(__LEVELS_AND_TIER):
        vm = ValidationMerger()
        vm.tested(
            keys=[SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE],
            codes=[ValidationCode.SOMATIC_MISMATCHED_LEVEL]
        )

        for level_key, level in SpecialEKeys.AMP_LEVELS_TO_LEVEL.items():
            if patch_meta.get(level_key, fallback_existing=True):
                max_level = level
                break

        e_key_somatic_clin_sig = key_map.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
        max_level_tier = __LEVELS_TO_TIER.get(max_level)
        somatic_clin_sig = patch_meta.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE, fallback_existing=True)
        if max_level_tier and somatic_clin_sig:
            if actual_tier_level := e_key_somatic_clin_sig.option_dictionary_property("tier").get(somatic_clin_sig, "X"):
                if actual_tier_level != max_level_tier:
                    somatic_clin_sig_pretty_value = e_key_somatic_clin_sig.pretty_value(somatic_clin_sig)
                    vm.add_message(
                        SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE,
                        code=ValidationCode.SOMATIC_MISMATCHED_LEVEL,
                        severity='warning',
                        message=f"Level {max_level} indicates {e_key_somatic_clin_sig.pretty_label} should be \"Tier {__TIER_TO_ROMAN.get(max_level_tier, max_level_tier)}\" not \"{somatic_clin_sig_pretty_value}\""
                    )
        return vm


@receiver(classification_validation_signal, sender=Classification)
def validate_variant_fields(sender, patch_meta: PatchMeta, key_map: EvidenceKeyMap, **kwargs) -> Optional[ValidationMerger]:  # pylint: disable=unused-argument
    """
    Ensures hgvs fields are valid, and that they match each other and the variant coordinate
    """
    vm = ValidationMerger()

    if not settings.CLASSIFICATION_MATCH_VARIANTS:
        # while this technically isn't about matching a variant, much of the work in get_variant_tuple -> pyhgvs
        # requires
        return vm

    VARIANT_LINKING_KEYS_SET = set(SpecialEKeys.VARIANT_LINKING_KEYS)
    # Only do this validation if at least variant linking key value has changed
    variant_values = {}
    variant_map = {}

    build_name = patch_meta.get(SpecialEKeys.GENOME_BUILD)
    if not build_name:
        return None

    build_changed = patch_meta.is_modified(SpecialEKeys.GENOME_BUILD)
    variant_key_changed = patch_meta.intersection_modified(VARIANT_LINKING_KEYS_SET)
    if build_changed or variant_key_changed:
        for evidence_key in SpecialEKeys.VARIANT_LINKING_KEYS:
            value = patch_meta.get(evidence_key)
            if value:
                # only consider it if it doesn't have other outstanding validation messages
                variant_values[evidence_key] = value

    if len(variant_values) >= 1:
        vm.tested(
            keys=VARIANT_LINKING_KEYS_SET,
            codes=[ValidationCode.MATCHING_ERROR, ValidationCode.INCONSISTENT_VARIANT, ValidationCode.MATCHING_ERROR]
        )

        try:
            # Take values from patch (not saved yet)
            genome_build = GenomeBuild.get_name_or_alias(build_name)
        except GenomeBuild.DoesNotExist:
            genome_build = None

        if genome_build:
            hgvs_matcher = HGVSMatcher(genome_build)
            for evidence_key, variant_value in variant_values.items():
                try:
                    if evidence_key in SpecialEKeys.VARIANT_LINKING_HGVS_KEYS:

                        # Not unsupported transcript type is handled by Imported Allele validation (not classification)
                        # if evidence_key == SpecialEKeys.C_HGVS and not Classification.is_supported_transcript(variant_value):
                        #     vm.add_message(evidence_key, code=ValidationCode.MATCHING_ERROR, severity='error', message='Transcript type not yet supported')
                        #     continue

                        variant_coordinate = hgvs_matcher.get_variant_coordinate(variant_value)
                    elif evidence_key == SpecialEKeys.VARIANT_COORDINATE:
                        variant_coordinate = VariantCoordinate.from_string(variant_value)
                    else:
                        raise ValueError(f'Unexpected evidence_key for variant matching {evidence_key}')

                    variant_str = Variant.format_tuple(*variant_coordinate)
                    variant_map[variant_str] = variant_map.get(variant_str, []) + [evidence_key]
                except NoTranscript:
                    # this should find its way into Validation Matching flag
                    pass
                except:
                    # this should be handled by ImportedAlleleInfo
                    pass
                # except ValueError as ve:
                #     # shouldn't happen as other validation should have already done this
                #     vm.add_message(evidence_key, code=ValidationCode.MATCHING_ERROR, severity='error', message=f'Error attempting to parse variant coordinate ({ve})')
                # except:
                #     vm.add_message(evidence_key, code=ValidationCode.MATCHING_ERROR, severity='error', message='Error attempting to parse variant coordinate')

    # check for this because if there were parsing errors we wont have any variant_values
    if len(variant_map) > 1:
        # generate message like:
        # c HGVS resolves to 4:53454A>T but g HGVS resolves to 4:4:53454A>C
        resolutions = ['%s resolves to (%s)' % (' and '.join([key_map.get(key).pretty_label for key in keys]), variant_str) for variant_str, keys in variant_map.items()]
        message = ' but '.join(resolutions)

        vm.add_message(SpecialEKeys.VARIANT_COORDINATE, code=ValidationCode.INCONSISTENT_VARIANT, severity='warning', message=message)

    return vm


SUB_CLIN_SIG = re.compile(r'this variant .{2,20}? classified .*?[.]', RegexFlag.IGNORECASE | RegexFlag.DOTALL)


@receiver(classification_validation_signal, sender=Classification)
def validate_clinical_significance(sender, patch_meta: PatchMeta, **kwargs) -> Optional[ValidationMerger]:  # pylint: disable=unused-argument
    """
    Warns if the text description seems to imply one clinical significance wheras the actua value states another
    """

    vm = ValidationMerger()
    vm.tested(
        keys=[SpecialEKeys.CLINICAL_SIGNIFICANCE],
        codes=['cs_mismatch']
    )

    if patch_meta.is_modified(SpecialEKeys.CLINICAL_SIGNIFICANCE) or patch_meta.is_modified(SpecialEKeys.INTERPRETATION_SUMMARY):
        cs: str = patch_meta.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        inter_summary: str = patch_meta.get(SpecialEKeys.INTERPRETATION_SUMMARY)

        if cs and inter_summary:
            attempt_clin_sig = SUB_CLIN_SIG.search(inter_summary)
            if attempt_clin_sig:
                clin_text = attempt_clin_sig.group(0).lower()

                actual_cs_parts: set[str] = set()

                if cs.startswith('VUS'):  # we don't complain about VUS A,B,C differences
                    actual_cs_parts.add('VUS')
                else:
                    if cs.startswith('L'):
                        actual_cs_parts.add('likely')
                    if cs.endswith('B'):
                        actual_cs_parts.add('benign')
                    elif cs.endswith('P'):
                        actual_cs_parts.add('pathogenic')
                    elif cs == 'R':
                        actual_cs_parts.add('risk factory')
                    elif cs == 'D':
                        actual_cs_parts.add('drug response')
                    elif cs == 'A':
                        actual_cs_parts.add('artefact')

                detected_cs_parts: set[str] = set()
                if 'vous' in clin_text or 'vus' in clin_text or 'variant of uncertain' in clin_text:
                    detected_cs_parts.add('VUS')
                if 'like' in clin_text:
                    detected_cs_parts.add('likely')
                if 'benign' in clin_text:
                    detected_cs_parts.add('benign')
                if 'pathogenic' in clin_text:
                    detected_cs_parts.add('pathogenic')
                if 'risk factor' in clin_text:
                    detected_cs_parts.add('risk factor')
                if 'drug response' in clin_text:
                    detected_cs_parts.add('drug response')
                if 'artefact' in clin_text:
                    detected_cs_parts.add('artefact')

                if detected_cs_parts and actual_cs_parts != detected_cs_parts:
                    detected_text = ' '.join(detected_cs_parts)
                    vm.add_message(
                        key=SpecialEKeys.CLINICAL_SIGNIFICANCE,
                        code='cs_mismatch',
                        severity='warning',
                        message=f'Classified as {EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(cs)} but detected terms ({detected_text}) in {EvidenceKeyMap.cached_key(SpecialEKeys.INTERPRETATION_SUMMARY).pretty_label}'
                    )
    return vm
