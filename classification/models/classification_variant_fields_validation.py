import re
from re import RegexFlag
from typing import Optional, Set

from django.conf import settings
from django.dispatch.dispatcher import receiver

from classification.enums import ValidationCode, SpecialEKeys
from classification.models import EvidenceKeyMap, PatchMeta
from classification.models.classification import Classification, \
    classification_validation_signal
from classification.models.classification_utils import ValidationMerger
from genes.hgvs import HGVSMatcher
from genes.models import NoTranscript
from snpdb.models import GenomeBuild, Variant

VARIANT_VALIDATING_CODES = {ValidationCode.MATCHING_ERROR, ValidationCode.INCONSISTENT_VARIANT, ValidationCode.MATCHING_ERROR}


@receiver(classification_validation_signal, sender=Classification)
def validate_variant_fields(sender, **kwargs) -> Optional[ValidationMerger]:  # pylint: disable=unused-argument
    """
    Ensures hgvs fields are valid, and that they match each other and the variant coordinate
    """
    vm = ValidationMerger()

    if not settings.VARIANT_CLASSIFICATION_MATCH_VARIANTS:
        # while this technically isn't about matching a variant, much of the work in get_variant_tuple -> pyhgvs
        # requires
        return vm

    patch_meta: PatchMeta = kwargs.get('patch_meta')
    evidence_key_map: EvidenceKeyMap = kwargs.get('key_map')

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

                        variant_tuple = hgvs_matcher.get_variant_tuple(variant_value)
                    elif evidence_key == SpecialEKeys.VARIANT_COORDINATE:
                        variant_tuple = Variant.get_tuple_from_string(variant_value, genome_build)
                    else:
                        raise ValueError(f'Unexpected evidence_key for variant matching {evidence_key}')

                    variant_str = Variant.format_tuple(*variant_tuple)
                    variant_map[variant_str] = variant_map.get(variant_str, []) + [evidence_key]
                except NoTranscript:
                    # this should find its way into Validation Matching flag
                    pass
                except ValueError as ve:
                    # shouldn't happen as other validation should have already done this
                    vm.add_message(evidence_key, code=ValidationCode.MATCHING_ERROR, severity='error', message=f'Error attempting to parse variant coordinate ({ve})')
                except:
                    vm.add_message(evidence_key, code=ValidationCode.MATCHING_ERROR, severity='error', message='Error attempting to parse variant coordinate')

    # check for this because if there were parsing errors we wont have any variant_values
    if len(variant_map) > 1:
        # generate message like:
        # c HGVS resolves to 4:53454A>T but g HGVS resolves to 4:4:53454A>C
        resolutions = ['%s resolves to (%s)' % (' and '.join([evidence_key_map.get(key).pretty_label for key in keys]), variant_str) for variant_str, keys in variant_map.items()]
        message = ' but '.join(resolutions)

        vm.add_message(SpecialEKeys.VARIANT_COORDINATE, code=ValidationCode.INCONSISTENT_VARIANT, severity='warning', message=message)

    return vm


# SUB_VUS_1 = re.compile(r'this variant .{2,20}? classified .{0,5}?3(A|B|C).{0,10}?V(?:O)?US', RegexFlag.IGNORECASE)
# SUB_VUS_2 = re.compile(r'this variant .{2,20}? classified .{0,5}?V(?:O)?US.{0,10}?3(A|B|C)', RegexFlag.IGNORECASE)
SUB_CLIN_SIG = re.compile(r'this variant .{2,20}? classified .*?[.]', RegexFlag.IGNORECASE | RegexFlag.DOTALL)


@receiver(classification_validation_signal, sender=Classification)
def validate_clinical_significance(sender, **kwargs) -> Optional[ValidationMerger]:  # pylint: disable=unused-argument
    vm = ValidationMerger()
    vm.tested(
        keys=[SpecialEKeys.CLINICAL_SIGNIFICANCE],
        codes=['cs_mismatch']
    )

    patch_meta: PatchMeta = kwargs.get('patch_meta')

    if patch_meta.is_modified(SpecialEKeys.CLINICAL_SIGNIFICANCE) or patch_meta.is_modified(SpecialEKeys.INTERPRETATION_SUMMARY):
        cs: str = patch_meta.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        inter_summary: str = patch_meta.get(SpecialEKeys.INTERPRETATION_SUMMARY)

        if cs and inter_summary:
            attempt_clin_sig = SUB_CLIN_SIG.search(inter_summary)
            if attempt_clin_sig:
                clin_text = attempt_clin_sig.group(0).lower()

                actual_cs_parts: Set[str] = set()

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

                detected_cs_parts: Set[str] = set()
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
