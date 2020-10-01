from lazy import lazy
from flags.models.models import FlagType, FlagTypeContext

class VariantClassificationFlagTypes:

    @lazy
    def variant_classification_flag_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk='variant_classification')

    @lazy
    def submitted_flag(self) -> FlagType:
        return FlagType.objects.get(pk='variant_classification_submitted')

    @lazy
    def matching_variant_flag(self) -> FlagType:
        return FlagType.objects.get(pk='variant_classification_matching_variant')

    @lazy
    def matching_variant_warning_flag(self) -> FlagType:
        return FlagType.objects.get(pk='variant_classification_matching_variant_warning')

    @lazy
    def transcript_version_change_flag(self) -> FlagType:
        return FlagType.objects.get(pk='variant_classification_transcript_version_change')

    @lazy
    def unshared_flag(self) -> FlagType:
        return FlagType.objects.get(pk='variant_classification_unshared')

    @lazy
    def suggestion(self) -> FlagType:
        return FlagType.objects.get(pk="variant_classification_suggestion")

    @lazy
    def significance_change(self) -> FlagType:
        return FlagType.objects.get(pk="variant_classification_significance_change")

    @lazy
    def variant_classification_outstanding_edits(self) -> FlagType:
        return FlagType.objects.get(pk="variant_classification_outstanding_edits")

    @lazy
    def variant_classification_clinical_context_changed(self) -> FlagType:
        return FlagType.objects.get(pk="variant_classification_clinical_context_change")

    @lazy
    def variant_classification_withdrawn(self) -> FlagType:
        return FlagType.objects.get(pk="variant_classification_withdrawn")

    # the below are all clinical context flags
    # should probably be put in a different class

    @lazy
    def internal_review(self) -> FlagType:
        return FlagType.objects.get(pk='variant_classification_internal_review')

    @lazy
    def discordant(self) -> FlagType:
        return FlagType.objects.get(pk='variant_classification_discordant')

    @lazy
    def clinical_context_discordance(self) -> FlagType:
        return FlagType.objects.get(pk='clinical_context_discordance')

variant_classification_flag_types = VariantClassificationFlagTypes()
