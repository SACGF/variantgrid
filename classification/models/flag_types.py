from lazy import lazy

from flags.models.models import FlagType, FlagTypeContext


class ClassificationFlagTypes:

    @lazy
    def classification_flag_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk='classification')

    @lazy
    def submitted_flag(self) -> FlagType:
        return FlagType.objects.get(pk='classification_submitted')

    @lazy
    def matching_variant_flag(self) -> FlagType:
        return FlagType.objects.get(pk='classification_matching_variant')

    @lazy
    def matching_variant_warning_flag(self) -> FlagType:
        return FlagType.objects.get(pk='classification_matching_variant_warning')

    @lazy
    def transcript_version_change_flag(self) -> FlagType:
        return FlagType.objects.get(pk='classification_transcript_version_change')

    @lazy
    def unshared_flag(self) -> FlagType:
        return FlagType.objects.get(pk='classification_unshared')

    @lazy
    def suggestion(self) -> FlagType:
        return FlagType.objects.get(pk="classification_suggestion")

    @lazy
    def significance_change(self) -> FlagType:
        return FlagType.objects.get(pk="classification_significance_change")

    @lazy
    def classification_outstanding_edits(self) -> FlagType:
        return FlagType.objects.get(pk="classification_outstanding_edits")

    @lazy
    def classification_clinical_context_changed(self) -> FlagType:
        return FlagType.objects.get(pk="classification_clinical_context_change")

    @lazy
    def classification_withdrawn(self) -> FlagType:
        return FlagType.objects.get(pk="classification_withdrawn")

    # the below are all clinical context flags
    # should probably be put in a different class

    @lazy
    def internal_review(self) -> FlagType:
        return FlagType.objects.get(pk='classification_internal_review')

    @lazy
    def discordant(self) -> FlagType:
        return FlagType.objects.get(pk='classification_discordant')

    @lazy
    def clinical_context_discordance(self) -> FlagType:
        return FlagType.objects.get(pk='clinical_context_discordance')

    @lazy
    def classification_not_public(self) -> FlagType:
        return FlagType.objects.get(pk='classification_not_public')

    @lazy
    def classification_pending_changes(self) -> FlagType:
        return FlagType.objects.get(pk='classification_pending_changes')

    CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY = "to_clin_sig"


classification_flag_types = ClassificationFlagTypes()
