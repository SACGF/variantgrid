from functools import cached_property

from flags.models.models import FlagType, FlagTypeContext


class ClassificationFlagTypes:

    @cached_property
    def classification_flag_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk='classification')

    @cached_property
    def submitted_flag(self) -> FlagType:
        return FlagType.objects.get(pk='classification_submitted')

    @cached_property
    def unshared_flag(self) -> FlagType:
        return FlagType.objects.get(pk='classification_unshared')

    @cached_property
    def suggestion(self) -> FlagType:
        return FlagType.objects.get(pk="classification_suggestion")

    @cached_property
    def significance_change(self) -> FlagType:
        return FlagType.objects.get(pk="classification_significance_change")

    @cached_property
    def classification_outstanding_edits(self) -> FlagType:
        return FlagType.objects.get(pk="classification_outstanding_edits")

    @cached_property
    def classification_clinical_context_changed(self) -> FlagType:
        return FlagType.objects.get(pk="classification_clinical_context_change")

    @cached_property
    def classification_withdrawn(self) -> FlagType:
        return FlagType.objects.get(pk="classification_withdrawn")

    # the below are all clinical context flags
    # should probably be put in a different class

    @cached_property
    def internal_review(self) -> FlagType:
        return FlagType.objects.get(pk='classification_internal_review')

    @cached_property
    def discordant(self) -> FlagType:
        return FlagType.objects.get(pk='classification_discordant')

    @cached_property
    def clinical_context_discordance(self) -> FlagType:
        return FlagType.objects.get(pk='clinical_context_discordance')

    @cached_property
    def classification_not_public(self) -> FlagType:
        return FlagType.objects.get(pk='classification_not_public')

    @cached_property
    def classification_pending_changes(self) -> FlagType:
        return FlagType.objects.get(pk='classification_pending_changes')

    CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY = "to_clin_sig"


classification_flag_types = ClassificationFlagTypes()
