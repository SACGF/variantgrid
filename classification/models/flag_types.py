from functools import cached_property

from flags.models.models import FlagType, FlagTypeContext


class ClassificationFlagTypes:

    """
    All flag types for classifications
    """

    @property
    def classification_flag_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk='classification')

    @property
    def submitted_flag(self) -> FlagType:
        return FlagType.objects.get(pk='classification_submitted')

    @property
    def condition_resolution(self) -> FlagType:
        return FlagType.objects.get(pk='condition_resolution')

    @property
    def unshared_flag(self) -> FlagType:
        return FlagType.objects.get(pk='classification_unshared')

    @property
    def suggestion(self) -> FlagType:
        return FlagType.objects.get(pk="classification_suggestion")

    @property
    def significance_change(self) -> FlagType:
        return FlagType.objects.get(pk="classification_significance_change")

    @property
    def classification_outstanding_edits(self) -> FlagType:
        return FlagType.objects.get(pk="classification_outstanding_edits")

    @property
    def classification_clinical_context_changed(self) -> FlagType:
        return FlagType.objects.get(pk="classification_clinical_context_change")

    @property
    def classification_withdrawn(self) -> FlagType:
        return FlagType.objects.get(pk="classification_withdrawn")

    @property
    def internal_review(self) -> FlagType:
        return FlagType.objects.get(pk='classification_internal_review')

    @property
    def discordant(self) -> FlagType:
        return FlagType.objects.get(pk='classification_discordant')

    @property
    def classification_not_public(self) -> FlagType:
        return FlagType.objects.get(pk='classification_not_public')

    @property
    def classification_pending_changes(self) -> FlagType:
        return FlagType.objects.get(pk='classification_pending_changes')

    CLASSIFICATION_PENDING_CHANGES_CLIN_SIG_KEY = "to_clin_sig"
    """
    This is the key that appears in the flag JSON data for a pending change, e.g. {"to_clin_sig": "VUS"}
    """

    # Technically a clinical_context_discordance FlagType

    @cached_property
    def clinical_context_discordance(self) -> FlagType:
        return FlagType.objects.get(pk='clinical_context_discordance')


classification_flag_types = ClassificationFlagTypes()
