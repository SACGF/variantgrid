from lazy import lazy

from flags.models.models import FlagTypeContext, FlagType


class AlleleFlagTypes:

    @lazy
    def allele_flag_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk='allele')

    @lazy
    def missing_37(self) -> FlagType:
        return FlagType.objects.get(pk='allele_missing_37')

    @lazy
    def missing_38(self) -> FlagType:
        return FlagType.objects.get(pk='allele_missing_38')

    @lazy
    def allele_37_not_38(self) -> FlagType:
        return FlagType.objects.get(pk='allele_37_not_38')

allele_flag_types = AlleleFlagTypes()
