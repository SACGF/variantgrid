from functools import cached_property

from flags.models.models import FlagTypeContext, FlagType


class AlleleFlagTypes:

    @cached_property
    def allele_flag_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk='allele')

    # @cached_property
    # def missing_37(self) -> FlagType:
    #     return FlagType.objects.get(pk='allele_missing_37')

    # @cached_property
    # def missing_38(self) -> FlagType:
    #     return FlagType.objects.get(pk='allele_missing_38')

    # @cached_property
    # def allele_37_not_38(self) -> FlagType:
    #     return FlagType.objects.get(pk='allele_37_not_38')

allele_flag_types = AlleleFlagTypes()
