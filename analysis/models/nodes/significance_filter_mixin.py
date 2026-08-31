from analysis.models.enums import NodeMatchInput
from snpdb.models.models_enums import AlleleOriginFilterDefault


class SignificanceFilterNodeMixin:
    """ Nodes that filter on clinical significance pills, split by allele origin, over either a parent's
        variants or everything matching. Expects `node_input` and `allele_origin` fields. """

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        if self.node_input == NodeMatchInput.MATCHING_VARIANTS:
            return 0
        return 1

    @property
    def allele_origin_filter(self) -> AlleleOriginFilterDefault:
        return AlleleOriginFilterDefault(self.allele_origin)

    @property
    def germline_enabled(self) -> bool:
        return self.allele_origin_filter != AlleleOriginFilterDefault.SOMATIC

    @property
    def somatic_enabled(self) -> bool:
        return self.allele_origin_filter != AlleleOriginFilterDefault.GERMLINE

    def _selected_values(self, field_values: dict[str, str]) -> list[str]:
        return [value for field, value in field_values.items() if getattr(self, field)]
