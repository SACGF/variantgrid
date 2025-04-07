import itertools
from typing import Optional
from classification.models import ConditionResolved, ClinVarExport, ClinVarAllele, ClassificationModification
from dataclasses import dataclass
from classification.models.clinvar_export_prepare import ClassificationModificationCandidate


@dataclass
class ClinVarExportStub:
    new_condition_umbrella: ConditionResolved
    new_classification_modification: ClassificationModification
    clinvar_export: Optional[ClinVarExport] = None

    @property
    def condition_umbrell(self) -> ConditionResolved:
        if condition_umbrella := self.new_condition_umbrella:
            return condition_umbrella
        if clinvar_export := self.clinvar_export:
            return clinvar_export.condition_resolved
        raise ValueError("ClinVarExportStub without condition")

    def assign_if_newer(self, classification_modification: ClassificationModification):
        if existing := self.new_classification_modification:
            if existing.curated_date > classification_modification:
                return
        self.new_classification_modification = classification_modification


    @staticmethod
    def from_clinvar_export(clinvar_export: ClinVarExport) -> 'ClinVarExportStub':
        return ClinVarExportStub(
            condition_umbrella=None,
            classification=None,
            clinvar_export=clinvar_export
        )


class ClinVarExportManager:

    def __init__(self, classification_modifications: list[ClassificationModification], existing_clinvar_exports: list[ClinVarExport]):
        self.stubs = [ClinVarExportStub.from_clinvar_export(clinvar_export) for clinvar_export in existing_clinvar_exports]

        # group classification modifications by condition, so and get the most up to date for each condition
        sorted_by_condition = classification_modifications.sort(key=lambda c: c.classification.condition_resolution)
        condition: ConditionResolved
        for condition, grouped_by_condition in itertools.groupby(sorted_by_condition, lambda c: c.classification.condition_resolution):
            latest_for_condition = max(grouped_by_condition, key=lambda c: c.curated_date_check)

            # do a quick check to see if the condition matches an umbrella exactly
            use_stub = None
            for stub in self.stubs:
                if stub.condition_umbrella == condition:
                    use_stub = stub

            if not use_stub:
                # now see which umbrella if any this falls under
                min_score = None
                for stub in self.stubs:
                    stub_score = condition.same_or_more_specific_step_count(stub.condition_umbrella)
                    if stub_score is not None and stub_score < min_score:
                        min_score = stub_score
                        use_stub = stub

            if use_stub:
                use_stub.assign_if_newer(latest_for_condition)
            else:
                # now see if we extend a stub
                min_score = None
                for stub in self.stubs:
                    stub_score = stub.condition_umbrella.same_or_more_specific_step_count(condition)
                    if stub_score is not None and stub_score < min_score:
                        min_score = stub_score
                        use_stub = stub
                if use_stub:
                    use_stub.new_condition_umbrella = condition
                    use_stub.assign_if_newer(latest_for_condition)
                else:
                    # wont fit into existing stub or even when extending it, have to create a new one
                    new_stub = ClinVarExportStub(new_condition_umbrella=condition, new_classification_modification=latest_for_condition)
                    self.stubs.append(new_stub)
