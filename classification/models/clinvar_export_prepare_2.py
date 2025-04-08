import itertools
from typing import Optional
from django.utils.timezone import now
from django.db.models import QuerySet

from classification.enums import ShareLevel, AlleleOriginBucket
from classification.models import ConditionResolved, ClinVarExport, ClinVarAllele, ClassificationModification, \
    ClinVarExportTypeBucket
from dataclasses import dataclass
from snpdb.models import Allele, ClinVarKey


@dataclass
class ClinVarExportStub:
    clinvar_allele: ClinVarAllele
    new_condition_umbrella: ConditionResolved
    new_classification_modification: ClassificationModification
    clinvar_export: Optional[ClinVarExport] = None

    @property
    def condition_umbrella(self) -> ConditionResolved:
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

    def apply(self) -> list[str]:
        log = []
        if clinvar_export := self.clinvar_export:
            if not self.new_classification_modification:
                if not clinvar_export.has_submission:
                    log.append(f"Deleting CE_{clinvar_export.pk} as no submissions and no classification attached")
                    clinvar_export.delete()
                    return log

            has_changes = False
            if new_condition_umbrella := self.new_condition_umbrella:
                if new_condition_umbrella != clinvar_export.condition:
                    log.append(f"Updating CE_{clinvar_export.pk} with new condition {clinvar_export.condition.summary} -> {new_condition_umbrella.summary}")
                    clinvar_export.condition = new_condition_umbrella
                    has_changes = True
            if new_classification_modification := self.new_classification_modification:
                if new_classification_modification != clinvar_export.classification_based_on:
                    log.append(f"Updating CE_{clinvar_export.pk} with new classification modification {clinvar_export.classification_based_on} -> {new_classification_modification}")
                    clinvar_export.update_classification(new_classification_modification)
                    has_changes = True
            if has_changes:
                clinvar_export.save()
            else:
                log.append(f"No changes for CE_{clinvar_export.pk}")

        else:
            if self.new_classification_modification and self.condition_umbrella:
                ce = ClinVarExport.new_condition(self.clinvar_allele, condition=self.condition_umbrella, candidate=self.new_classification_modification)
                ce.save()
                log.append(f"Creating new clinvar export CE_{ce.pk} {self.condition_umbrella.summary} {self.new_classification_modification}")
        return log

    @staticmethod
    def from_clinvar_export(clinvar_allele: ClinVarAllele, clinvar_export: ClinVarExport) -> 'ClinVarExportStub':
        return ClinVarExportStub(
            clinvar_allele=clinvar_allele,
            new_condition_umbrella=None,
            new_classification_modification=None,
            clinvar_export=clinvar_export
        )


@dataclass
class ClinVarAlleleLog:
    clinvar_allele: ClinVarAllele
    log: list[str]


class ClinVarAlleleExportManager:

    def __init__(self, clinvar_allele: ClinVarAllele):
        self.clinvar_allele = clinvar_allele
        existing_alleles = ClinVarExport.objects.filter(clinvar_allele=clinvar_allele).order_by('-pk')
        self.stubs = [ClinVarExportStub.from_clinvar_export(self.clinvar_allele, clinvar_export) for clinvar_export in existing_alleles]

    def add_classification_modification(self, classification_modification: ClassificationModification):
        condition = classification_modification.classification.condition_resolution_obj
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
                if stub_score is not None and (min_score is None or (stub_score < min_score)):
                    min_score = stub_score
                    use_stub = stub

        if use_stub:
            use_stub.assign_if_newer(classification_modification)
        else:
            # now see if we extend a stub
            min_score = None
            for stub in self.stubs:
                stub_score = stub.condition_umbrella.same_or_more_specific_step_count(condition)
                if stub_score is not None and (min_score is None or (stub_score < min_score)):
                    min_score = stub_score
                    use_stub = stub
            if use_stub:
                use_stub.new_condition_umbrella = condition
                use_stub.assign_if_newer(classification_modification)
            else:
                # won't fit into existing stub or even when extending it, have to create a new one
                new_stub = ClinVarExportStub(clinvar_allele=self.clinvar_allele, new_condition_umbrella=condition, new_classification_modification=classification_modification)
                self.stubs.append(new_stub)

    def apply(self) -> ClinVarAlleleLog:
        complete_log = []
        for stub in self.stubs:
            complete_log += stub.apply()
        return ClinVarAlleleLog(clinvar_allele=self.clinvar_allele, log=complete_log)


def allele_origin_to_clinvar_export_types(bucket: AlleleOriginBucket) -> list[ClinVarExportTypeBucket]:
    if bucket == AlleleOriginBucket.GERMLINE:
        return [ClinVarExportTypeBucket.GERMLINE]
    elif bucket == AlleleOriginBucket.SOMATIC:
        return [ClinVarExportTypeBucket.ONCOGENIC, ClinVarExportTypeBucket.CLINICAL_IMPACT]


def clinvar_export_type_to_allele_origin(bucket: ClinVarExportTypeBucket) -> AlleleOriginBucket:
    if bucket == ClinVarExportTypeBucket.GERMLINE:
        return AlleleOriginBucket.GERMLINE
    else:
        return AlleleOriginBucket.SOMATIC


class ClinVarExportManager:

    def __init__(self, clinvar_key: ClinVarKey):
        self.clinvar_key = clinvar_key

    def base_queryset(self) -> QuerySet[ClassificationModification]:
        labs = list(self.clinvar_key.lab_set.all())
        return ClassificationModification.objects.filter(
            classification__withdrawn=False,
            classification__share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
            classification__lab__in=labs,
            is_last_published=True,
            classification__allele_info__allele__isnull=False
        ).select_related(
            'classification',
            'classification__allele_info__allele',
            'classification__lab',
        ).order_by('classification__allele_info__allele_id', 'classification__lab__clinvar_key_id')

    def run_whole_key(self) -> list[str]:
        complete_log = []
        processed_clinvar_alleles: set[int] = set()
        for allele_origin_bucket in {AlleleOriginBucket.GERMLINE, AlleleOriginBucket.SOMATIC}:
            query = self.base_queryset().filter(classification__allele_origin_bucket=allele_origin_bucket).iterator()
            for allele, classifications_by_allele in itertools.groupby(query, lambda cm: cm.classification.allele_object):
                classifications_by_allele = list(classifications_by_allele)
                for export_type_bucket in allele_origin_to_clinvar_export_types(allele_origin_bucket):
                    clinvar_allele_log = self.run_for(allele, export_type_bucket, classifications_by_allele)
                    if clinvar_allele_log:
                        clinvar_allele = clinvar_allele_log.clinvar_allele
                        complete_log.append(f"Processing ClinVarAllele {clinvar_allele.pk} {clinvar_allele.allele} {clinvar_allele.get_clinvar_export_bucket_display()}")
                        processed_clinvar_alleles.add(clinvar_allele_log.clinvar_allele.pk)
                        complete_log += clinvar_allele_log.log

        no_valid_candidates = ClinVarAllele.objects.filter(clinvar_key=self.clinvar_key).exclude(pk__in=processed_clinvar_alleles)
        for clinvar_allele in no_valid_candidates:
            complete_log.append(f"Found zero classifications for ClinVarAllele {clinvar_allele.pk}")
            complete_log += ClinVarAlleleExportManager(clinvar_allele).apply().log
        return complete_log

    def run_for_single(self, clinvar_allele: ClinVarAllele):
        if clinvar_allele.clinvar_key != self.clinvar_key:
            raise ValueError("Provided clinvar_allele is for a different ClinVarKey than this export manager")

        clinvar_allele_log = self.run_for(allele=clinvar_allele.allele, export_type=clinvar_allele.clinvar_export_bucket)
        return clinvar_allele_log.log or []
        # it's a little inefficient going from ClinVarAllele to its unique parts back to a ClinVarAllele

    def run_for(self, allele: Allele, export_type: ClinVarExportTypeBucket, modifications: Optional[list[ClassificationModification]] = None) -> Optional[ClinVarAlleleLog]:
        if modifications is None:
            modifications = list(self.base_queryset().filter(classification__allele_info__allele=allele, classification__allele_origin_bucket=clinvar_export_type_to_allele_origin(export_type)))

        if modifications:
            modifications = [mod for mod in modifications if mod.classification.condition_resolution_obj]
            sorted_by_condition = sorted(modifications, key=lambda c: c.classification.condition_resolution_obj)
            condition: ConditionResolved
            clinvar_allele_manager: Optional[ClinVarAlleleExportManager] = None
            for condition, grouped_by_condition in itertools.groupby(sorted_by_condition,
                                                                     lambda c: c.classification.condition_resolution_obj):
                if condition.terms:
                    if not clinvar_allele_manager:
                        clinvar_allele, _ = ClinVarAllele.objects.get_or_create(
                            allele=allele,
                            clinvar_export_bucket=export_type
                        )
                        clinvar_allele_manager = ClinVarAlleleExportManager(clinvar_allele)

                    latest_for_condition = max(grouped_by_condition, key=lambda c: c.curated_date_check)
                    clinvar_allele_manager.add_classification_modification(latest_for_condition)

            if clinvar_allele_manager:
                output = clinvar_allele_manager.apply()
                clinvar_allele = clinvar_allele_manager.clinvar_allele
                clinvar_allele.last_evaluated = now()
                clinvar_allele.save()
                return output
