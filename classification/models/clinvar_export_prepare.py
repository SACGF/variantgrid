import itertools
from collections import defaultdict
from typing import List, Optional, Set, Iterable

from django.utils import timezone
from django.utils.timezone import now

from classification.enums import ShareLevel, AlleleOriginBucket
from classification.models import ClinVarAllele, ClassificationModification, ClinVarExport, \
    ConditionResolved, ClinVarExportStatus, ClinVarExportTypeBucket, ClinVarExportDeleteStatus
from classification.models.abstract_utils import ConsolidatingMerger
from library.utils import pretty_collection
from snpdb.lab_picker import LabPickerData
from snpdb.models import Allele, ClinVarKey, Lab

"""
This code creates or updates ClinVarExports with the latest data.
It finds the best candidate classification for each lab / allele / condition combo,
sees if there's already a ClinVar entry for it, if so updates, if not create
"""


class ClassificationModificationCandidate:

    def __init__(self,
                 modification: ClassificationModification,
                 allele_origin_bucket: Optional[AlleleOriginBucket] = None,
                 condition_umbrella: Optional[ConditionResolved] = None,
                 failed_candidates: Optional[Set[ClassificationModification]] = None):
        self.modification = modification
        self.allele_origin_bucket = allele_origin_bucket or modification.classification.allele_origin_bucket
        self.condition_umbrella: ConditionResolved = condition_umbrella or modification.classification.condition_resolution_obj.as_mondo_if_possible()
        self.failed_candidates: Set[ClassificationModification] = failed_candidates or set()

        if self.condition_umbrella is None or not bool(self.condition_umbrella.terms):
            raise ValueError("Candidate must have a resolved condition associated with it")

        if self.allele_origin_bucket not in (AlleleOriginBucket.GERMLINE, AlleleOriginBucket.SOMATIC):
            raise ValueError(f"Candidate must have an allele origin of germline or somatic, not {modification.classification.allele_origin_bucket}")


class ClinVarConsolidatingMerger(ConsolidatingMerger[ClinVarExport, ClassificationModificationCandidate]):

    def __init__(self,
                 clinvar_allele: ClinVarAllele,
                 allele_origin_bucket: AlleleOriginBucket,
                 force_update: bool = True):
        self.clinvar_allele = clinvar_allele
        self.allele_origin_bucket = allele_origin_bucket
        self.log: List[str] = []
        self.force_update = force_update
        super().__init__()

    def retrieve_established(self) -> Set[ClinVarExport]:
        # retrieve existing ClinVarExport records excluding those deleted or marked for deletion
        # as we want to delete them based on their SCVs more than anything else
        return set(ClinVarExport.objects.filter(clinvar_allele=self.clinvar_allele).filter(delete_status=ClinVarExportDeleteStatus.LIVE_RECORD))

    def establish_new_candidate(self, new_candidate: ClassificationModificationCandidate) -> ClinVarExport:
        self.log.append(f"Created export record for {new_candidate.condition_umbrella.summary} : {new_candidate.modification.id_str}")
        return ClinVarExport.new_condition(
            clinvar_allele=self.clinvar_allele,
            condition=new_candidate.condition_umbrella,
            candidate=new_candidate.modification
        )

    def combine_candidates_if_possible(self, candidate_1: ClassificationModificationCandidate, candidate_2: ClassificationModificationCandidate) -> Optional[ClassificationModificationCandidate]:
        if general_condition := ConditionResolved.more_general_term_if_related(candidate_1.condition_umbrella, candidate_2.condition_umbrella):
            mod_1 = candidate_1.modification
            mod_2 = candidate_2.modification
            failed_candidates = candidate_1.failed_candidates.union(candidate_2.failed_candidates)

            best_modification: ClassificationModification
            # we want the biggest date
            if mod_1.curated_date_check < mod_2.curated_date_check:
                best_modification = mod_2
                failed_candidates.add(mod_1)
            else:
                best_modification = mod_1
                failed_candidates.add(mod_2)
            return ClassificationModificationCandidate(
                modification=best_modification,
                allele_origin_bucket=self.allele_origin_bucket,
                condition_umbrella=general_condition,
                failed_candidates=failed_candidates)
        else:
            # these candidates aren't compatible, keep them separate
            return None

    def merge_into_established_if_possible(self, established: ClinVarExport, new_candidate: Optional[ClassificationModificationCandidate]) -> bool:
        """
        Maps an existing group to a condition group
        """
        if new_candidate := new_candidate:
            if new_candidate.condition_umbrella.is_same_or_more_specific(established.condition_resolved):
                if established.classification_based_on == new_candidate.modification:
                    # no change, but still return True to indicate we've found a match
                    if self.force_update:
                        established.update()  # clinvar_key config could have changed, or method for generating JSON could have changed
                    self.log.append(
                        f"No change to existing export record for {established.condition_resolved.summary} : {new_candidate.modification.id_str}")
                    return True
                else:
                    #  can merge, so lets do it
                    # no need to update condition in an hour
                    # retrieve_established.condition_resolved = new_condition.condition_umbrella
                    #  TODO, include (or print) debug info about other candidates
                    self.log.append(
                        f"Updating export record for {established.condition_resolved.summary} : {new_candidate.modification.id_str}")
                    established.update_classification(new_candidate.modification)
                    return True
        else:
            # new_candidate is None
            if not established.clinvarexportsubmission_set.exists():
                self.log.append(f"Export record for {established.condition_resolved.summary} has no classification or history, deleting")
                established.delete()

            elif established.classification_based_on is None:
                self.log.append(f"No change to existing export record for {established.condition_resolved.summary} : None")
                # no change
            else:
                self.log.append(f"Updating export record for {established.condition_resolved.summary} : None")
                established.update_classification(None)

        return False


ClinVarAlleleExportLog = List[str]


class ClinvarExportPrepare:

    @staticmethod
    def update_export_records(perspective: Optional[LabPickerData] = None):
        clinvar_keys: Set[ClinVarKey]
        if perspective:
            clinvar_keys = {lab.clinvar_key for lab in perspective.selected_labs if lab.clinvar_key}
        else:
            clinvar_keys = set(ClinVarKey.objects.all())
        if clinvar_keys:
            return ClinvarExportPrepare.update_export_records_for_keys(clinvar_keys)

    @staticmethod
    def _has_condition(cm: ClassificationModification):
        if resolved_condition := cm.classification.condition_resolution_obj:
            if len(resolved_condition.terms) >= 1:
                return True
        return False

    @staticmethod
    def process_allele(
            clinvar_key: ClinVarKey,
            allele: Allele,
            clinvar_export_bucket: ClinVarExportTypeBucket,
            modifications: Iterable[ClassificationModification]) -> List[str]:

        clinvar_allele, _ = ClinVarAllele.objects.get_or_create(
            clinvar_key=clinvar_key,
            allele=allele,
            clinvar_export_bucket=clinvar_export_bucket
        )
        clinvar_merger = ClinVarConsolidatingMerger(clinvar_allele, allele_origin_bucket=clinvar_export_bucket)

        def is_valid_candidate_for(mod: ClassificationModification, bucket: ClinVarExportTypeBucket):
            if bucket == ClinVarExportTypeBucket.GERMLINE:
                return mod.classification.allele_origin_bucket == AlleleOriginBucket.GERMLINE
            elif bucket in {ClinVarExportTypeBucket.ONCOGENIC, ClinVarExportTypeBucket.CLINICAL_IMPACT}:
                return mod.classification.allele_origin_bucket == AlleleOriginBucket.SOMATIC

        no_condition_count = 0
        for mod in modifications:
            if is_valid_candidate_for(mod, clinvar_export_bucket):
                if ClinvarExportPrepare._has_condition(mod):
                    clinvar_merger.add_new_candidate(ClassificationModificationCandidate(mod))
                else:
                    no_condition_count += 1
        clinvar_merger.consolidate()

        total = clinvar_allele.clinvarexport_set.count()
        in_error = clinvar_allele.clinvarexport_set.filter(status=ClinVarExportStatus.IN_ERROR).count()
        valid = total - in_error

        clinvar_allele.classifications_missing_condition = no_condition_count
        clinvar_allele.submissions_valid = valid
        clinvar_allele.submissions_invalid = in_error
        clinvar_allele.last_evaluated = now()
        clinvar_allele.save()

        log = clinvar_merger.log
        log = [f"{clinvar_key} : {entry}" for entry in log if not entry.startswith("No change")]

        if no_condition_count:
            if no_condition_count == 1:
                log.append(f"{clinvar_key} : 1 shared classification for {clinvar_export_bucket} {allele} doesn't have resolved conditions")
            else:
                log.append(f"{clinvar_key} : {no_condition_count} shared classifications for {clinvar_export_bucket} {allele} don't have resolved conditions")

        return log

    @staticmethod
    def update_export_records_for_keys(clinvar_keys: Set[ClinVarKey]) -> ClinVarAlleleExportLog:
        # work on clinvar keys, not on labs, as a user could have access to one lab but the clinvar key might be for 2
        # and a clinvar key has to get all labs updated or none, can't deal with partial
        clinvar_labs = Lab.objects.filter(clinvar_key__in=clinvar_keys)

        all_classifications = ClassificationModification.objects.filter(
            classification__withdrawn=False,
            classification__share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS,
            classification__lab__in=clinvar_labs,
            is_last_published=True,
            classification__allele__isnull=False
        ).select_related(
            'classification',
            'classification__allele_info__allele',
            'classification__lab',
            'classification__lab__clinvar_key'
        ).order_by('classification__allele_id', 'classification__lab__clinvar_key_id').iterator(chunk_size=1000)

        combined_log = []

        # need to keep track of which allele, clinvar_key combos we've seen,
        # so we can check all other alleles to make sure they haven't lost candidates
        clinvar_export_buckets = (ClinVarExportTypeBucket.GERMLINE, ClinVarExportTypeBucket.ONCOGENIC, ClinVarExportTypeBucket.CLINICAL_IMPACT)

        allele_origin_clinvar_key_processed_alleles = {}
        for clinvar_export_bucket in clinvar_export_buckets:
            allele_origin_clinvar_key_processed_alleles[clinvar_export_bucket] = defaultdict(list)

        # FIND ALL THE GOOD CANDIDATES FOR clinvar key / allele combo

        # loop through classifications for the lab grouping by allele
        for allele, classifications_for_allele in itertools.groupby(all_classifications, lambda cm: cm.classification.allele_object):
            # now loop that grouping by clinvar key
            for clinvar_key, alleles_for_clinvar_key in itertools.groupby(classifications_for_allele, lambda cm: cm.classification.lab.clinvar_key):
                alleles_for_clinvar_key = list(alleles_for_clinvar_key)
                for clinvar_export_bucket in clinvar_export_buckets:
                    clinvar_allele_classifications_allele_origin = [cac for cac in alleles_for_clinvar_key]
                    combined_log += ClinvarExportPrepare.process_allele(
                        clinvar_key=clinvar_key,
                        allele=allele,
                        clinvar_export_bucket=clinvar_export_bucket,
                        modifications=clinvar_allele_classifications_allele_origin
                    )
                    allele_origin_clinvar_key_processed_alleles[clinvar_export_bucket][clinvar_key].append(allele)

        # FIND ALL THE PRE-EXISTING RECORDS THAT NO LONGER HAVE ANY CANDIDATES

        # loop through ClinVarAlleles for clinvar key that we didn't find by looking at all the non-withdrawn classifications for that lab
        # i.e. these will be the alleles that
        for clinvar_key in clinvar_keys:
            for clinvar_export_bucket in clinvar_export_buckets:
                clinvar_key_to_processed_alleles = allele_origin_clinvar_key_processed_alleles[clinvar_export_bucket]
                for cva in ClinVarAllele.objects.filter(
                        clinvar_key=clinvar_key,
                        clinvar_export_bucket=clinvar_export_bucket
                    ).exclude(allele__in=clinvar_key_to_processed_alleles.get(clinvar_key, [])).select_related('allele').iterator():
                    combined_log += ClinvarExportPrepare.process_allele(
                        clinvar_key=clinvar_key,
                        allele=cva.allele,
                        clinvar_export_bucket=clinvar_export_bucket,
                        modifications=[])

        # WRAP THINGS UP

        completed_date = timezone.now()
        for clinvar_key in clinvar_keys:
            clinvar_key.last_full_run = completed_date
            clinvar_key.save()

        if not combined_log:
            combined_log.append(f"No changes detected to ClinVarExports for {pretty_collection(clinvar_keys)}")
        return combined_log
