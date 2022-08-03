import itertools
from collections import defaultdict
from typing import List, Optional, Set, Iterable

from django.utils import timezone
from django.utils.timezone import now

from classification.enums import ShareLevel
from classification.models import ClinVarAllele, ClassificationModification, ClinVarExport, \
    ConditionResolved, ClinVarExportStatus
from classification.models.abstract_utils import ConsolidatingMerger
from library.utils import pretty_collection
from snpdb.lab_picker import LabPickerData
from snpdb.models import Allele, ClinVarKey, Lab


class ClassificationModificationCandidate:

    def __init__(self,
                 modification: ClassificationModification,
                 condition_umbrella: Optional[ConditionResolved] = None,
                 failed_candidates: Optional[Set[ClassificationModification]] = None):
        self.modification = modification
        self.condition_umbrella: ConditionResolved = condition_umbrella or modification.classification.condition_resolution_obj.as_mondo_if_possible()
        self.failed_candidates: Set[ClassificationModification] = failed_candidates or set()

        if self.condition_umbrella is None or not bool(self.condition_umbrella.terms):
            raise ValueError("Candidate must have a resolved condition associated with it")


class ClinVarConsolidatingMerger(ConsolidatingMerger[ClinVarExport, ClassificationModificationCandidate]):

    def __init__(self, clinvar_allele: ClinVarAllele, force_update: bool = True):
        self.clinvar_allele = clinvar_allele
        self.log: List[str] = []
        self.force_update = force_update
        super().__init__()

    def retrieve_established(self) -> Set[ClinVarExport]:
        return set(ClinVarExport.objects.filter(clinvar_allele=self.clinvar_allele))

    def establish_new_candidate(self, new_candidate: ClassificationModificationCandidate) -> ClinVarExport:
        self.log.append(f"Created export record for {new_candidate.condition_umbrella.summary} : {new_candidate.modification.id_str}")
        return ClinVarExport.new_condition(clinvar_allele=self.clinvar_allele, condition=new_candidate.condition_umbrella,
                                           candidate=new_candidate.modification)

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
            return ClassificationModificationCandidate(modification=best_modification, condition_umbrella=general_condition, failed_candidates=failed_candidates)
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
    def process_allele(clinvar_key: ClinVarKey, allele: Allele, modifications: Iterable[ClassificationModification]) -> List[str]:
        clinvar_allele, _ = ClinVarAllele.objects.get_or_create(clinvar_key=clinvar_key, allele=allele)
        clinvar_merger = ClinVarConsolidatingMerger(clinvar_allele)

        no_condition_count = 0
        for mod in modifications:
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
                log.append(f"{clinvar_key} : 1 shared classification for {allele} doesn't have resolved conditions")
            else:
                log.append(f"{clinvar_key} : {no_condition_count} shared classifications for {allele} don't have resolved conditions")

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
            allele__isnull=False
        ).select_related('classification', 'classification__allele', 'classification__lab', 'classification__lab__clinvar_key').order_by('classification__allele_id', 'classification__lab__clinvar_key_id')

        combined_log = list()

        # need to keep track of which allele, clinvar_key combos we've seen
        # so we can check all other alleles to make sure they haven't lost candidates
        clinvar_key_to_processed_alleles = defaultdict(list)
        for allele, allele_classifications in itertools.groupby(all_classifications, lambda cm: cm.classification.allele):
            for clinvar_key, clinvar_allele_classifications in itertools.groupby(allele_classifications, lambda cm: cm.classification.lab.clinvar_key):
                combined_log += ClinvarExportPrepare.process_allele(clinvar_key=clinvar_key, allele=allele, modifications=clinvar_allele_classifications)
                clinvar_key_to_processed_alleles[clinvar_key].append(allele)

        # loop through ClinVarAlleles for clinvar key that we didn't find by looking at all the non-withdrawn classifications for that lab
        # i.e. these will be the alleles that
        for clinvar_key in clinvar_keys:
            for allele_id in ClinVarAllele.objects.filter(clinvar_key=clinvar_key).exclude(allele__in=clinvar_key_to_processed_alleles.get(clinvar_key, list())).values_list('allele', flat=True):
                allele = Allele.objects.get(pk=allele_id)
                combined_log += ClinvarExportPrepare.process_allele(clinvar_key=clinvar_key, allele=allele, modifications=list())

        completed_date = timezone.now()
        for clinvar_key in clinvar_keys:
            clinvar_key.last_full_run = completed_date
            clinvar_key.save()

        if not combined_log:
            combined_log.append(f"No changes detected to ClinVarExports for {pretty_collection(clinvar_keys)}")
        return combined_log
