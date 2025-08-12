from dataclasses import dataclass, field
from functools import cached_property
from typing import Iterable, Optional

from django.core.management import BaseCommand

from classification.enums import TestingContextBucket, AlleleOriginBucket, ShareLevel, SpecialEKeys, ConflictType
from classification.models import classification_flag_types, Classification, ClassificationModification, \
    ClassificationSummaryCacheDict, ClassificationSummaryCalculator, ConflictKey, Conflict, \
    AlleleOriginGrouping
from classification.services.conflict_services import ConflictDataRow, ClinSigCalculator, OncPathCalculator
from flags.models import FlagComment
from datetime import datetime, timedelta

from library.utils import IterableStitcher
from snpdb.models import Allele, Lab


@dataclass
class ClassificationChange:
    """
    Represents a change or state of a Classification, be it a new submission or a withdrawn/unwithdraw flag
    The withdraw/unwithdraw will only have partial information but a submission should have everything required.
    The withdraw/unwithdraw ClassificationChange will have to be amended to a publish
    """
    classification_id: int
    change_date: datetime
    lab: Optional[Lab] = None
    allele: Optional[Allele] = None
    withdrawn: Optional[bool] = None
    summary_data: Optional[ClassificationSummaryCacheDict] = None
    share_level: Optional[ShareLevel] = None

    current_grouping: Optional['AlleleOriginGroupingKey'] = None
    """
    The grouping this was last put into, there may have been subsequence calls to merge which will change what the new grouping should be
    """

    def merge(self, other: 'ClassificationChange') -> None:
        """
        Updates the existing ClassificationChange with data from "other" whenever "other has data
        :param other: A more up to date ClassificationChange
        """
        if self.classification_id != other.classification_id:
            raise ValueError("Can't merge records with different IDs")

        self.change_date = other.change_date
        if other.withdrawn is not None:
            self.withdrawn = other.withdrawn
        if other.allele is not None:
            self.allele = other.allele
        if other.summary_data is not None:
            self.summary_data = other.summary_data
        if other.share_level is not None:
            self.share_level = other.share_level
        if other.lab is not None:
            self.lab = other.lab

    def __lt__(self, other) -> bool:
        # note this date is based on when the action happened within VariantGrid, not the dates manually entered on Classifications
        return self.change_date < other.change_date

    def calc_grouping_key(self) -> 'AlleleOriginGroupingKey':
        """
        Which AlleleOriginGroupingKey should this go into in its current state.
        Note that withdrawn records don't go into any grouping
        """
        if self.withdrawn:
            return None
        if not (self.allele and self.summary_data and self.share_level):
            # occasionally there might still be withdraws for
            return None
        allele_origin_bucket = self.summary_data.get("allele_origin_bucket")
        testing_context_bucket = TestingContextBucket(self.summary_data.get("somatic", {}).get("testing_context_bucket") or "G")
        return AlleleOriginGroupingKey(
            allele=self.allele,
            allele_origin_bucket=AlleleOriginBucket(allele_origin_bucket),
            testing_context_bucket=testing_context_bucket,
            tumor_type_category=self.summary_data.get("tumor_type_category")
        )


class ClassificationSummaryCalculatorHistoric(ClassificationSummaryCalculator):
    """
    Extend the default ClassificationSummaryCalculator so this works better on historic values
    """

    @property
    def allele_origin_bucket(self) -> AlleleOriginBucket:
        # normally this is just taken from the Classification
        return AlleleOriginBucket.bucket_for_allele_origin(self.cm.get(SpecialEKeys.ALLELE_ORIGIN))

    @cached_property
    def pending_classification_value(self) -> Optional[str]:
        # default pending calculation works on latest flag value rather than current value
        return None


@dataclass(frozen=True)
class AlleleOriginGroupingKey:
    """
    Used to identify the context that a classification belongs to at any point
    See AlleleGroupingInfo for how records that share the grouping key are treated
    """
    allele: Allele
    allele_origin_bucket: AlleleOriginBucket
    testing_context_bucket: TestingContextBucket
    tumor_type_category: Optional[str]


@dataclass(frozen=True)
class LabShare:
    lab: Lab
    share_level: ShareLevel


@dataclass
class AlleleOriginGroupingInfo:
    grouping_key: AlleleOriginGroupingKey
    classification_changes: dict[int, ClassificationChange] = field(default_factory=dict)

    def remove_classification(self, classification_id: int):
        self.classification_changes.pop(classification_id)

    def add_classification(self, classification_change: ClassificationChange):
        self.classification_changes[classification_change.classification_id] = classification_change

    def calculate_and_log_history(self, override_date: datetime):
        """
        Calculate the would be classification groupings, then have them see if they need to add to the conflict history
        """

        # grab the latest classification from each lab/share level combo
        lab_id_to_latest: dict[LabShare, ClassificationChange] = {}
        for value in self.classification_changes.values():
            lab_share = LabShare(value.lab, value.share_level)
            existing_value_for_lab = lab_id_to_latest.get(lab_share)
            if existing_value_for_lab is None or existing_value_for_lab < value:
                lab_id_to_latest[lab_share] = value

        oncpath_data: list[ConflictDataRow] = []
        clinsig_data: list[ConflictDataRow] = []
        for latest_summary in lab_id_to_latest.values():
            lab_id = latest_summary.lab.pk
            summary_info = latest_summary.summary_data
            share_level = latest_summary.share_level

            oncpath_data.append(ConflictDataRow.from_data(
                row=summary_info.get("pathogenicity", {}),
                lab_id=lab_id,
                share_level=share_level))
            clinsig_data.append(ConflictDataRow.from_data(
                row=summary_info.get("somatic", {}),
                lab_id=lab_id,
                share_level=share_level))

            if self.grouping_key.allele_origin_bucket == AlleleOriginBucket.SOMATIC:
                conflict_key = ConflictKey(
                    conflict_type=ConflictType.CLIN_SIG,
                    allele_origin_bucket=self.grouping_key.allele_origin_bucket,
                    testing_context_bucket=self.grouping_key.testing_context_bucket,
                    tumor_type_category=self.grouping_key.tumor_type_category
                )
                ClinSigCalculator(self.grouping_key.allele, conflict_key, clinsig_data).log_history(override_date)

        conflict_key = ConflictKey(
            conflict_type=ConflictType.ONCPATH,
            allele_origin_bucket=self.grouping_key.allele_origin_bucket,
            testing_context_bucket=self.grouping_key.testing_context_bucket,
            tumor_type_category=self.grouping_key.tumor_type_category
        )
        OncPathCalculator(self.grouping_key.allele, conflict_key, oncpath_data).log_history(override_date)

    @property
    def latest_change_data(self) -> Optional[ClassificationSummaryCacheDict]:
        return max(self.classification_changes.values(), default=None)


@dataclass
class Populate:

    classification_id_to_change: dict[int, ClassificationChange] = field(default_factory=dict)
    """
    A complete list of classifications current state keyed by the classification ID
    """

    grouping_key_to_grouping: dict[AlleleOriginGroupingKey, AlleleOriginGroupingInfo] = field(default_factory=dict)
    """
    Keeps track of all the state of the AlleleOriginGroups
    """

    classification_id_dirty: set[int] = field(default_factory=set)
    """
    All classification IDs that have had a chance since log_and_update was last called
    """

    def apply_update(self, change: ClassificationChange):
        if existing := self.classification_id_to_change.get(change.classification_id):
            existing.merge(change)
        else:
            self.classification_id_to_change[change.classification_id] = change
        self.classification_id_dirty.add(change.classification_id)

    def grouping_for_key(self, allele_origin_grouping_key: AlleleOriginGroupingKey) -> AlleleOriginGroupingInfo:
        origin_grouping_info = self.grouping_key_to_grouping.get(allele_origin_grouping_key)
        if not origin_grouping_info:

            origin_grouping_info = AlleleOriginGroupingInfo(grouping_key=allele_origin_grouping_key)
            self.grouping_key_to_grouping[allele_origin_grouping_key] = origin_grouping_info
        return origin_grouping_info

    def log_and_update(self, override_date: datetime):
        """
        Work out which AlleleOriginGroupings have changed (via applying the changes of the modified classifications
        since this method was last called - then call calculate_and_log_history to make new entries (where apprioriate)
        to the Conflict history
        :param override_date: The date the conflict history record should be changed to
        :return:
        """
        dirty_groups: set[AlleleOriginGroupingKey] = set()
        for dirty_id in self.classification_id_dirty:
            cc = self.classification_id_to_change[dirty_id]
            last_grouping = cc.current_grouping
            next_grouping = cc.calc_grouping_key()
            if last_grouping == next_grouping:
                if next_grouping:
                    dirty_groups.add(next_grouping)
            else:
                if last_grouping:
                    dirty_groups.add(last_grouping)
                    self.grouping_for_key(last_grouping).remove_classification(dirty_id)
                cc.current_grouping = cc.calc_grouping_key()
                if cc.current_grouping:
                    self.grouping_for_key(cc.current_grouping).add_classification(cc)
                    dirty_groups.add(cc.current_grouping)

        self.classification_id_dirty.clear()
        for dirty_group in dirty_groups:
            self.grouping_for_key(dirty_group).calculate_and_log_history(override_date)

    @cached_property
    def flag_collection_id_to_classification(self) -> dict[int, int]:
        """
        Creates a dictory that links flag_collection_id back to classification_id (so we can map classification flags back to this)
        """
        flag_collection_map = {}
        for classification_id, flag_collection_id in Classification.objects.values_list("pk", "flag_collection_id").iterator():
            flag_collection_map[flag_collection_id] = classification_id
        return flag_collection_map

    def flag_collection_to_classification_id(self, flag_collection_id: int) -> Optional[int]:
        return self.flag_collection_id_to_classification.get(flag_collection_id)

    def withdrawn_iterable(self) -> Iterable[ClassificationChange]:
        flag_qs = FlagComment.objects.filter(flag__flag_type=classification_flag_types.classification_withdrawn) \
            .order_by("created") \
            .values_list("resolution__status", "created", "flag__collection_id")

        for res_status, created, flag_collection_id in flag_qs.iterator():
            if classification_id := self.flag_collection_to_classification_id(flag_collection_id):
                withdrawn = res_status == "O"
                yield ClassificationChange(
                    classification_id=classification_id,
                    change_date=created,
                    withdrawn=withdrawn
                )

    def classification_iterable(self) -> Iterable[ClassificationChange]:
        for cm in ClassificationModification.objects.filter(
            published=True,
            classification__allele_info__allele__isnull=False
        ).order_by("created").iterator():
            summary = ClassificationSummaryCalculator(cm).cache_dict()
            yield ClassificationChange(
                classification_id=cm.classification_id,
                change_date=cm.modified,
                lab=cm.lab,
                allele=cm.classification.allele_object,
                summary_data=summary,
                share_level=cm.share_level_enum
            )

    def stitch_iterable(self):
        return IterableStitcher[ClassificationChange](
            iterables=[
                self.withdrawn_iterable(),
                self.classification_iterable()
            ]
        )

    def full_clean(self):
        Conflict.objects.all().delete()

    def run(self, time_delta: timedelta = timedelta(days=1)):
        end_time: Optional[datetime] = None
        for cc in self.stitch_iterable():
            if end_time is None:
                end_time = cc.change_date.replace(hour=0, minute=0, second=0, microsecond=0) + time_delta
            if cc.change_date > end_time:
                self.log_and_update(end_time)
                end_time = end_time + time_delta
            self.apply_update(cc)
        self.log_and_update(end_time)

    def run_latest(self):
        # now do latest
        for aog in AlleleOriginGrouping.objects.iterator():
            aog.update()


class Command(BaseCommand):

    def add_arguments(self, parser):
        pass

    def handle(self, *args, **options):
        populate = Populate()
        populate.full_clean()
        populate.run()
        populate.run_latest()


