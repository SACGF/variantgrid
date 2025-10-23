import itertools
from abc import ABC, abstractmethod
from collections import defaultdict, Counter
from dataclasses import dataclass
from enum import IntEnum
from functools import cached_property
from html import escape
from typing import Optional, Iterable, Callable, Collection, Tuple
from datetime import datetime
from django.contrib.auth.models import User
from django.db import transaction
from django.urls import reverse
from django.utils.safestring import SafeString
from more_itertools.recipes import partition
from classification.enums import AlleleOriginBucket, ConflictSeverity, ShareLevel, ConflictType, TestingContextBucket, \
    SpecialEKeys
from classification.models import AlleleOriginGrouping, EvidenceKeyMap, Conflict, ConflictKey, ConflictHistory, \
    ConflictLab, ClassificationGrouping, ClassificationModification, ConflictNotification, ConflictNotificationStatus, \
    ConflictNotificationRun, AlleleGrouping, DiscordanceReportTriageStatus
from genes.hgvs import CHGVS
from library.utils import strip_json, first, sort_and_group_by, RowSpanTable, RowSpanCellValue, LinkType
from snpdb.models import Allele, Lab

"""

@dataclass(frozen=True)
class ClinSigData(ConflictDataRow):
    clin_sig: str
    amp_level: int

    @staticmethod
    def from_data(lab_id: int, share_level: ShareLevel, row: dict) -> 'ClinSigData':
        return ClinSigData(lab_id=lab_id, share_level=share_level, clin_sig=row.get("clinical_significance"), amp_level=row.get("amp_level"))

    def to_json(self) -> dict:
        return {
            "lab_id": self.lab_id,
            "clinical_significance": self.clin_sig,
            "amp_level": self.amp_level
        }
"""


@dataclass
class ConflictDataRow:
    lab_id: int
    share_level: ShareLevel
    classification: Optional[str]
    clinical_significance: Optional[str]
    amp_level: Optional[str]
    exclude: bool = False
    c_hgvs: Optional[CHGVS] = None
    message: Optional[str] = None

    @cached_property
    def lab(self):
        # Lab uses Caching Object Manager
        return Lab.objects.get(id=self.lab_id)

    @cached_property
    def value_label(self) -> str:
        value_parts = []
        if classification := self.classification:
            value_parts.append(EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(classification))
        if clinical_significance := self.clinical_significance:
            value_parts.append(EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value(clinical_significance))
        if amp_level := self.amp_level:
            value_parts.append(amp_level)
        return " ".join(value_parts)

    def can_view(self, user: User) -> bool:
        return self.share_level.is_discordant_level or (user and self.share_level.has_access(self.lab, user))

    def with_exclude(self, exclude: bool) -> 'ConflictDataRow':
        if self.exclude == exclude:
            return self
        return ConflictDataRow(
            lab_id=self.lab_id,
            share_level=self.share_level,
            classification=self.classification,
            clinical_significance=self.clinical_significance,
            amp_level=self.amp_level,
            exclude=exclude
        )

    def to_json(self) -> dict:
        return strip_json({
            "lab_id": self.lab_id,
            "share_level": self.share_level.value,
            "classification": self.classification,
            "clinical_significance": self.clinical_significance,
            "amp_level": self.amp_level,
            "exclude": self.exclude
        })

    @staticmethod
    def from_json(row: dict) -> 'ConflictDataRow':
        return ConflictDataRow(
            lab_id=row.get("lab_id"),
            share_level=ShareLevel(row.get("share_level")),
            classification=row.get("classification"),
            clinical_significance=row.get("clinical_significance"),
            amp_level=row.get("amp_level"),
            exclude=row.get("exclude")
        )

    @staticmethod
    def from_data(row: dict, lab_id: int, share_level: ShareLevel, c_hgvs: Optional[CHGVS] = None) -> 'ConflictDataRow':
        return ConflictDataRow(
            lab_id=lab_id,
            share_level=share_level,
            classification=row.get("classification"),
            clinical_significance=row.get("clinical_significance"),
            amp_level=row.get("amp_level"),
            c_hgvs=c_hgvs
        )

    @staticmethod
    def to_not_excluded_dict(data_rows: list['ConflictDataRow']) -> dict[Lab, 'ConflictDataRow']:
        conflict_dict: dict[Lab, ConflictDataRow] = {}
        for data_row in data_rows:
            if not data_row.exclude:
                lab = data_row.lab
                if existing := conflict_dict.get(lab):
                    if data_row < existing:
                        continue
                conflict_dict[lab] = data_row
        return conflict_dict

    @property
    def bucket(self) -> Optional[int]:
        if classification := self.classification:
            return EvidenceKeyMap.clinical_significance_to_bucket().get(classification)
        return None

    @staticmethod
    def pathogenic_buckets() -> set[int]:
        # Note, that all of these should be the same bucket
        # but just in case
        path_buckets: set[int] = set()
        clin_sig_buckets = EvidenceKeyMap.clinical_significance_to_bucket()
        for value in ("LP", "P", "LO", "O"):
            if bucket := clin_sig_buckets.get(value):
                path_buckets.add(bucket)
        return path_buckets

    def __lt__(self, other: 'ConflictDataRow') -> bool:
        # sort by labs by default to keep ordering consistent over entries
        if self.lab_id != other.lab_id:
            return self.lab < other.lab
        return self.share_level < other.share_level


class ConflictCalculator(ABC):

    def __init__(self, allele: Allele, conflict_key: ConflictKey, conflict_data: list[ConflictDataRow]):
        self.allele = allele
        self.conflict_key = conflict_key
        self.conflict_data = self.exclude_duplicates(conflict_data)

    @property
    def included_conflict_data(self) -> list[ConflictDataRow]:
        return [cd for cd in self.conflict_data if not cd.exclude]

    @abstractmethod
    def calculate_severity(self) -> ConflictSeverity:
        raise NotImplementedError()

    def should_exclude(self, row: ConflictDataRow) -> bool:
        return not row.share_level.is_discordant_level

    def exclude_duplicates(self, data_list: list[ConflictDataRow]) -> list[ConflictDataRow]:
        by_lab_id: dict[int, list[ConflictDataRow]] = defaultdict(list)
        for row in data_list:
            by_lab_id[row.lab_id].append(row)

        updated_list: list[ConflictDataRow] = []
        # TODO use Lab proper ordering
        # but then need to work out efficient way to load Lab object
        for key, values in by_lab_id.items():
            for index, data_row in enumerate(sorted(list(values), reverse=True)):
                if index == 0 and not self.should_exclude(data_row):
                    updated_list.append(data_row)
                else:
                    updated_list.append(data_row.with_exclude(True))
        return list(sorted(updated_list))

    @cached_property
    def involved_lab_ids(self) -> set[int]:
        lab_ids = set()
        for row in self.conflict_data:
            lab_ids.add(row.lab_id)
        return lab_ids

    @cached_property
    def active_labs(self) -> dict[Lab, bool]:
        lab_ids = defaultdict(lambda: False)
        for row in self.conflict_data:
            lab_ids[row.lab] |= not row.exclude
        return dict(lab_ids)

    def data_json(self) -> dict:
        data_rows = []
        included, excluded = partition(lambda x: x.exclude, sorted(self.conflict_data))
        for row in list(included) + list(excluded):
            data_rows.append(row.to_json())
        return {
            "rows": data_rows
        }

    @transaction.atomic
    def log_history(self, override_date: Optional[datetime] = None):
        allele = self.allele
        conflict_key = self.conflict_key
        conflict: Conflict
        created: bool
        if not bool(self.included_conflict_data):
            # if our calculation is empty (e.g. we worked out there's no contributing labs)
            # don't make the conflict just for it to be empty
            # BUT still need to call log_history in case there was data once and now we need to log
            # it's no longer there
            conflict = Conflict.objects.filter(
                allele=allele,
                conflict_type=conflict_key.conflict_type,
                allele_origin_bucket=conflict_key.allele_origin_bucket,
                testing_context_bucket=conflict_key.testing_context_bucket,
                tumor_type_category=conflict_key.tumor_type_category,
            ).first()
            if not conflict:
                return
            created = False
        else:
            conflict, created = Conflict.objects.get_or_create(
                allele=allele,
                conflict_type=conflict_key.conflict_type,
                allele_origin_bucket=conflict_key.allele_origin_bucket,
                testing_context_bucket=conflict_key.testing_context_bucket,
                tumor_type_category=conflict_key.tumor_type_category,
            )

        latest: Optional[ConflictHistory] = None
        if not created:
            latest = ConflictHistory.objects.filter(conflict=conflict, is_latest=True).first()

        # also need to create ConflictLab
        active_labs = self.active_labs
        ConflictLab.objects.bulk_create([
            ConflictLab(
                conflict=conflict,
                classification_grouping=classification_grouping_for_conflict(conflict, lab),
                lab=lab,
                active=active,
            ) for lab, active in active_labs.items()
        ], update_conflicts=True, update_fields=["active"], unique_fields=["conflict", "lab_id"])
        # TODO have something different from withdrawn
        ConflictLab.objects.filter(conflict=conflict).exclude(lab__in=active_labs.keys()).update(active=False)

        c_hgvs_source = self.included_conflict_data
        if not c_hgvs_source:
            c_hgvs_source = self.conflict_data

        c_hgvs_set = [row.c_hgvs for row in c_hgvs_source]
        if c_hgvs_list := [c_hgvs.to_json_short() for c_hgvs in sorted(set(c_hgvs_set)) if c_hgvs]:
            conflict.meta_data["c_hgvs"] = c_hgvs_list
            conflict.save(update_fields=["meta_data"], update_modified=False)

        data = self.data_json()
        if latest:
            if latest.data == data:
                return  # data is already up to date
            else:
                latest.is_latest = False
                latest.save(update_fields=["is_latest"])

        ch = ConflictHistory(
            conflict=conflict,
            data=data,
            severity=self.calculate_severity(),
            is_latest=True
        )
        ch.save()
        if override_date:
            ch.created = override_date
            ch.modified = override_date
            ch.save(update_fields=["created", "modified"], update_modified=False)

        if not override_date:
            # don't log notifications when filling in the backlog of data
            ConflictCompare(latest, ch).store_to_db()


class ConflictCompareType(IntEnum):
    NO_CHANGE = 0
    BECAME_CONFLICT = 1
    BECAME_AGREED = 2
    CHANGE_OF_LABS = 3

    @property
    def label(self):
        match self:
            case ConflictCompareType.NO_CHANGE: return "No change"
            case ConflictCompareType.CHANGE_OF_LABS: return "Change of labs"
            case ConflictCompareType.BECAME_CONFLICT: return "Now in conflict"
            case ConflictCompareType.BECAME_AGREED: return "Now resolved"


# @dataclass
# class ConflictLabChange:
#     lab: Lab
#     old_value: Optional[ConflictDataRow]
#     new_value: Optional[ConflictDataRow]
#     curation_date: Optional[ClassificationDate]
#
#     @cached_property
#     def label(self):
#         new_value_str: Optional[str] = None
#         old_value_str: Optional[str] = None
#         if old_value := self.old_value:
#             old_value_str = old_value.value_label
#         if new_value := self.new_value:
#             new_value_str = new_value.value_label
#
#         curation_date_str = None
#         if curation_date := self.curation_date:
#             curation_date_str = f" _Curated on {curation_date.relevant_date.date_str}_"
#
#         if old_value_str == new_value_str:
#             return f"**{new_value_str}**{curation_date_str}"
#         elif not old_value_str:
#             return f"**{new_value_str}** (Newly Uploaded)" + curation_date_str
#         elif not new_value_str:
#             return f"{old_value_str} -> Withdrawn"
#         else:
#             return f"{old_value_str} -> **{new_value_str}**" + curation_date_str



@dataclass
class ConflictCompare:
    previous_state: Optional[ConflictHistory]
    current_state: ConflictHistory

    # def conflict_lab_changes(self):
    #     old_data_rows = {}
    #     new_data_rows = ConflictDataRow.to_not_excluded_dict(self.current_state.data_rows())
    #     if previous_state := self.previous_state:
    #         old_data_rows = ConflictDataRow.to_not_excluded_dict(previous_state.data_rows())
    #
    #     all_labs = list(sorted(old_data_rows.keys() | new_data_rows.keys()))
    #     allele_origin_grouping = self.current_state.conflict.allele_origin_grouping()
    #     for lab in all_labs:
    #         old_value = old_data_rows.get(lab)
    #         new_value = new_data_rows.get(lab)
    #         share_level = new_value.share_level if new_value else old_value.share_level
    #
    #         curation_date: ClassificationDate
    #         if allele_origin_grouping:
    #             if classification_grouping := ClassificationGrouping.objects.filter(
    #                 allele_origin_grouping=allele_origin_grouping,
    #                 lab=lab,
    #                 share_level=share_level,
    #             ).first():
    #                 curation_date = classification_grouping.latest_classification_modification.curated_date_check
    #         yield ConflictLabChange(
    #             lab=lab,
    #             old_value=old_data_rows.get(lab),
    #             new_value=new_data_rows.get(lab),
    #             curation_date=curation_date
    #         )

    @property
    def conflict(self) -> Conflict:
        return self.current_state.conflict

    @staticmethod
    def from_conflict_notification(cn: ConflictNotification):
        return ConflictCompare(previous_state=cn.previous_state, current_state=cn.current_state)

    @cached_property
    def notification_change_type(self) -> ConflictCompareType:

        was_reportable = False
        if self.previous_state:
            was_reportable = self.previous_state.severity >= ConflictSeverity.MAJOR

        is_reportable_now = self.current_state.severity >= ConflictSeverity.MAJOR

        if was_reportable != is_reportable_now:
            if was_reportable:
                return ConflictCompareType.BECAME_AGREED
            else:
                return ConflictCompareType.BECAME_CONFLICT
        else:
            if self.previous_state:
                if self.previous_state.involved_lab_ids != self.current_state.involved_lab_ids:
                    return ConflictCompareType.CHANGE_OF_LABS

            return ConflictCompareType.NO_CHANGE

    @property
    def is_notifiable_difference(self):
        return bool(self.notification_change_type)

    def store_to_db(self):
        if ConflictNotification.objects.filter(
                conflict=self.current_state.conflict,
                notification_run__isnull=True).update(current_state=self.current_state):
            # ConflictNotification already existed for this conflict, update its state
            return

        if self.is_notifiable_difference:
            # didn't already exist, only create it if it'll be a notifiable difference
            # as an optimisation
            ConflictNotification.objects.create(
                conflict=self.current_state.conflict,
                previous_state=self.previous_state,
                current_state=self.current_state,
            )

    @property
    def involved_lab_ids(self) -> set[int]:
        all_lab_ids = self.current_state.involved_lab_ids
        if previous_state := self.previous_state:
            all_lab_ids |= previous_state.involved_lab_ids
        return all_lab_ids

    def __lt__(self, other):
        return self.current_state.conflict < other.current_state.conflict


def combine_outstanding_conflict_notifications() -> Optional[ConflictNotificationRun]:
    if ConflictNotification.objects.filter(notification_run__isnull=True).exists():
        notification_run = ConflictNotificationRun.objects.create()
        ConflictNotification.objects.filter(notification_run__isnull=True).update(notification_run=notification_run)
        return notification_run
    return None


class OncPathCalculator(ConflictCalculator):

    def should_exclude(self, row: ConflictDataRow):
        return super().should_exclude(row) or row.bucket is None

    def calculate_severity(self) -> ConflictSeverity:
        oncpath_data = self.included_conflict_data
        if len(oncpath_data) == 0:
            return ConflictSeverity.NO_SUBMISSIONS
        elif len(oncpath_data) == 1:
            return ConflictSeverity.SINGLE_SUBMISSION
        else:
            buckets: set = set()
            values: set = set()
            for row in oncpath_data:
                if bucket := row.bucket:
                    buckets.add(bucket)
                values.add(row.classification)

            if len(buckets) > 1:
                if buckets.intersection(ConflictDataRow.pathogenic_buckets()):
                    return ConflictSeverity.MEDICALLY_SIGNIFICANT
                else:
                    return ConflictSeverity.MAJOR
            if len(values) > 1:
                if len(values) == 2 and "VUS" in values and all(x.startswith("VUS") for x in values):
                    return ConflictSeverity.MINOR # Is VUS_A vs VUS
                else:
                    return ConflictSeverity.MEDIUM
            else:
                return ConflictSeverity.SAME


class ClinSigCalculator(ConflictCalculator):

    def should_exclude(self, row: ConflictDataRow):
        return super().should_exclude(row) or row.clinical_significance is None

    def calculate_severity(self) -> ConflictSeverity:
        clin_sig_data = self.included_conflict_data
        if len(clin_sig_data) == 0:
            return ConflictSeverity.NO_SUBMISSIONS
        elif len(clin_sig_data) == 1:
            return ConflictSeverity.SINGLE_SUBMISSION
        else:
            has_tier_1_and_2 = False
            tiers = set()
            for row in self.conflict_data:
                # TODO this is hardcoding values expected in evidence keys
                cs = row.clinical_significance
                if cs == "tier_1_or_2":
                    has_tier_1_and_2 = True
                else:
                    tiers.add(cs)

            # for now ever difference in tier is a MAJOR diff
            # the only exception being tier_1 vs tier_1_or_2 or tier_2 vs tier_1_or_2 which
            # is considered minor
            if len(tiers) > 1:
                return ConflictSeverity.MEDICALLY_SIGNIFICANT
            elif len(tiers) == 1 and has_tier_1_and_2:
                only_tier = first(tiers)
                if only_tier in {"tier_1", "tier_2"}:
                    return ConflictSeverity.MINOR # tier 1 vs tier 1_or_2
                else:
                    return ConflictSeverity.MEDICALLY_SIGNIFICANT
            return ConflictSeverity.SAME


def calculate_and_apply_conflicts_for(allele_origin_grouping: AlleleOriginGrouping):
    oncpath_data: list[ConflictDataRow] = []
    clinsig_data: list[ConflictDataRow] = []
    allele_origin_bucket = allele_origin_grouping.allele_origin_bucket

    for cg in allele_origin_grouping.classificationgrouping_set.select_related("latest_classification_modification__classification"):
        lab_id = cg.latest_classification_modification.classification.lab_id
        summary_info = cg.latest_classification_modification.classification.summary_typed
        c_hgvs = cg.latest_classification_modification.imported_c_hgvs_obj

        oncpath_data.append(ConflictDataRow.from_data(row=summary_info.get("pathogenicity", {}), lab_id=lab_id, share_level=cg.share_level_obj, c_hgvs=c_hgvs))
        clinsig_data.append(ConflictDataRow.from_data(row=summary_info.get("somatic", {}), lab_id=lab_id, share_level=cg.share_level_obj, c_hgvs=c_hgvs))

    allele = allele_origin_grouping.allele_grouping.allele
    if allele_origin_bucket == AlleleOriginBucket.SOMATIC:
        conflict_key = ConflictKey(
            conflict_type=ConflictType.CLIN_SIG,
            allele_origin_bucket=allele_origin_bucket,
            testing_context_bucket=allele_origin_grouping.testing_context_bucket,
            tumor_type_category=allele_origin_grouping.tumor_type_category
        )
        ClinSigCalculator(allele, conflict_key, clinsig_data).log_history()

    conflict_key = ConflictKey(
        conflict_type=ConflictType.ONCPATH,
        allele_origin_bucket=allele_origin_bucket,
        testing_context_bucket=allele_origin_grouping.testing_context_bucket,
        tumor_type_category=allele_origin_grouping.tumor_type_category
    )
    OncPathCalculator(allele, conflict_key, oncpath_data).log_history()


def latest_classifications_for_conflict(conflict: Conflict) -> Optional[ClassificationModification]:
    # easy case before we get into
    if grouping := AlleleOriginGrouping.objects.filter(
        allele_grouping__allele=conflict.allele,
        allele_origin_bucket=conflict.allele_origin_bucket,
        testing_context_bucket=conflict.testing_context_bucket,
        tumor_type_category=conflict.tumor_type_category
    ).first():
        # TODO filter out records that don't have a relevant value? Still probably pretty useful to include them
        return [cg.latest_classification_modification for cg in ClassificationGrouping.objects.filter(allele_origin_grouping=grouping)]
    else:
        raise ValueError(f"Can't find AlleleOriginGrouping for conflict {conflict}")


@dataclass(frozen=True)
class ConflictMerge:
    conflicts: list[Conflict]

    @cached_property
    def conflict(self) -> Conflict:
        return first(self.conflicts)

    @cached_property
    def conflict_types(self) -> list[ConflictType]:
        conflict_types = set()
        for conflict in self.conflicts:
            conflict_types.add(ConflictType(conflict.conflict_type))
        return list(sorted(conflict_types))

    @cached_property
    def conflict_type_labels(self) -> list[str]:
        # conflict type label can change depending on allele origin (even if calculated the same)
        return [conflict_type.label_for_context(self.conflict.allele_origin_bucket) for conflict_type in self.conflict_types]

    @staticmethod
    def merge(conflicts: Iterable[Conflict]) -> list['ConflictMerge']:
        # groupby conflict key but with testing context hardcoded to None, so records with the same everything else will be grouped together
        all_merged: list[ConflictMerge] = []
        for group, conflicts in itertools.groupby(sorted(conflicts), lambda c: ConflictKey(
                None,
                c.allele_origin_bucket,
                c.testing_context_bucket,
                c.tumor_type_category,
                c.latest.severity
            )):
            conflict_merge = ConflictMerge(list(conflicts))
            all_merged.append(conflict_merge)
        return all_merged

    def __getattr__(self, method_name):
        # delegate to first conflict object if we get a method we're not familiar with
        def method(*args, **kwargs):
            result = getattr(self.conflict, method_name)
            if isinstance(result, Callable):
                return result(*args, **kwargs)
            else:
                return result

        return method


class GroupConflictsLinks(IntEnum):
    NO_LINKS = 0
    OVERLAP_LINK = 1
    FILTER_LINK = 2


@dataclass
class ConflictSeverityIconAndText:
    html: str
    css_classes: list[str]
    unique_classifications: set[str]
    unique_clin_sigs: set[str]


def conflict_severity_icon_and_text(user_lab_ids: Collection[int], severity: ConflictSeverity, data_rows: list[ConflictDataRow]) -> ConflictSeverityIconAndText:
    # severity_type_css = [allele_origin_css]
    # if currently_viewing == conflict:
    #     severity_type_css += ["strong"]

    severity_type_css = []
    if severity < ConflictSeverity.SAME:
        severity_type_css += ["no-value"]
    elif severity == ConflictSeverity.SAME:
        pass
    elif severity < ConflictSeverity.MAJOR:
        severity_type_css += ["strong", "text-warning"]
    else:
        severity_type_css += ["strong", "text-danger"]

    unique_classifications = set()
    unique_clin_sigs = set()
    your_lab_ids = set()
    other_lab_ids = set()
    for row in data_rows:
        if not row.exclude:
            if user_lab_ids is not None:
                if row.lab_id in user_lab_ids:
                    your_lab_ids.add(row.lab_id)
                else:
                    other_lab_ids.add(row.lab_id)

            if classification := row.classification:
                unique_classifications.add(classification)
            if clinical_sig := row.clinical_significance:
                unique_clin_sigs.add(clinical_sig)

    # now check the labs
    severity_parts = []
    if user_lab_ids is not None:
        your_lab_count = len(your_lab_ids)
        other_lab_count = len(other_lab_ids)

        if your_lab_count:
            severity_parts.append('<i class="fa-solid fa-user-doctor text-success" title="Your lab is involved"></i>')
        if other_lab_count == 1:
            severity_parts.append('<i class="fa-solid fa-user text-secondary" title="Another lab is involved"></i>')
        elif other_lab_count > 1:
            severity_parts.append(
                '<i class="fa-solid fa-users text-secondary" title="Multiple other labs are involved"></i>')

    severity_parts.append(escape(severity.label))
    severity_str = SafeString("".join(severity_parts))
    return ConflictSeverityIconAndText(
        html=severity_str,
        css_classes=severity_type_css,
        unique_classifications=unique_classifications,
        unique_clin_sigs=unique_clin_sigs,
    )


def group_conflicts(
        conflicts: Iterable[Conflict],
        link_types: GroupConflictsLinks = GroupConflictsLinks.NO_LINKS,
        currently_viewing: Optional[Conflict] = None,
        user: Optional[User] = None,
    ) -> 'RowSpanTable':

    def href_for_parts(parts: list[str]):
        nonlocal link_types
        if link_types == GroupConflictsLinks.FILTER_LINK:
            parts_array = ",".join(f"'{p}'" for p in parts)
            return f"overlap_filter([{ parts_array }])"
        return None

    user_lab_ids = None
    if user:
        user_lab_ids = set(lab.pk for lab in Lab.valid_labs_qs(user, admin_check=False))

    conflicts = sorted(conflicts, key=lambda c: c.conflict_type)

    table = RowSpanTable(6)
    classification_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
    clin_sig_key = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)

    for allele_origin, sub_conflicts in sort_and_group_by(conflicts, lambda c: c.allele_origin_bucket):
        allele_origin_bucket = AlleleOriginBucket(allele_origin)
        link_parts = [allele_origin_bucket.value, allele_origin_bucket.label]

        allele_origin_css = f"allele-origin-{allele_origin_bucket.value}"
        allele_origin_full_css = [allele_origin_css]
        allele_origin_matches = False
        if currently_viewing and currently_viewing.allele_origin_bucket == allele_origin_bucket:
            allele_origin_matches = True
            allele_origin_full_css += ["currently-viewing"]

        table.add_cell(0, RowSpanCellValue(
            allele_origin_bucket.label,
            css_classes=allele_origin_full_css,
            href=href_for_parts(link_parts),
            link_type=LinkType.JAVASCRIPT
       ))

        for testing_context_str, sub_sub_conflicts in sort_and_group_by(sub_conflicts, lambda c: c.testing_context_bucket):

            testing_context = TestingContextBucket(testing_context_str)
            testing_link_parts = link_parts + [testing_context.value, testing_context.label]

            testing_context_css = [allele_origin_css]
            testing_context_matches = False
            if allele_origin_matches and currently_viewing and currently_viewing.testing_context_bucket == testing_context:
                testing_context_matches = True
                testing_context_css += ["currently-viewing"]

            href: Optional[str] = None
            if allele_origin_bucket != AlleleOriginBucket.GERMLINE:
                # TODO what does germline or not have to do with it??
                href = href_for_parts(testing_link_parts)

            table.add_cell(1, RowSpanCellValue(testing_context.label, css_classes=testing_context_css, href=href, link_type=LinkType.JAVASCRIPT))

            for condition, sub_sub_sub_conflicts in sort_and_group_by(sub_sub_conflicts, lambda c: c.tumor_type_category or "PLACEHOLDER"):
                condition_css = [allele_origin_css]
                condition_matches = False
                if testing_context_matches and currently_viewing and (currently_viewing.tumor_type_category or "PLACEHOLDER") == condition:
                    condition_matches = True
                    condition_css += ["currently-viewing"]

                if condition == "PLACEHOLDER":
                    condition_css += ["no-value"]
                    condition = "N/A"
                    if testing_context.should_have_subdivide:
                        condition = "Indeterminate condition"

                table.add_cell(2, RowSpanCellValue(condition, css_classes=condition_css))

                for value_type, sub_sub_sub_sub_conflicts in sort_and_group_by(sub_sub_sub_conflicts, lambda c: ConflictType(c.conflict_type)):
                    sub_sub_sub_sub_conflicts_list = list(sub_sub_sub_sub_conflicts)
                    if len(sub_sub_sub_sub_conflicts_list) != 1:
                        raise ValueError(f"Multiple conflicts found for {allele_origin} {testing_context} {condition} {sub_sub_sub_sub_conflicts_list}")

                    conflict = sub_sub_sub_sub_conflicts_list[0]
                    value_type_label = ConflictType(value_type).label_for_context(allele_origin_bucket)
                    value_type_css = [allele_origin_css]
                    if currently_viewing:
                        if condition_matches and currently_viewing and (currently_viewing.conflict_type == value_type):
                            value_type_css += ["currently-viewing"]
                        else:
                            value_type_css += ["not-viewing"]

                    value_type_css += ["value-type"]
                    href = None
                    link_type = LinkType.JAVASCRIPT
                    if link_types == GroupConflictsLinks.NO_LINKS:
                        pass
                    # the below provides modal links on the allele page to a conflict history - but it makes the experience
                    # a little too cluttered
                    # elif not currently_viewing:
                    #     href = reverse("conflict_feed", kwargs={"conflict_id": conflict.pk})
                    #     link_type = LinkType.MODAL
                    elif currently_viewing != conflict:
                        href = reverse("conflict", kwargs={"conflict_id": conflict.pk})
                        link_type = LinkType.LINK

                    table.add_cell(3, RowSpanCellValue(value_type_label, css_classes=value_type_css, href=href, link_type=link_type))

                    # severity = ConflictSeverity(conflict.latest.severity)
                    # severity_type_css = [allele_origin_css]
                    # if severity < ConflictSeverity.SAME:
                    #     severity_type_css += ["no-value"]
                    # elif severity == ConflictSeverity.SAME:
                    #     pass
                    # elif severity < ConflictSeverity.MAJOR:
                    #     severity_type_css += ["strong", "text-warning"]
                    # else:
                    #     severity_type_css += ["strong", "text-danger"]
                    #
                    # if currently_viewing == conflict:
                    #     severity_type_css += ["strong"]
                    #
                    # conflict_history: ConflictHistory = conflict.latest
                    # unique_classifications = set()
                    # unique_clin_sigs = set()
                    # your_lab_ids = set()
                    # other_lab_ids = set()
                    # for row in conflict_history.data_rows():
                    #     if not row.exclude:
                    #         if user_lab_ids is not None:
                    #             if row.lab_id in user_lab_ids:
                    #                 your_lab_ids.add(row.lab_id)
                    #             else:
                    #                 other_lab_ids.add(row.lab_id)
                    #
                    #         if classification := row.classification:
                    #             unique_classifications.add(classification)
                    #         if clinical_sig := row.clinical_significance:
                    #             unique_clin_sigs.add(clinical_sig)
                    #
                    #
                    # # now check the labs
                    #
                    # severity_parts = []
                    # if user_lab_ids is not None:
                    #     your_lab_count = len(your_lab_ids)
                    #     other_lab_count = len(other_lab_ids)
                    #
                    #     if your_lab_count:
                    #         severity_parts.append('<i class="fa-solid fa-user-doctor text-success" title="Your lab is involved"></i>')
                    #     if other_lab_count == 1:
                    #         severity_parts.append('<i class="fa-solid fa-user text-secondary" title="Another lab is involved"></i>')
                    #     elif other_lab_count > 1:
                    #         severity_parts.append('<i class="fa-solid fa-users text-secondary" title="Multiple other labs are involved"></i>')
                    #
                    # severity_parts.append(escape(severity.label))
                    # severity_str = SafeString("".join(severity_parts))

                    parts = conflict_severity_icon_and_text(user_lab_ids, severity=ConflictSeverity(conflict.latest.severity), data_rows=conflict.latest.data_rows())
                    severity_type_css = [] + parts.css_classes
                    severity_type_css += [allele_origin_css]
                    if currently_viewing == conflict:
                        severity_type_css += ["strong"]
                    table.add_cell(4, RowSpanCellValue(parts.html, css_classes=severity_type_css))

                    values = []
                    value_css = [allele_origin_css]
                    if currently_viewing == conflict:
                        value_css += ["strong"]
                    if unique_classifications := parts.unique_classifications:
                        values.extend(classification_key.sort_values(unique_classifications))

                    if unique_clin_sigs := parts.unique_clin_sigs:
                        values.extend(clin_sig_key.pretty_value(val) for val in clin_sig_key.sort_values(unique_clin_sigs))

                    values_str: str
                    if values:
                        values_str = ", ".join(values)
                        table.add_cell(5, RowSpanCellValue(values_str, css_classes=value_css))
                    else:
                        table.add_cell(5, RowSpanCellValue("-", css_classes=value_css + ["no-value"]))

    return table


def conflict_lab_for_grouping(classification_grouping: ClassificationGrouping, conflict_type: ConflictType) -> Optional[ConflictLab]:
    # TODO could do this in a single filter but easier to debug this way
    if conflict := Conflict.objects.filter(
            allele=classification_grouping.allele,
            conflict_type=conflict_type,
            allele_origin_bucket=classification_grouping.allele_origin_bucket,
            testing_context_bucket=classification_grouping.allele_origin_grouping.testing_context_bucket,
            tumor_type_category=classification_grouping.allele_origin_grouping.tumor_type_category,
        ).first():
        if conflict_lab := ConflictLab.objects.filter(
            conflict=conflict,
            lab=classification_grouping.lab
        ).first():
            return conflict_lab
    return None


def classification_grouping_for_conflict(conflict: Conflict, lab: Lab) -> Optional[ClassificationGrouping]:
    # TODO could do this in a single filter but easier to debug this way
    if ag := AlleleGrouping.objects.filter(allele=conflict.allele).first():
        if aog := AlleleOriginGrouping.objects.filter(
            allele_grouping=ag,
            allele_origin_bucket=conflict.allele_origin_bucket,
            testing_context_bucket=conflict.testing_context_bucket,
            tumor_type_category=conflict.tumor_type_category
        ).first():
            if cg := ClassificationGrouping.objects.filter(allele_origin_grouping=aog, lab=lab).first():
                return cg
    print(f"Couldn't find classification grouping for {conflict} {lab}")
    return None


# def apply_conflict_lab_to_grouping(conflict_lab: ConflictLab):
#     if classification_grouping := classification_grouping_for_conflict_lab(conflict_lab):
#         will_amend = conflict_lab.status == DiscordanceReportTriageStatus.REVIEWED_WILL_FIX
#         match conflict_lab.conflict.conflict_type:
#             case ConflictType.ONCPATH:
#                 if classification_grouping.pending_change_onc_path != will_amend:
#                     classification_grouping.pending_change_onc_path = will_amend
#                     classification_grouping.save()
#             case ConflictType.CLIN_SIG:
#                 if classification_grouping.pending_change_clin_sig != will_amend:
#                     classification_grouping.pending_change_clin_sig = will_amend
#                     classification_grouping.save()
#             case _:
#                 pass
