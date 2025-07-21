from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, TypeVar, Generic

from classification.enums import AlleleOriginBucket, ConflictSeverity, ShareLevel, ConflictType
from classification.models import AlleleOriginGrouping, EvidenceKeyMap, Conflict, ConflictKey

"""
    classification: Optional[str]
    sort: int
    bucket: Optional[int]
    pending: Optional[str]
"""


@dataclass(frozen=True)
class ConflictDataRow(ABC):
    lab_id: int
    share_level: ShareLevel

    @property
    @abstractmethod
    def to_json_value(self) -> dict:
        raise NotImplementedError()


T = TypeVar("T", bound=ConflictDataRow)


class ConflictCalculator(ABC, Generic[T]):

    @abstractmethod
    def calculate_severity(self) -> ConflictSeverity:
        raise NotImplementedError()

    @property
    @abstractmethod
    def data_list(self) -> list[T]:
        raise NotImplementedError()

    @staticmethod
    def remove_duplicate_data(data_list: list[T]):
        by_lab_id = {}
        for row in data_list:
            if existing := by_lab_id.get(row.lab_id):
                if existing.share_level > row.share_level:
                    by_lab_id[row.lab_id] = row
                elif existing.share_level == row.share_level:
                    raise ValueError(f"Lab ID {row.lab_id} {existing.share_level} found multiple times within a context")
            else:
                by_lab_id[row.lab_id] = row
        return list(by_lab_id.values())

    def data_json(self) -> dict:
        data_json = {}
        for row in self.data_list:
            data_json[row.lab_id] = row.to_json_value()
        return data_json


@dataclass(frozen=True)
class OncPathData(ConflictDataRow):
    classification: str

    @staticmethod
    def from_data(lab_id: int, share_level: ShareLevel, row: dict) -> 'OncPathData':
        return OncPathData(lab_id=lab_id, share_level=share_level, classification=row.get("classification"))

    @property
    def bucket(self) -> Optional[int]:
        return EvidenceKeyMap.clinical_significance_to_bucket().get(self.classification)

    def to_json_value(self) -> dict:
        return {
            "classification": self.classification
        }


class OncPathCalculator(ConflictCalculator):

    def __init__(self, conflict_key: ConflictKey, oncpath_data: list[OncPathData]):
        self.conflict_key = conflict_key
        self.oncpath_data = ConflictCalculator.remove_duplicate_data(oncpath_data)

    def calculate_severity(self) -> ConflictSeverity:
        oncpath_data = self.oncpath_data
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
                if 3 in buckets:  # don't like this hardcoding, 3 = LP, P
                    return ConflictSeverity.MAJOR
                else:
                    return ConflictSeverity.MEDIUM
            if len(values) > 1:
                return ConflictSeverity.MINOR
            else:
                return ConflictSeverity.SAME

    @property
    def data_list(self) -> list[T]:
        return self.oncpath_data


@dataclass(frozen=True)
class ClinSigData(ConflictDataRow):
    clin_sig: str
    amp_level: int

    @staticmethod
    def from_data(lab_id: int, share_level: ShareLevel, row: dict) -> 'ClinSigData':
        return ClinSigData(lab_id=lab_id, share_level=share_level, clin_sig=row.get("clinical_significance"), amp_level=row.get("amp_level"))

    def to_json_value(self) -> dict:
        return {
            "clinical_significance": self.clin_sig,
            "amp_level": self.amp_level
        }


class ClinSigCalculator(ConflictCalculator):

    def __init__(self, conflicy_key: ConflictKey, clinsig_data: list[ClinSigData]):
        self.conflicy_key = conflicy_key
        self.clinsig_data = ConflictCalculator.remove_duplicate_data(clinsig_data)

    def calculate_severity(self) -> ConflictSeverity:
        if len(self.clinsig_data) == 0:
            return ConflictSeverity.NO_SUBMISSIONS
        elif len(self.clinsig_data) == 1:
            return ConflictSeverity.SINGLE_SUBMISSION
        else:
            return ConflictSeverity.MEDIUM

    @property
    def data_list(self) -> list[T]:
        return self.clinsig_data


def calculate_and_apply_conflicts_for(allele_origin_grouping: AlleleOriginGrouping):
    oncpath_data: list[OncPathData] = []
    clinsig_data: list[ClinSigData] = []
    allele_origin_bucket = allele_origin_grouping.allele_origin_bucket

    for cg in allele_origin_grouping.classificationgrouping_set.select_related("latest_classification_modification__classification").filter(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).all():
        lab_id = cg.latest_classification_modification.classification.lab_id
        summary_info = cg.latest_classification_modification.classification.summary_typed

        oncpath_data.append(OncPathData.from_data(lab_id=lab_id, share_level=cg.share_level, row=summary_info.get("pathogenicity", {})))
        clinsig_data.append(ClinSigData.from_data(lab_id=lab_id, share_level=cg.share_level, row=summary_info.get("clinical_significance", {})))

    try:
        if allele_origin_bucket == AlleleOriginBucket.SOMATIC:
            conflict_key = ConflictKey(
                conflict_type=ConflictType.CLIN_SIG,
                allele_origin_bucket=allele_origin_bucket,
                testing_context_bucket=allele_origin_grouping.testing_context_bucket,
                tumor_type_category=allele_origin_grouping.tumor_type_category
            )
            clin_sig_calc = ClinSigCalculator(conflict_key, clinsig_data)
            Conflict.log_history(allele=allele_origin_grouping.allele_grouping.allele, conflict_key=conflict_key, data=clin_sig_calc.data_json(), severity=clin_sig_calc.calculate_severity())

        conflict_key = ConflictKey(
            conflict_type=ConflictType.ONCPATH,
            allele_origin_bucket=allele_origin_bucket,
            testing_context_bucket=allele_origin_grouping.testing_context_bucket,
            tumor_type_category=allele_origin_grouping.tumor_type_category
        )
        onc_path_calc = OncPathCalculator(conflict_key, oncpath_data)
        Conflict.log_history(allele=allele_origin_grouping.allele_grouping.allele, conflict_key=conflict_key,
                             data=onc_path_calc.data_json(), severity=onc_path_calc.calculate_severity())
    except ValueError as ve:
        print(oncpath_data)
        print(clinsig_data)
        raise ValueError(f"Multiple entries for lab found under {allele_origin_grouping.allele_grouping.allele} - {conflict_key} ({ve})")
