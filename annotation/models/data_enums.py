from dataclasses import dataclass
from datetime import datetime, date
from enum import StrEnum
from typing import Optional, Union
from dataclasses_json import DataClassJsonMixin
from django.db.models import TextChoices


class ClassificationDateType(StrEnum):
    # To be deprecated for more general EffectiveDateType

    CURATION = ""  # default
    VERIFIED = "Verified"
    SAMPLE_DATE = "Sample Date"
    CREATED = "Created"


class EffectiveDateType(TextChoices):
    CREATED = "created", "Created"
    CURATED = "curated", "Curated"
    SAMPLE_DATE = "sample", "Sample"
    VERIFIED = "verified", "Verified"
    UNKNOWN = "unknown", "Unknown"

    @staticmethod
    def from_classification_date_type(classification_date_type: 'ClassificationDateType') -> 'EffectiveDateType':
        match classification_date_type:
            case ClassificationDateType.CREATED: return EffectiveDateType.CREATED
            case ClassificationDateType.CURATION: return EffectiveDateType.CURATED
            case ClassificationDateType.SAMPLE_DATE: return EffectiveDateType.SAMPLE_DATE
            case ClassificationDateType.VERIFIED: return EffectiveDateType.VERIFIED
            case _: return EffectiveDateType.UNKNOWN


@dataclass
class EffectiveDate(DataClassJsonMixin):
    date: Optional[str] = None
    date_type: EffectiveDateType = EffectiveDateType.UNKNOWN

    @staticmethod
    def from_datetime(value: Union[datetime, date], date_type: EffectiveDateType = EffectiveDateType.UNKNOWN):
        date_str: Optional[str] = None
        value_date: Optional[date]
        if isinstance(value, date):
            value_date = value
        elif isinstance(value, datetime):
            value_date = datetime.date()
        else:
            raise ValueError(f"Not datetime or date {value}")
        if value_date:
            date_str = f"{value_date.year:04}-{value_date.month:02}-{value_date.day:02}"
        return EffectiveDate(date=date_str, date_type=date_type)

    @staticmethod
    def default_json():
        return EffectiveDate().to_dict()

    def __str__(self):
        if date_val := self.date:
            return f"{date_val} {self.date_type.label}"
        return "Date Unknown"


@dataclass(frozen=True)
class ClassificationDate:
    # Deprecated, migrate to EffectiveDate

    date_type: ClassificationDateType
    datetime: Optional[datetime] = None
    date: Optional[date] = None

    def __post_init__(self):
        if not self.date and not self.datetime:
            raise ValueError("Either datetime or date must be provided")

    @property
    def name(self):
        return self.date_type.value

    @property
    def value(self) -> Union[date, datetime]:
        if self.date:
            return self.date
        return self.datetime

    @property
    def date_str(self) -> str:
        from classification.models import Classification
        return Classification.to_date_str(self.value)

    def __lt__(self, other: 'ClassificationDate'):
        if self.datetime and other.datetime:
            return self.datetime < other.datetime
        else:
            return (self.date or date.min) < (other.date or date.min)
