import re
from dataclasses import dataclass
from typing import Optional

from django.http import HttpRequest
from django.urls import reverse

from classification.enums import SpecialEKeys
from classification.models import ClassificationModification, EvidenceKeyMap
from classification.views.classification_export_view import InvalidExportParameter
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, delimited_row, ExportTweak
from snpdb.models import Lab, Allele


class ClassificationLab(ExportRow):

    def __init__(self, lab: Lab, cms: list[ClassificationModification]):
        self.lab = lab
        self.cms = cms

    def list_of(self, *args: list[str]) -> Optional[str]:
        """
        Return a single string that's the combination of all the evidence keys found in args
        :param args: strings that match EvidenceKey keys
        :return: "\n" seperated unique combined values for the given evidence keys
        """
        value_set = set()
        for cm in self.cms:
            value_tuple = [cm.get(e_key) for e_key in args]
            value_tuple = tuple(str(x) if x else None for x in value_tuple)
            if any(value_tuple):
                value_set.add(value_tuple)
        if value_set:
            value_set_strs = []
            for vs in sorted(value_set):
                value_set_strs.append(" ".join(v or "#" for v in vs))
            return "\n".join(value_set_strs)

    @export_column("Record Count")
    def record_count(self):
        return len(self.cms)

    @export_column("Lab Record ID")
    def lab_record_id(self):
        if len(self.cms) <= 5:
            return ", ".join(cms.classification.lab_record_id for cms in self.cms)
        else:
            return "multiple"

    @export_column("c.HGVSs")
    def c_hgvs(self):
        return self.list_of(SpecialEKeys.C_HGVS)

    @export_column("Classification")
    def classification_value(self):
        return self.list_of(SpecialEKeys.CLINICAL_SIGNIFICANCE)

    @export_column("Patient IDs")
    def patient_ids(self):
        return self.list_of(SpecialEKeys.PATIENT_ID)

    @export_column("Curation Details")
    def curation_details(self):
        return self.list_of(SpecialEKeys.CURATION_DATE, SpecialEKeys.CURATED_BY)

    @export_column("Curation Verified Details")
    def curation_verified_details(self):
        return self.list_of(SpecialEKeys.CURATION_VERIFIED_DATE, SpecialEKeys.CURATION_VERIFIED_BY)

    @export_column("Interpretation Summary")
    def interpretation_summary(self):
        return self.list_of(SpecialEKeys.INTERPRETATION_SUMMARY)


@dataclass
class ClassificationLabCompare(ExportRow):

    allele: Allele
    lab_a_data: ClassificationLab
    lab_b_data: ClassificationLab
    e_keys: EvidenceKeyMap

    @classmethod
    def zip_sub_data(cls) -> set[type]:
        # zip sub data means Lab A and Lab B data will be combined rather than 1 right after the other
        # e.g. Lab A.Record Count, Lab B.Record Count, Lab A.c.HGVS, Lab B.c.HGVS etc
        return {ClassificationLab}

    @export_column("Allele")
    def _allele(self):
        if self.allele:
            return get_url_from_view_path(reverse('view_allele', kwargs={"allele_id": self.allele.pk}))

    @export_column("Lab A", sub_data=ClassificationLab)
    def _lab_a(self):
        return self.lab_a_data

    @export_column("Lab B", sub_data=ClassificationLab)
    def _lab_b(self):
        return self.lab_b_data

    @staticmethod
    def tokenize(data: Optional[str]):
        if data is None:
            return set()
        return set(re.findall(r'\b\w+\b', data.lower()))

    @export_column("Fields with Differences")
    def differences(self):
        if not self.lab_a_data.cms:
            return "ONLY " + self.lab_b_data.lab.name
        elif not self.lab_b_data.cms:
            return "ONLY " + self.lab_a_data.lab.name

        different_keys = list()
        for key in self.e_keys.all_keys:
            a_values = ClassificationLabCompare.tokenize(self.lab_a_data.list_of(key.key))
            b_values = ClassificationLabCompare.tokenize(self.lab_b_data.list_of(key.key))
            if a_values != b_values:
                different_keys.append(key.pretty_label)

        if different_keys:
            return ", ".join(different_keys)
        else:
            return "SAME"


@register_classification_exporter("lab_compare")
class ClassificationExportInternalCompare(ClassificationExportFormatter):
    """
    Formatter typically used to compare old lab data to new lab data (after the labs have been inserted into test as
    different labs). Identifies what has accidentally changed
    """

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportInternalCompare2':
        return ClassificationExportInternalCompare(
            classification_filter=ClassificationFilter.from_request(request),
        )

    def __init__(self, classification_filter: ClassificationFilter):
        super().__init__(classification_filter)
        if not self.classification_filter.include_sources or len(self.classification_filter.include_sources) != 2:
            raise InvalidExportParameter("Lab Compare requires the <b>inclusion</b> of 2 labs.")
        two_labs = list(sorted(self.classification_filter.include_sources))
        self.lab_a = two_labs[0]
        self.lab_b = two_labs[1]
        self.e_keys = EvidenceKeyMap.cached()

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def lab_label(self, lab_key: str, field_name: str):
        if lab_key == "Lab A":
            lab_key = self.lab_a.name
        elif lab_key == "Lab B":
            lab_key = self.lab_b.name
        return f"{field_name} ({lab_key})"

    def header(self) -> list[str]:
        return [
            delimited_row(
                ClassificationLabCompare.csv_header(export_tweak=ExportTweak(sub_label_joiner=self.lab_label))
            )
        ]

    def row(self, allele_data: AlleleData) -> list[str]:
        if not allele_data.allele:
            return []

        lab_a_cs: list[ClassificationModification] = []
        lab_b_cs: list[ClassificationModification] = []
        for cm in allele_data.cms_regardless_of_issues:
            if cm.lab == self.lab_a:
                lab_a_cs.append(cm)
            elif cm.lab == self.lab_b:
                lab_b_cs.append(cm)

        if not lab_a_cs and not lab_b_cs:
            return []

        lab_a_data = ClassificationLab(self.lab_a, lab_a_cs)
        lab_b_data = ClassificationLab(self.lab_b, lab_b_cs)

        return [delimited_row(
            ClassificationLabCompare(allele_data.allele, lab_a_data, lab_b_data, self.e_keys).to_csv())
        ]
