from typing import List

from django.http import HttpRequest
from django.urls.base import reverse

from classification.enums import SpecialEKeys
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, delimited_row


class ClassificationlabCompareRow(ExportRow):

    def __init__(self,  AlleleData, lab1, lab2):
        self.AlleleData = AlleleData
        self.lab1 = lab1
        self.lab2 = lab2

    @export_column("Allele URL")
    def allele_url(self):
        url = get_url_from_view_path(reverse('view_allele', kwargs={"allele_id": self.AlleleData}))
        return url

    @export_column("Lab 1")
    def lab_clinical_significance(self):
        return self.lab1

    @export_column("Lab 2")
    def lab2_clinical_significance(self):
        return self.lab2

    @export_column("Difference")
    def difference(self, lab_name1, lab_name2):
        if (self.lab1 != '' and self.lab2 != '') and self.lab1 == self.lab2:
            return 'Same'
        elif self.lab1 == '' and self.lab2 != '':
            return lab_name2
        elif self.lab1 != '' and self.lab2 == '':
            return lab_name1
        elif self.lab1 != self.lab2:
            return 'Different'
        elif self.lab1 == '' and self.lab2 == '':
            return 'No data'


@register_classification_exporter("lab_compare")
class ClassificationExportInternalCompare(ClassificationExportFormatter):

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportInternalCompare':
        return ClassificationExportInternalCompare(
            classification_filter=ClassificationFilter.from_request(request),
        )

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> List[str]:
        if self.classification_filter.include_sources and len(self.classification_filter.include_sources) == 2:
            lab_names = sorted([str(lab) for lab in self.classification_filter.include_sources])
            return [delimited_row(['Allele URL', lab_names[0], lab_names[1], 'Difference'], ',')]
        else:
            raise ValueError("Must specify 2 labs to compare")

    def footer(self) -> List[str]:
        return []

    def row(self, allele_data: AlleleData) -> List[str]:
        rows: List[str] = []
        lab1 = set()
        lab2 = set()

        if not allele_data.allele_id:
            return []
        if self.classification_filter.include_sources:
            lab_name1, lab_name2 = sorted(self.classification_filter.include_sources)
            for cm in allele_data.cms_regardless_of_issues:
                if cm.lab == lab_name1 and cm.clinical_significance is not None:
                    lab1.add(cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
                if cm.lab == lab_name2 and cm.clinical_significance is not None:
                    lab2.add(cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))

            lab1_set = ','.join(lab1)
            lab2_set = ','.join(lab2)
            row = ClassificationlabCompareRow(AlleleData=allele_data.allele_id, lab1=lab1_set, lab2=lab2_set)
            rows.append(delimited_row([row.allele_url(), row.lab_clinical_significance(),
                                       row.lab2_clinical_significance(), row.difference(lab_name1, lab_name2)], ','))
        else:
            raise ValueError("Error comparing labs")
        return rows
