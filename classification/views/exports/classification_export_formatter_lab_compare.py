import re
from collections import defaultdict
from typing import List, Optional

from django.http import HttpRequest
from django.urls.base import reverse

from classification.enums import SpecialEKeys
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, delimited_row


class ClassificationlabCompareRow(ExportRow):

    def __init__(self,  allele_data: AlleleData, lab_1: Optional[str] = None, lab_2: Optional[str] = None,
                 comment: Optional[str] = None, lab_1_authorised_date: Optional[str] = None,
                 lab_2_authorised_date: Optional[str] = None, lab_1_curated_date: Optional[str] = None,
                 lab_2_curated_date: Optional[str] = None, lab_1_interpretation_summary: Optional[str] = None,
                 lab_2_interpretation_summary: Optional[str] = None,
                 patient_id: Optional[str] = None):
        self.AlleleData = allele_data.allele_id
        self.lab1 = lab_1
        self.lab2 = lab_2
        self.comment = comment
        self.lab_1_authorised_date = lab_1_authorised_date
        self.lab_2_authorised_date = lab_2_authorised_date
        self.lab_1_curated_date = lab_1_curated_date
        self.lab_2_curated_date = lab_2_curated_date
        self.Patient_id = patient_id
        self.lab_1_interpretation_summary = lab_1_interpretation_summary
        self.lab_2_interpretation_summary = lab_2_interpretation_summary

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

    @export_column("Classification")
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

    @export_column("Authorised Date Difference")
    def authorised_date_diff(self):
        if not self.lab_1_authorised_date or not self.lab_2_authorised_date:
            return 'No Overlap'
        elif self.lab_1_authorised_date != self.lab_2_authorised_date:
            return 'Different'
        elif self.lab_1_authorised_date == self.lab_2_authorised_date:
            return 'Same'
        else:
            return ''

    @export_column("Created Date Difference")
    def created_date_diff(self):
        if not self.lab_1_curated_date or not self.lab_2_curated_date:
            return 'No Overlap'
        elif self.lab_1_curated_date != self.lab_2_curated_date:
            return 'Different'
        elif self.lab_1_curated_date == self.lab_2_curated_date:
            return 'Same'
        else:
            return ''

    @export_column("Fields With Differences")
    def comment(self):
        return ''

    @export_column("Patient_Id")
    def patient_id(self):
        return self.Patient_id


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
            return [delimited_row(['Allele URL', 'Patient_Id', lab_names[0], lab_names[1],
                                   'Classification',
                                   f'{lab_names[0]} Authorised Date',
                                   f'{lab_names[1]} Authorised Date',
                                   'Authorised Date Difference',
                                   f'{lab_names[0]} Curated Date',
                                   f'{lab_names[1]} Curated Date',
                                   'Curated Date Difference',
                                   'Fields With Differences',
                                   f'{lab_names[0]} Interpretation Summary',
                                   f'{lab_names[1]} Interpretation Summary'
                                   ], ',')]
        else:
            raise ValueError("Must specify 2 labs to compare")

    def footer(self) -> List[str]:
        return []

    def row(self, allele_data: AlleleData) -> List[str]:

        lab1_comment_data = defaultdict(list)
        lab2_comment_data = defaultdict(list)

        rows: List[str] = []
        lab1 = set()
        lab2 = set()
        message = ''
        patient_ids = set()
        authorised_date_lab1 = set()
        authorised_date_lab2 = set()
        created_date_lab1 = set()
        created_date_lab2 = set()
        interpretation_summary_lab1 = set()
        interpretation_summary_lab2 = set()

        if not allele_data.allele_id:
            return []
        if self.classification_filter.include_sources:
            lab_name1, lab_name2 = sorted(self.classification_filter.include_sources)
            for cm in allele_data.cms_regardless_of_issues:
                authorised_by = cm.get('curation_verified_by', '')
                authorised_date = cm.get(SpecialEKeys.CURATION_VERIFIED_DATE, '')
                created_by = cm.get('curated_by', '')
                created_date = cm.get(SpecialEKeys.CURATION_DATE, '')
                patient_id = cm.get(SpecialEKeys.PATIENT_ID, '')
                interpretation_summary = cm.get(SpecialEKeys.INTERPRETATION_SUMMARY, '')
                cs = cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
                if patient_id != '':
                    patient_ids.add(patient_id)
                if cm.lab == lab_name1:
                    if cs and cs != 'None':
                        lab1.add(cs)
                    for key, value in cm.evidence.items():
                        lab1_comment_data[key].append(value)
                    authorised_date_lab1.add(authorised_by + ' ' + authorised_date)
                    created_date_lab1.add(created_by + ' ' + created_date)
                    interpretation_summary_lab1.add(interpretation_summary)
                if cm.lab == lab_name2:
                    if cs and cs != 'None':
                        lab2.add(cs)
                    for key, value in cm.evidence.items():
                        lab2_comment_data[key].append(value)
                    authorised_date_lab2.add(authorised_by + ' ' + authorised_date)
                    created_date_lab2.add(created_by + ' ' + created_date)
                    interpretation_summary_lab2.add(interpretation_summary)

            lab1_set = ','.join(sorted(lab1))
            lab2_set = ','.join(sorted(lab2))
            patient_id_set = ','.join(sorted(patient_ids))
            authorised_date_lab1_set = ','.join(sorted(authorised_date_lab1))
            authorised_date_lab2_set = ','.join(sorted(authorised_date_lab2))
            created_date_lab1_set = ','.join(sorted(created_date_lab1))
            created_date_lab2_set = ','.join(sorted(created_date_lab2))
            interpretation_summary_lab1_set = ','.join(sorted(interpretation_summary_lab1))
            interpretation_summary_lab2_set = ','.join(sorted(interpretation_summary_lab2))

            if lab1_comment_data and lab2_comment_data:
                for key in lab1_comment_data.keys():
                    lab1_values = lab1_comment_data[key]
                    lab2_values = lab2_comment_data[key]

                    comment_data_same = self.compare_comment_data(lab1_values, lab2_values)

                    if not comment_data_same:
                        message += f'{key}, '
            else:
                message = 'No Overlap'

            row = ClassificationlabCompareRow(allele_data=allele_data, patient_id=patient_id_set,
                                              lab_1=lab1_set,
                                              lab_2=lab2_set,
                                              lab_1_authorised_date=authorised_date_lab1_set,
                                              lab_2_authorised_date=authorised_date_lab2_set,
                                              lab_1_curated_date=created_date_lab1_set,
                                              lab_2_curated_date=created_date_lab2_set,
                                              lab_1_interpretation_summary=interpretation_summary_lab1_set,
                                              lab_2_interpretation_summary=interpretation_summary_lab2_set,
                                              comment=message
                                              )
            rows.append(delimited_row([row.allele_url(), row.patient_id(),
                                       row.lab_clinical_significance(),
                                       row.lab2_clinical_significance(),
                                       row.difference(lab_name1, lab_name2),
                                       row.lab_1_authorised_date,
                                       row.lab_2_authorised_date,
                                       row.authorised_date_diff(),
                                       row.lab_1_curated_date,
                                       row.lab_2_curated_date,
                                       row.created_date_diff(),
                                       row.comment,
                                       row.lab_1_interpretation_summary,
                                       row.lab_2_interpretation_summary,
                                       ], ','))
        else:
            raise ValueError("Error comparing labs")
        return rows

    def tokenize_text(self, text):
        tokens = re.findall(r'\b\w+\b', text.lower())
        return set(tokens)

    def compare_comment_data(self, lab1_values, lab2_values):
        lab1_tokens = set()
        lab2_tokens = set()

        for value in lab1_values:
            lab1_tokens.update(self.tokenize_text(str(value).replace('\\n', ' ')))

        for value in lab2_values:
            lab2_tokens.update(self.tokenize_text(str(value).replace('\\n', ' ')))

        return lab1_tokens == lab2_tokens

