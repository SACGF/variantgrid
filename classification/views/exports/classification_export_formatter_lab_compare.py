import re
from collections import defaultdict
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

    def __init__(self,  AlleleData, lab1, lab2, comment=None, authorised_date=None, created_date=None, Patient_id=None):
        self.AlleleData = AlleleData
        self.lab1 = lab1
        self.lab2 = lab2
        self.comment = comment
        self.authorised_date = authorised_date
        self.created_date = created_date
        self.Patient_id = Patient_id

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

    @export_column("Clinical Significance")
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

    @export_column("Authorised Date")
    def authorised_date(self):
        return self.authorised_date

    @export_column("Created Date")
    def created_date(self):
        return self.created_date

    @export_column("Other Variables")
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
            return [delimited_row(['Allele URL', 'Patient_Id', lab_names[0], lab_names[1], 'Clinical Significance',
                                   'Authorised Date', 'Created Date', 'Other Variables'], ',')]
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

        if not allele_data.allele_id:
            return []
        if self.classification_filter.include_sources:
            lab_name1, lab_name2 = sorted(self.classification_filter.include_sources)
            for cm in allele_data.cms_regardless_of_issues:
                authorised_by = cm.get('curation_verified_by', '')
                authorised_date = cm.get('curation_verified_date', '')
                created_by = cm.get('curated_by', '')
                created_date = cm.get('curation_date', '')
                patient_id = cm.get('patient_id', '')
                if patient_id != '':
                    patient_ids.add(patient_id)
                if cm.lab == lab_name1:
                    if cm.clinical_significance is not None:
                        lab1.add(cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
                    for key, value in cm.evidence.items():
                        lab1_comment_data[key].append(value)
                    authorised_date_lab1.add(authorised_by + ' ' + authorised_date)
                    created_date_lab1.add(created_by + ' ' + created_date)
                if cm.lab == lab_name2:
                    if cm.clinical_significance is not None:
                        lab2.add(cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
                    for key, value in cm.evidence.items():
                        lab2_comment_data[key].append(value)
                    authorised_date_lab2.add(authorised_by + ' ' + authorised_date)
                    created_date_lab2.add(created_by + ' ' + created_date)

            lab1_set = ','.join(lab1)
            lab2_set = ','.join(lab2)
            patient_id_set = ','.join(patient_ids)
            authorised_date_lab1_set = ','.join(authorised_date_lab1)
            authorised_date_lab2_set = ','.join(authorised_date_lab2)
            created_date_lab1_set = ','.join(created_date_lab1)
            created_date_lab2_set = ','.join(created_date_lab2)

            if authorised_date_lab1_set != authorised_date_lab2_set:
                difference_authorised_by = 'Different'
            else:
                difference_authorised_by = 'Same'
            if created_date_lab1_set != created_date_lab2_set:
                difference_created_by = 'Different'
            else:
                difference_created_by = 'Same'

            for key in lab1_comment_data.keys():
                lab1_values = lab1_comment_data[key]
                lab2_values = lab2_comment_data[key]

                comment_data_same = self.compare_comment_data(lab1_values, lab2_values)

                if not comment_data_same:
                    message += f'{key}, '

            row = ClassificationlabCompareRow(AlleleData=allele_data.allele_id, Patient_id=patient_id_set, lab1=lab1_set, lab2=lab2_set,
                                              authorised_date=difference_authorised_by,
                                              created_date=difference_created_by, comment=message

                                              )
            rows.append(delimited_row([row.allele_url(), row.patient_id(), row.lab_clinical_significance(),
                                       row.lab2_clinical_significance(), row.difference(lab_name1, lab_name2),
                                       row.authorised_date, row.created_date,
                                       row.comment], ','))
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

