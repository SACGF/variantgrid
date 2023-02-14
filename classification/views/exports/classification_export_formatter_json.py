import json
from functools import cached_property
from typing import List

from django.http import HttpRequest

from classification.models import ClassificationJsonParams, ClassificationModification
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter


@register_classification_exporter("json")
class ClassificationExportFormatterJSON(ClassificationExportFormatter):

    def __init__(self, classification_filter: ClassificationFilter):
        self.first_row = True
        super().__init__(classification_filter=classification_filter)

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterJSON':
        return ClassificationExportFormatterJSON(
            classification_filter=ClassificationFilter.from_request(request)
        )

    @cached_property
    def json_params(self) -> ClassificationJsonParams:
        return ClassificationJsonParams(current_user=self.classification_filter.user,
                                        include_data=True,
                                        api_version=2,
                                        strip_complicated=True,
                                        include_messages=False)

    def header(self) -> List[str]:
        return ['{"records":[']

    def row(self, allele_data: AlleleData) -> List[str]:
        rows = []
        for ci in allele_data.all_cms:
            if row := self.to_row(ci.classification, withdrawn=ci.withdrawn):
                if not self.first_row:
                    row = f",{row}"
                else:
                    self.first_row = False
                rows.append(row)
        return rows

    def footer(self) -> List[str]:
        return ["]}"]

    def to_row(self, vcm: ClassificationModification, withdrawn: bool) -> str:
        json_values = vcm.as_json(self.json_params)
        if 'fatal_error' in json_values:
            return None

        if self.classification_filter.is_discordant(vcm):
            json_values['discordant'] = True
        if withdrawn:
            json_values = {
                "id": json_values["id"],
                "delete": True
            }

        return json.dumps(json_values)

    def extension(self) -> str:
        return "json"

    def content_type(self) -> str:
        return "application/json"
