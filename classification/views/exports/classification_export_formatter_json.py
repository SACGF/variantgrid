import json
from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from django.conf import settings
from django.http import HttpRequest

from classification.models import ClassificationJsonParams, ClassificationModification, EvidenceKeyMap
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData, \
    DiscordanceReportStatus
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter


@dataclass(frozen=True)
class FormatDetailsJSON:
    # format evidence keys to nice human labels or leave as raw codes easier handled by code
    full_detail: bool = False

    @staticmethod
    def from_request(request: HttpRequest) -> 'FormatDetailsCSV':
        full_detail = settings.CLASSIFICATION_DOWNLOADABLE_NOTES_AND_EXPLAINS or (request.query_params.get('full_detail') == 'true' and request.user.is_superuser)

        return FormatDetailsJSON(
            full_detail=full_detail
        )


@register_classification_exporter("json")
class ClassificationExportFormatterJSON(ClassificationExportFormatter):

    def __init__(self, classification_filter: ClassificationFilter, format_details: FormatDetailsJSON):
        self.e_keys = EvidenceKeyMap.cached()
        self.format_details = format_details
        super().__init__(classification_filter=classification_filter)



    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterJSON':

        return ClassificationExportFormatterJSON(
            classification_filter=ClassificationFilter.from_request(request),
            format_details=FormatDetailsJSON.from_request(request)
        )

    @cached_property
    def json_params(self) -> ClassificationJsonParams:
        include_data: tuple[bool, set[str]]
        if self.format_details.full_detail:
            include_data = True
        else:
            include_data = [e_key.key for e_key in self.e_keys.all_keys if e_key.is_vital_key]

        return ClassificationJsonParams(current_user=self.classification_filter.user,
                                        api_version=2,
                                        strip_complicated=True,
                                        strip_notes_and_explains=not self.format_details.full_detail,
                                        include_messages=False,
                                        include_data=include_data,
                                        #include_data = True,
                                        # FIXME don't want this to be the default
                                        remove_acmg_namespace=True)

    @property
    def delimiter_for_row(self):
        return ","

    def header(self) -> list[str]:
        return ['{"records":[']

    def row(self, allele_data: AlleleData) -> list[str]:
        rows = []
        for ci in allele_data.all_cms:
            if row := self.to_row(ci.classification, withdrawn=ci.withdrawn):
                rows.append(row)
        return rows

    def footer(self) -> list[str]:
        return ["]}"]

    def get_discordant_status(self, discordant_status):
        if not discordant_status:
            return ''
        elif discordant_status == DiscordanceReportStatus.ON_GOING:
            return "active discordance"
        elif discordant_status == DiscordanceReportStatus.CONTINUED:
            return "continued discordance"
        elif discordant_status == DiscordanceReportStatus.PENDING_CONCORDANCE:
            return "pending concordance"
        else:
            return discordant_status

    def to_row(self, vcm: ClassificationModification, withdrawn: bool) -> Optional[str]:
        json_values = vcm.as_json(self.json_params)
        if 'fatal_error' in json_values:
            return None

        if discordant_status := self.classification_filter.is_discordant(vcm):
            json_values['discordant'] = True
            json_values['discordant_status'] = self.get_discordant_status(discordant_status)
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
