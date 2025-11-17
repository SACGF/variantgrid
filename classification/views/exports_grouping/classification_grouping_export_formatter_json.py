from dataclasses import dataclass
from functools import cached_property
from typing import Iterator

from django.conf import settings

from classification.models import ClassificationJsonParams, EvidenceKeyMap
from classification.views.exports_grouping.classification_grouping_export_filter import \
    ClassificationGroupingExportFormat, ClassificationGroupingExportFormatProperties, ClassificationGroupingExportFilter


@dataclass(frozen=True)
class JSONFormatDetails:
    full_detail = False


class ClassificationGroupingExportFormatterJSON(ClassificationGroupingExportFormat):

    @classmethod
    def format_properties(cls) -> ClassificationGroupingExportFormatProperties:
        return ClassificationGroupingExportFormatProperties(
            http_content_type="text/json",
            extension="json",
            delimiter_for_row=",\n"
        )

    def __init__(self,
                 classification_grouping_filter: ClassificationGroupingExportFilter,
                 json_format_details: JSONFormatDetails = JSONFormatDetails()):
        self.json_format_details = json_format_details
        super().__init__(classification_grouping_filter)

    def header(self) -> list[str]:
        return ['{"records:[']

    @cached_property
    def json_params(self):
        include_data: tuple[bool, set[str]]
        populate_literature_with_citations = False
        if self.json_format_details.full_detail:
            include_data = True
        else:
            e_keys = EvidenceKeyMap.cached()
            include_data = [e_key.key for e_key in self.e_keys.all_keys if e_key.is_vital_key]
            populate_literature_with_citations = settings.CLASSIFICATION_DOWNLOADABLE_JSON_LITERATURE_CITATIONS

        return ClassificationJsonParams(
            current_user=self.classification_grouping_filter.user,
            api_version=2,
            include_data=include_data,
            populate_literature_with_citations=populate_literature_with_citations
        )

    def single_row_generator(self) -> Iterator[str]:
        for cg in self.queryset().iterator():
            yield cg.latest_classification_modification.as_json(
                self.json_params
            )

    def footer(self) -> list[str]:
        return [']}']