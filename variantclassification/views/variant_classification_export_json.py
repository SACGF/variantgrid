import json
from typing import Optional

from variantclassification.models.variant_classification import VariantClassificationModification
from variantclassification.models import VariantClassificationJsonParams
from variantclassification.views.variant_classification_export_utils import ExportFormatter,\
    AlleleGroup


class ExportFormatterJSON(ExportFormatter):
    """
    Formats as JSON
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.first_row = True

    def header(self) -> str:
        return '{"records":[\n'

    def footer(self) -> str:
        return '\n]}'

    def to_row(self, vcm: VariantClassificationModification) -> Optional[str]:
        params = VariantClassificationJsonParams(current_user=self.user,
                                                 include_data=True,
                                                 api_version=2,
                                                 strip_complicated=True,
                                                 include_messages=False)
        json_values = vcm.as_json(params)
        if 'fatal_error' in json_values:
            return None

        if self.is_discordant(vcm.variant_classification):
            json_values['discordant'] = True
        json_str = json.dumps(json_values)

        if self.first_row:
            self.first_row = False
            return json_str
        return ',\n' + json_str

    def row(self, group: AlleleGroup) -> str:
        row_str = ''
        for vcm in group.data:
            row_str += self.to_row(vcm)
        return row_str

    def content_type(self) -> str:
        return 'application/json'

    def filename(self) -> str:
        return self.generate_filename(extension='json')
