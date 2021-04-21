import json
from datetime import datetime
from typing import Optional
from lazy import lazy
from classification.models.classification import ClassificationModification
from classification.models import ClassificationJsonParams
from classification.views.classification_export_utils import ExportFormatter,\
    AlleleGroup


class ExportFormatterJSON(ExportFormatter):
    """
    Formats as JSON
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.first_row = True

    @property
    def supports_fully_withdrawn(self) -> bool:
        return True

    def header(self) -> str:
        return '{"records":[\n'

    def footer(self) -> str:
        return '\n]}'

    @lazy
    def json_params(self) -> ClassificationJsonParams:
        return ClassificationJsonParams(current_user=self.user,
                                                 include_data=True,
                                                 api_version=2,
                                                 strip_complicated=True,
                                                 include_messages=False)

    def to_row(self, vcm: ClassificationModification, withdrawn=False) -> Optional[str]:
        json_values = vcm.as_json(self.json_params)
        if 'fatal_error' in json_values:
            return None

        if self.is_discordant(vcm.classification):
            json_values['discordant'] = True
        if withdrawn:
            json_values['delete'] = True

        json_str = json.dumps(json_values)

        if self.first_row:
            self.first_row = False
            return json_str

        self.row_count += 1
        return ',\n' + json_str

    def row(self, group: AlleleGroup) -> str:
        row_str = ''
        for vcm in group.data:
            row_str += self.to_row(vcm)
        for vcm in group.withdrawn:
            row_str += self.to_row(vcm, withdrawn=True)
        return row_str

    def content_type(self) -> str:
        return 'application/json'

    def filename(self) -> str:
        return self.generate_filename(extension='json')
