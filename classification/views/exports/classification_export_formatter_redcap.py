import csv
import io
import re
from typing import Iterable

from django.db.models import OrderBy
from django.db.models.expressions import RawSQL
from django.http import HttpRequest, HttpResponse
from django.http.response import HttpResponseBase

from classification.enums import EvidenceKeyValueType
from classification.models import EvidenceKey, Classification, EvidenceKeyMap, ClassificationModification
from classification.views.classification_export_utils import KeyValueFormatter, UsedKeyTracker
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from classification.views.exports.classification_export_utils import CitationCounter
from library.utils import delimited_row


# WARNING: REDCap export ignores the since parameter
# WARNING: REDCap export ignores variant warnings

class RedcapGroup:
    """
    A bunch of records linked to variants with the same allele
    Also contains a "target_variant" from the desired genome_build
    """

    def __init__(self, group_on_value):
        self.group_on_value = group_on_value
        self.data: list[ClassificationModification] = []

    def __bool__(self):
        return bool(self.data)


@register_classification_exporter("redcap")
class ClassificationExportFormatterRedCap(ClassificationExportFormatter):
    """
    Exports data in the format that RedCap can import it.
    @see https://www.project-redcap.org/
    """

    def __init__(self, classification_filter: ClassificationFilter):
        self.used_key_arrays = []
        self.group_on = 'redcap_record_id'
        self.e_keys = EvidenceKeyMap.cached()
        super().__init__(classification_filter=classification_filter)
        for group in self.redcap_rows():
            for idx, vcm in enumerate(group.data):
                used_keys: UsedKeyTracker
                if len(self.used_key_arrays) < (idx + 1):
                    used_keys = UsedKeyTracker(self.classification_filter.user, self.e_keys, RedcapKeyValueFormatter(index=idx + 1))
                    self.used_key_arrays.append(used_keys)
                else:
                    used_keys = self.used_key_arrays[idx]
                used_keys.check_record(vcm)

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterRedCap':
        classification_filter = ClassificationFilter.from_request(request)
        return ClassificationExportFormatterRedCap(
            classification_filter=classification_filter
        )

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> list[str]:
        out = io.StringIO()
        writer = csv.writer(out, delimiter=',')

        header_array = [REDCAP_PREFIX + 'record_id', REDCAP_PREFIX + 'count']
        standard_keys = ['id', 'lab', 'lab_record_id', 'version', 'evidence_weights', 'citations']

        for used_keys in self.used_key_arrays:
            header_array = header_array + [REDCAP_PREFIX + key + '_' + str(used_keys.key_value_formatter.index) for key in standard_keys]
            header_array = header_array + used_keys.header()

        writer.writerow(header_array)
        return [out.getvalue()]

    def row(self, allele_data: AlleleData) -> list[str]:
        raise ValueError("row() is not supported in REDCap Export")

    def redcap_rows(self) -> Iterable[RedcapGroup]:
        qs = self.classification_filter.cms_qs
        qs = qs.exclude(**{f'published_evidence__{self.group_on}__value__isnull': True})
        # Have to use RawSQL to get to order by jsonb values, you can include them in a filter/exclude
        # but they just don't work in order_by normally
        qs = qs.order_by(OrderBy(RawSQL('cast(published_evidence->>%s as jsonb)->>%s', (self.group_on, 'value'))))

        group = None
        for vcm in qs.iterator():
            value = vcm.get(self.group_on, None)
            if group is None or group.group_on_value != value:
                if group:
                    yield group
                group = RedcapGroup(value)
            group.data.append(vcm)
        if group:
            yield group

    def row_generator(self) -> Iterable[list[str]]:
        for group in self.redcap_rows():
            row_values = [group.group_on_value, len(group.data)]
            for idx, vcm in enumerate(group.data):
                vc = vcm.classification

                evidence_weights = Classification.summarize_evidence_weights(vcm.evidence, self.e_keys)
                citation_count = CitationCounter()
                citation_count.reference_citations(vcm)
                citations = ', '.join(citation_count.citation_ids())

                row_values = row_values + [vc.id, vc.lab.name, vc.lab_record_id, vcm.created.timestamp(),
                                           evidence_weights, citations]
                row_values = row_values + self.used_key_arrays[idx].row(vcm)

            yield [delimited_row(row_values, delimiter=',')]


REDCAP_PREFIX = 'vc_'


class RedcapDefinition:

    @staticmethod
    def header_row() -> list[str]:
        return [
            'Variable / Field Name',
            'Form Name',
            'Section Header',
            'Field Type',
            'Field Label',
            'Choices, Calculations, OR Slider Labels',

            'Field Note',

            'Text Validation Type OR Show Slider Number',

            'Text Validation Min',
            'Text Validation Max',
            'Identifier',
            'Branching Logic (Show field only if...)',

            'Required Field?',

            'Custom Alignment',
            'Question Number (surveys only)',
            'Matrix Group Name',
            'Matrix Ranking?',
            'Field Annotation'
        ]

    def __init__(self, name=None, form_name=None, section_header=None, field_type=None, label=None, logic=None):
        self.name = name
        self.form_name = form_name
        self.section_header = section_header
        self.field_type = field_type
        self.label = label
        self.choices = None
        self.note = None
        self.text_validation = None
        self.text_validation_min = None
        self.text_validation_max = None
        self.identifier = None
        self.branching_logic = logic
        self.required_field = None
        self.custom_alignment = None
        self.question_number = None
        self.matrix_group_name = None
        self.matrix_ranking = None
        self.annotation = None

    @staticmethod
    def safe_key(key: str) -> str:
        # redcap has a limit of 26 characters for a key
        # we prefix vc_ bringing that down to 23
        # we can suffix with _1, _2, ... _10 and notes _n1 to _n10 bringing that down to 19
        return key[:19]

    def row(self):
        return [
            self.name,
            self.form_name,
            self.section_header,
            self.field_type,
            self.label,
            self.choices,
            self.note,
            self.text_validation,
            self.text_validation_min,
            self.text_validation_max,
            self.identifier,
            self.branching_logic,
            self.required_field,
            self.custom_alignment,
            self.question_number,
            self.matrix_group_name,
            self.matrix_ranking,
            self.annotation
        ]

def redcap_data_dictionary(request: HttpRequest) -> HttpResponseBase:
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="redcap_data_definition.csv"'
    export_redcap_definition(response)
    return response

def export_redcap_definition(output):

    repeat_count = 10
    form_name = 'shariant'

    csv_writer = csv.writer(output)

    csv_writer.writerow(RedcapDefinition.header_row())

    dont_group = {'redcap_record_id'}

    csv_writer.writerow(RedcapDefinition(name=REDCAP_PREFIX + 'record_id', form_name=form_name, section_header='Record', field_type='text', label='Redcap Record ID').row())
    csv_writer.writerow(RedcapDefinition(name=REDCAP_PREFIX + 'count', form_name=form_name, section_header='', field_type='text', label='Count').row())
    for idx in range(repeat_count):
        suffix = '_' + str(idx + 1)
        show_if = '[%scount] >= %s' % (REDCAP_PREFIX, str(idx + 1))

        label = '[' + REDCAP_PREFIX + 'c_hgvs' + suffix + '] ' + \
                '([' + REDCAP_PREFIX + 'clinical_significance' + suffix + '])'
        csv_writer.writerow(RedcapDefinition(name=REDCAP_PREFIX + 'var_des' + suffix, form_name=form_name, section_header='', field_type='descriptive', label='Variant' + str(idx + 1) + ': ' + label, logic=show_if).row())

    key_map = EvidenceKeyMap.instance()
    for idx in range(repeat_count):
        suffix = '_' + str(idx + 1)
        show_if = '[%scount] >= %s' % (REDCAP_PREFIX, str(idx + 1))

        csv_writer.writerow(RedcapDefinition(name=REDCAP_PREFIX + 'lab' + suffix, form_name=form_name, section_header='Ids' + suffix, field_type='text', label='Lab', logic=show_if).row())
        csv_writer.writerow(RedcapDefinition(name=REDCAP_PREFIX + 'lab_record_id' + suffix, form_name=form_name, field_type='text', label='Lab Record ID', logic=show_if).row())
        csv_writer.writerow(RedcapDefinition(name=REDCAP_PREFIX + 'id' + suffix, form_name=form_name, field_type='text', label='Shariant ID', logic=show_if).row())
        csv_writer.writerow(RedcapDefinition(name=REDCAP_PREFIX + 'version' + suffix, form_name=form_name, field_type='text', label='Record Version', logic=show_if).row())
        csv_writer.writerow(RedcapDefinition(name=REDCAP_PREFIX + 'evidence_weights' + suffix, form_name=form_name, field_type='text', label='Evidence Weights', logic=show_if).row())
        csv_writer.writerow(RedcapDefinition(name=REDCAP_PREFIX + 'citations' + suffix, form_name=form_name, field_type='text', label='Citations', logic=show_if).row())

        last_evidence_category = None
        for key_entry in key_map.all_keys:
            if key_entry.key in dont_group:
                continue

            row = RedcapDefinition(
                name=REDCAP_PREFIX + key_entry.redcap_key(suffix=idx + 1),
                form_name=form_name,
                label=key_entry.pretty_label,
                logic=show_if)

            row_note = RedcapDefinition(
                name=REDCAP_PREFIX + key_entry.redcap_key(suffix=idx + 1, is_note=True),
                form_name=form_name,
                label=key_entry.pretty_label + ' Note',
                field_type='notes',
                logic=show_if)

            if key_entry.evidence_category != last_evidence_category:
                last_evidence_category = key_entry.evidence_category
                row.section_header = key_entry.get_evidence_category_display() + suffix

            validation = ''

            value_type = key_entry.value_type
            redcap_type = 'text'
            if value_type == EvidenceKeyValueType.BOOLEAN:
                redcap_type = 'yesno'
            elif value_type == EvidenceKeyValueType.SELECT or value_type == EvidenceKeyValueType.CRITERIA:
                redcap_type = 'dropdown'
            elif value_type == EvidenceKeyValueType.TEXT_AREA:
                redcap_type = 'notes'
            elif value_type == EvidenceKeyValueType.DATE:
                redcap_type = 'text'
                validation = 'date_dmy'

            row.field_type = redcap_type

            option_text = ''
            if value_type == EvidenceKeyValueType.SELECT or value_type == EvidenceKeyValueType.CRITERIA:
                options = key_entry.virtual_options
                if options:
                    for option in options:
                        if option.get('index') is not None:
                            if option_text:
                                option_text += ' | '
                            option_text += str(option.get('index')) + ', '
                            if option.get('label'):
                                # FIXME escape any | in the label
                                option_text += option.get('label')
                            else:
                                option_text += EvidenceKey.pretty_label_from_string(option.get('key'))

            row.choices = option_text
            row.text_validation = validation

            if key_entry.mandatory:
                row.required_field = 'y'

            csv_writer.writerow(row.row())
            csv_writer.writerow(row_note.row())


class RedcapKeyValueFormatter(KeyValueFormatter):

    DATE_FORMAT_RE = re.compile('([0-9]+)-([0-9]+)-([0-9]+)')

    def __init__(self, index: int):
        self.index = index

    def header_for(self, ekey: EvidenceKey, is_note: bool = False, pretty=False) -> str:
        return REDCAP_PREFIX + ekey.redcap_key(self.index, is_note=is_note)

    def value_for(self, ekey: EvidenceKey, value, pretty=False):
        if ekey.value_type in {EvidenceKeyValueType.SELECT, EvidenceKeyValueType.CRITERIA}:
            if value:
                return Classification._export_value(value, key_data=ekey, export_key='index')

        elif ekey.value_type == EvidenceKeyValueType.MULTISELECT:
            return ekey.pretty_value(value)

        elif ekey.value_type == EvidenceKeyValueType.DATE:
            if value:
                match = RedcapKeyValueFormatter.DATE_FORMAT_RE.match(value)
                if match:
                    return match[3] + '/' + match[2] + '/' + match[1]
            return None

        elif value is True:
            return 'TRUE'

        elif value is False:
            return 'FALSE'

        elif value is not None:
            return str(value)

        return None
