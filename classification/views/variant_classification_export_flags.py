from collections import namedtuple
import csv
import io

from flags.models.models import Flag, FlagComment
from library.django_utils import get_url_from_view_path
from classification.enums.variant_classification_enums import SpecialEKeys
from classification.models.evidence_key import EvidenceKeyMap
from classification.models.variant_classification import VariantClassificationModification
from classification.views.variant_classification_export_utils import ExportFormatter


Problem = namedtuple('Problem', 'code message')


class ExportFormatterFlags(ExportFormatter):

    CODE_SUBSTITUTE = {'mandatory': 'Missing value'}

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def header(self) -> str:
        # hard-code columns that aren't in evidence

        header = ['url', 'lab', 'lab_record_id', 'imported_chgvs', 'last updated', 'problem_code', 'problem']
        return ExportFormatter.write_single_row(header, delimiter=',')

    def iterate_problems(self, vcm: VariantClassificationModification):
        evidence = vcm.evidence
        # validation messages
        for key, blob in evidence.items():
            validations = blob.get('validation')
            if validations:
                for validation in validations:
                    key_label = self.ekeys.get(key).pretty_label
                    code = validation.get('code', '')
                    code = self.CODE_SUBSTITUTE.get(code, code)
                    code = (code[0].upper() + code[1:]).replace('_', ' ')
                    message = validation.get('message')
                    yield Problem(code=code, message=f'{key_label}: {message}')
        # flags
        open_flags = vcm.variant_classification.flag_collection_safe.flags(only_open=True)\
            .select_related('flag_type', 'resolution')
        flag: Flag
        for flag in open_flags:
            flag_type = flag.flag_type
            parts = [flag_type.label]
            if flag.resolution_id != 'open':
                parts.append(flag.resolution.label)
            last_comment = FlagComment.last(flag)
            if last_comment and last_comment.text:
                parts.append(last_comment.text)

            yield Problem(code=flag_type.label, message='\n'.join(parts))

    def row(self, vcm: VariantClassificationModification):
        vc = vcm.variant_classification
        url = get_url_from_view_path(vc.get_absolute_url())
        lab = vc.lab.name
        lab_record_id = vc.lab_record_id
        imported_chgvs = vc.get(SpecialEKeys.C_HGVS, None)
        version = vcm.created.strftime("%Y-%m-%d %H:%M")

        out = io.StringIO()
        writer = csv.writer(out, delimiter=',')

        for problem in self.iterate_problems(vcm):
            row = [url, lab, lab_record_id, imported_chgvs, version, problem.code, problem.message]
            writer.writerow(row)

        return out.getvalue()

    def prepare_groups(self):
        pass

    def row_iterator(self):
        return self.qs

    def filename(self) -> str:
        return self.generate_filename(include_genome_build=False, suffix='issues')
