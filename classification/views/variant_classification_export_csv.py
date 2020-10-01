import csv
import io
from typing import Optional

from genes.hgvs import CHGVS
from classification.models.evidence_key import EvidenceKeyMap
from classification.models.variant_classification import VariantClassification, \
    VariantClassificationModification
from classification.views.variant_classification_export_utils import ExportFormatter, \
    AlleleGroup, UsedKeyTracker, \
    KeyValueFormatter


class ExportFormatterCSV(ExportFormatter):
    """
    Formats as a general CSV
    """

    def __init__(self, filename_override:str = None, pretty=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.used_keys = UsedKeyTracker(self.user, self.ekeys, KeyValueFormatter(), pretty=pretty)
        self.filename_override = filename_override

    def header(self) -> str:
        qs = self.qs
        for vcm in qs.all():
            self.used_keys.check_record(vcm)

        # hard-code columns that aren't in evidence

        header = ['id', 'lab', 'lab_record_id', 'share_level', 'version', 'liftover_error', 'target_genome_build', 'target_c_hgvs', 'acmg_criteria', 'evidence_weights', 'citations', 'in_discordance'] + self.used_keys.header()
        return ExportFormatter.write_single_row(header, delimiter=',')

    def to_row(self, c_parts: Optional[CHGVS], vcm: VariantClassificationModification, message=None) -> list:
        vc = vcm.variant_classification

        acmg_criteria = vcm.criteria_strength_summary(self.ekeys)
        evidence_weights = VariantClassification.summarize_evidence_weights(vcm.evidence, self.ekeys)
        citations = ', '.join([c.ref_id() for c in vcm.citations])

        full_chgvs = None
        if c_parts:
            full_chgvs = c_parts.full_c_hgvs

        row = [
            vc.id,
            vc.lab.name,  # row
            vc.lab_record_id,
            vcm.share_level_enum.label,
            vcm.created.timestamp(),
            message,
            self.genome_build.name,
            full_chgvs,
            acmg_criteria,
            evidence_weights,
            citations,
            'TRUE' if self.is_discordant(vc) else 'FALSE',
        ] + self.used_keys.row(variant_classification_modification=vcm)
        return row

    def row(self, group: AlleleGroup) -> str:
        out = io.StringIO()
        writer = csv.writer(out, delimiter=',')

        #url = get_url_from_view_path( group.target_variant.get_absolute_url() )
        for c_parts, vcms in group.iter_c_hgvs():

            for vcm in vcms:
                writer.writerow(self.to_row(c_parts, vcm))

        return out.getvalue()

    def footer(self):
        if not self.error_message_ids:
            return None

        out = io.StringIO()
        writer = csv.writer(out, delimiter=',')
        for vcm, message in self.record_and_error():
            writer.writerow(self.to_row(None, vcm, message))

        return out.getvalue()

    def filename(self) -> str:
        if self.filename_override:
            return self.filename_override
        return self.generate_filename()
