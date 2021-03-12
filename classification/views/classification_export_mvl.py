import csv
import io
from typing import Dict, List, Tuple, Union

from django.conf import settings
from django.utils.timezone import now
from pytz import timezone

from annotation.citations import get_citations, CITATION_COULD_NOT_LOAD_TEXT, CitationDetails
from annotation.views import simple_citation_html
from genes.hgvs import CHGVS
from library.django_utils import get_url_from_view_path
from classification.enums.classification_enums import SpecialEKeys
from classification.models import EvidenceKey
from classification.models.evidence_key import EvidenceKeyMap
from classification.regexes import db_ref_regexes
from classification.views.classification_export_utils import ExportFormatter, \
    AlleleGroup, ConflictStrategy, VariantWithChgvs


class ExportFormatterMVL(ExportFormatter):
    """
    Formats classifications for Agilent 5.2 MVL usage
    """

    @property
    def version(self):
        return '2.1'

    @property
    def use_full_chgvs(self):
        return True

    def __init__(self, conflict_strategy: str, cs_override_labels: Dict[str, str], *args, **kwargs):
        super().__init__( *args, **kwargs)
        self.conflict_strategy = conflict_strategy
        self.cs_translator = {
            'B': 'BENIGN',
            'LB': 'LIKELY_BENIGN',
            'VUS': 'VOUS',
            'LP': 'LIKELY_PATHOGENIC',
            'P': 'PATHOGENIC'
        }
        self.cs_translator.update(cs_override_labels)
        # VUS label will be used as the default for any unknown values
        # useful for VUS_A, VUS_B, VUS_C
        # slightly dodgy (but only option) for Risk Factor, Drug Response etc
        self.vous_label = self.cs_translator.get('VUS')

    RAW_SCORE = {
        'B': 1,
        'LB': 2,
        'VUS': 3,
        'VUS_A': 3,
        'VUS_B': 3,
        'VUS_C': 3,
        'LP': 4,
        'P': 5
    }
    DEFAULT_SCORE = 3

    def header(self) -> str:
        return '\t'.join(['transcript', 'c_nomen', 'classification', 'variant information', 'report abstract']) + '\n'

    @staticmethod
    def mvl_safe(value: str) -> str:
        if value is None:
            return ''

        return str(value).replace('<', '&lt;').replace('\n', '<br>').replace('\t', '&emsp')

    def summarise_vcm(self, vcmc: VariantWithChgvs, row_chgvs: CHGVS) -> Tuple[str, str]:
        vcm = vcmc.vcm
        vc = vcm.classification

        def format_label(label: str) -> str:
            return f'<span style="color:#888">{label} </span> '

        def get_value(ekey: EvidenceKey) -> str:
            return ExportFormatterMVL.mvl_safe(ekey.pretty_value(vcm.get(ekey.key))) or 'None'

        def format_key(ekey: EvidenceKey) -> str:
            nonlocal vcm

            label = ekey.pretty_label
            value = get_value(ekey)
            value = db_ref_regexes.link_html(value)
            return f'{format_label(label)} {value}'

        def citation_title(cd: CitationDetails) -> str:
            if cd.title == CITATION_COULD_NOT_LOAD_TEXT:
                return ''
            else:
                return cd.title

        # no longer provide link to individual classifications (as you might then miss the overall changes)
        record_label = vc.friendly_label

        notes: List[str] = list()

        if self.is_discordant(vc):
            notes.append('<span style="color:#c44">This record is in discordance</span>')

        my_raw_c = vcmc.chgvs.raw_c
        row_raw_c = row_chgvs.raw_c
        if my_raw_c != row_raw_c:
            notes.append('<span style="color:#c44">Warning c.hgvs representation is different across transcript versions</span>')

        # we don't want to display 1,000,000 etc nucleotides to the user, show the user friendly short version
        short_hgvs = format_label('c.hgvs') + vc.get_c_hgvs(self.genome_build, use_full=False)

        cs = format_key(self.ekeys.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
        curated = format_key(self.ekeys.get(SpecialEKeys.CURATION_DATE))
        condition = vcmc.vcm.condition_text

        criteria_summary = format_label('Criteria met') + (vcm.criteria_strength_summary(self.ekeys) or 'None')
        #citation_anchors = [html_link(cd.citation_link, f'{cd.source}:{cd.citation_id}') + ' ' + citation_title(cd) for cd in get_citations(vcm.citations)]
        citation_anchors = [simple_citation_html(cd) for cd in get_citations(vcm.citations)]
        # interpretation summary is too large to include for every record
        #summary = format_key(self.ekeys.get(SpecialEKeys.INTERPRETATION_SUMMARY))
        citations_string = format_label('Citations') + (('<br>' + '<br>'.join(citation_anchors)) if citation_anchors else 'None')

        short_summary = f"{format_label(record_label)} {get_value(self.ekeys.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))}"
        record_label_bold = f'<b>{record_label}</b>'
        note = '<br>'.join(notes)
        long_summary = '<br>'.join(x for x in [record_label_bold, short_hgvs, note, cs, curated, condition, criteria_summary, citations_string] if x is not None)
        return short_summary, long_summary

    def row(self, group: AlleleGroup) -> str:
        out = io.StringIO()
        writer = csv.writer(out, delimiter='\t')

        date_str = now().astimezone(tz=timezone(settings.TIME_ZONE)).strftime("%Y-%m-%d")
        url = get_url_from_view_path(group.target_variant.get_absolute_url()) + f'?refer=mvl&seen={date_str}'
        variant_details = f'<a href="{url}" target="_blank">Click here for up-to-date classifications on this variant.</a>'

        for c_parts, vcms_w_chgvs in group.iter_c_hgvs_versionless_transcripts():

            transcript = c_parts.transcript
            c_hgvs = c_parts.raw_c
            # could probably do <th style="text-align:left">
            # but
            summaries = []
            abstracts = []
            warnings: List[str] = []
            using_classification_score = None
            classification = ''
            different_strengths = set()
            has_special_cs = False
            discordant_count = 0
            has_diff_chgvs = False

            for vcm_w_chgvs in vcms_w_chgvs:

                vcm = vcm_w_chgvs.vcm
                vc = vcm.classification
                if self.is_discordant(vc):
                    discordant_count += 1

                raw_classification = vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
                label = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(raw_classification) or 'Unclassified'
                different_strengths.add(label)
                this_classification = self.cs_translator.get(raw_classification, self.vous_label)
                this_classification_score = ExportFormatterMVL.RAW_SCORE.get(raw_classification, ExportFormatterMVL.DEFAULT_SCORE)
                has_special_cs = has_special_cs or raw_classification not in ExportFormatterMVL.RAW_SCORE

                if using_classification_score is None or \
                    (self.conflict_strategy == ConflictStrategy.MOST_BENIGN and this_classification_score < using_classification_score) or \
                    (self.conflict_strategy == ConflictStrategy.MOST_PATHOGENIC and this_classification_score > using_classification_score):
                    using_classification_score = this_classification_score
                    classification = this_classification

                summary, abstract = self.summarise_vcm(vcm_w_chgvs, c_parts)
                summaries.append(summary)
                abstracts.append(abstract)

                if vcm_w_chgvs.chgvs.raw_c != c_hgvs:
                    has_diff_chgvs = True

            if has_diff_chgvs:
                warnings.append('Warning <b>c.hgvs representations are different across transcript versions</b>')
            if has_special_cs:
                warnings.append('Warning <b>Contains non-standard clinical significance</b>')
            if discordant_count:
                if discordant_count == 1:
                    warnings.append(f'Warning <b>1 record is in discordance</b>')
                else:
                    warnings.append(f'Warning <b>{discordant_count} records are in discordance</b>')
            if len(different_strengths) > 1:
                strength_list = list(different_strengths)
                strength_list.sort()
                strength_list = ', '.join(strength_list)
                warnings.append(f'Warning <b>Multiple clinical significances recorded for this transcript : {strength_list}</b>')

            warning_text = '<br>'.join(warnings)
            if warning_text:
                warning_text = warning_text + '<br/>'

            complete_summary = '<ul>' + ''.join(f'<li>{summary}</li>' for summary in summaries) + '</ul><br>'
            complete_abstract = '<br><br>'.join(abstracts)

            combined_data = warning_text + f'Data as of {date_str}, <a href="{url}" target="_blank">Click here for up-to-date classifications on this variant.</a><br>' + complete_summary + complete_abstract

            writer.writerow([transcript, c_hgvs, classification, combined_data, variant_details])

        return out.getvalue()

    def filename(self) -> str:
        return self.generate_filename(suffix='mvl', extension='tsv')
