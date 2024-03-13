import json
from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from typing import Any

from django.http import HttpRequest
from django.template.loader import render_to_string
from django.urls import reverse

from annotation.views import simple_citation_html
from classification.enums import SpecialEKeys
from classification.models import ClassificationModification, EvidenceKeyMap
from classification.models.classification_groups import ClassificationGroups, ClassificationGroupUtils
from classification.views.classification_export_utils import ConflictStrategy
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter, \
    DiscordanceReportStatus
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from classification.views.exports.classification_export_utils import CHGVSData, CitationCounter
from library.django_utils import get_url_from_view_path
from library.utils import delimited_row, export_column, ExportRow
from snpdb.models import Allele


class FormatDetailsMVLFileFormat(str, Enum):
    TSV = "tsv"
    JSON = "json"
    HTML = "html"


class FormatDetailsMVL:
    """
    Object to track how specific Alissa instance wants to format data.
    Specifically some understand the term "LIKELY_BENIGN" and some "LIKELY BENIGN"
    """

    def __init__(self):
        """
        :var classification_mapping: VariantGrid -> Alissa of ekey clinical_significance
        :var conflict_strategy: Alissa only allows a single classification for a variant, so if there's conflicting use this
        :var is_shell: If true, just used for testing of c.hgvs imports and omits all other data (safe to share around)
        """
        self.classification_mapping = {
            'B': 'BENIGN',
            'LB': 'LIKELY_BENIGN',
            'VUS': 'VOUS',
            'LP': 'LIKELY_PATHOGENIC',
            'P': 'PATHOGENIC'
        }
        self.compatability_mode = False
        self.conflict_strategy = ConflictStrategy.MOST_PATHOGENIC
        self.is_shell = False
        self.format: FormatDetailsMVLFileFormat = FormatDetailsMVLFileFormat.TSV

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

    def alissa_label_for(self, raw_classification: str) -> str:
        """
        :param raw_classification: The value of the ekey clinical_significance as known by variantgrid
        :return: the corresponding Alissa value (or Alissa's VUS for any non standard values)
        """
        return self.classification_mapping.get(raw_classification) or \
            self.classification_mapping.get('VUS')

    @staticmethod
    def from_request(request: HttpRequest) -> 'FormatDetailsMVL':
        """
        Create a MVL formatted based on the Export classification_export.html
        Expects the following
        cs_<b|lb|vus|lp|p> : override values for Alissa for corresponding classification
        mvl_detail : set to 'shell' if just testing
        conflict_strategy : most_benign / most_pathogenic
        :param request: Request from classification_export.html
        :return: A populated FormatDetailsMVL
        """
        format_details = FormatDetailsMVL()

        format_details.conflict_strategy = ConflictStrategy(request.query_params.get('conflict_strategy', ConflictStrategy.MOST_PATHOGENIC))
        for key in ['b', 'lb', 'vus', 'lp', 'p']:
            cs_label = request.query_params.get(f'cs_{key}')
            if cs_label:
                format_details.classification_mapping[key.upper()] = cs_label

        if request.query_params.get('mvl_detail', 'standard') == 'shell':
            format_details.is_shell = True

        if request.query_params.get('mode') == 'compatibility':
            format_details.compatability_mode = True

        if file_format := request.query_params.get('file_format'):
            format_details.format = FormatDetailsMVLFileFormat(file_format)
            # TODO allow the JSON extra to be populated

        return format_details


@dataclass
class MVLCHGVSData:
    """
    Provides all the data needed to create a row in a MVL
    """
    data: CHGVSData
    format: FormatDetailsMVL
    group_utils: ClassificationGroupUtils

    @property
    def cms(self) -> [ClassificationModification]:
        return self.data.cms


@dataclass
class MVLClinicalSignificance:
    """
    This is the equivilent of the evidence key 'clinical_significance'
    :var alissa: The overall value to be reported to Alissa
    """
    alissa: str
    special: list[str]
    all: list[str]

    @staticmethod
    def from_data(mvl_data: MVLCHGVSData):
        special_classifications = set()
        all_classifications = set()
        alissa_clasification: str = ""
        best_score = None
        strategy = mvl_data.format.conflict_strategy

        for cm in mvl_data.cms:
            raw_classification = cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
            label = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(raw_classification) or 'Unclassified'
            all_classifications.add(label)
            score: int
            if lookup_score := FormatDetailsMVL.RAW_SCORE.get(raw_classification):
                score = lookup_score
            else:
                score = FormatDetailsMVL.DEFAULT_SCORE
                special_classifications.add(label)
            if best_score is None or \
                    (score < best_score and strategy == ConflictStrategy.MOST_BENIGN) or \
                    (score > best_score and strategy == ConflictStrategy.MOST_PATHOGENIC):
                best_score = score
                alissa_clasification = mvl_data.format.alissa_label_for(raw_classification)

        return MVLClinicalSignificance(
            alissa=alissa_clasification,
            special=sorted(special_classifications),
            all=sorted(all_classifications)
        )


class MVLEntry(ExportRow):
    """
    Generates a single row for Alissa
    """

    def __init__(self, mvl_data: MVLCHGVSData):
        self.mvl_data = mvl_data

    # @export_column(categories={"format": "json"})
    # def id(self):
    # allele ID isn't unique, as a single Allele might still have different transcripts
    #     return self.mvl_data.data.allele.allele_id

    @property
    def formatter(self):
        return self.mvl_data.format

    @property
    def data(self):
        return self.mvl_data.data

    @property
    def _cm(self) -> ClassificationModification:
        return self.data.cms[0]

    @cached_property
    def classifications_values(self) -> MVLClinicalSignificance:
        return MVLClinicalSignificance.from_data(self.mvl_data)

    @export_column(categories={"format": {"tsv", "json"}})
    def transcript(self):
        return self.data.chgvs.transcript

    @export_column("c_nomen", categories={"format": "tsv"})
    def c_nomen_tsv(self):
        return self.data.chgvs.raw_c

    @export_column("cNomen", categories={"format": "json"})
    def c_nomen_json(self):
        return self.data.chgvs.raw_c

    @export_column(categories={"format": "json"})
    def gene(self):
        return self.data.chgvs.gene_symbol

    @export_column("lastUpdatedOn", categories={"format": "json"})
    def last_updated_on(self):
        return self.data.last_updated.isoformat()

    @export_column(categories={"format": {"tsv", "json"}})
    def classification(self):
        if self.formatter.is_shell:
            return 'VOUS'
        return self.classifications_values.alissa

    @cached_property
    def variant_anchor_tag(self):
        # if we want to produce the same URLs as before for comparison, at the cost of a lot of speed
        if self.formatter.compatability_mode:
            v = Allele.objects.get(pk=self._cm.classification.allele_id).variant_for_any_build(self.mvl_data.data.source.genome_build)
            url = v.get_absolute_url()
            url = get_url_from_view_path(url) + f'?refer=mvl&seen={self.data.source.date_str}'
            return f'<a href="{url}" target="_blank">Click here for up-to-date classifications on this variant.</a>'
        else:
            # restore this after a comparison
            url = reverse('view_allele', kwargs={'allele_id': self._cm.classification.allele_id})
            url = get_url_from_view_path(url) + f'?refer=mvl&seen={self.data.source.date_str}'
            return f'<a href="{url}" target="_blank">Click here for up-to-date classifications on this variant.</a>'

    @cached_property
    def groups(self) -> ClassificationGroups:
        return ClassificationGroups(
            classification_modifications=self.data.cms,
            genome_build=self.data.source.genome_build,
            group_utils=self.mvl_data.group_utils
        )

    def warnings(self) -> list[str]:
        warnings: list[str] = []

        if self.data.different_chgvs:
            warnings.append('Warning <b>c.hgvs representations are different across transcript versions</b>')
        if self.classifications_values.special:
            warnings.append('Warning <b>Contains non-standard clinical significance</b>')

        discordant_count = 0
        continued_discordance_count = 0
        pending_discordance_count = 0
        for cms in self.data.cms:
            if discordance_status := self.data.source.is_discordant(cms):
                if discordance_status == DiscordanceReportStatus.CONTINUED:
                    continued_discordance_count += 1
                elif discordance_status == DiscordanceReportStatus.PENDING_CONCORDANCE:
                    pending_discordance_count += 1
                else:
                    discordant_count += 1

        if discordant_count:
            warning = f'Warning <b>{discordant_count} ' + ('records are' if discordant_count > 1 else 'record is') + ' in active discordance</b>'
            warnings.append(warning)
        if continued_discordance_count:
            warning = f'Warning <b>{continued_discordance_count} ' + ('records are' if continued_discordance_count > 1 else 'record is') + ' in continued discordance</b>'
            warnings.append(warning)
        if pending_discordance_count:
            warning = f'Warning <b>{pending_discordance_count} ' + (
                'records are' if pending_discordance_count > 1 else 'record is') + ' in pending concordance (discordance resolution has been agreed but not yet applied)</b>'
            warnings.append(warning)

        if len(self.classifications_values.all) > 1:
            strength_list = ', '.join(self.classifications_values.all)
            warnings.append(
                f'Warning <b>Multiple clinical significances recorded for this transcript : {strength_list}</b>')
        return warnings

    def citations_html(self) -> str:
        citation_counter = CitationCounter()
        for group in self.groups:
            citation_counter.reference_citations(group.most_recent)
        citations_html = "<br><b>Citations Latest</b>:<br>"
        has_citation = False
        for citation, lab_names in citation_counter.ordered_references():
            has_citation = True
            references = ", ".join(lab_names)
            citations_html += f"<p>{simple_citation_html(citation)}<br><i>Referenced by</i>: {references}</p>"

        if not has_citation:
            citations_html += "No citations provided"
        return citations_html

    @export_column('variant information', categories={"format": "tsv"})
    def variant_information_tsv(self):
        return self._variant_information()

    @export_column('variantInfo', categories={"format": "json"})
    def variant_information_json(self):
        return self._variant_information()

    def _variant_information(self):
        if self.formatter.is_shell:
            return 'This is a test'
        groups_html = render_to_string('classification/classification_groups_mvl.html', {"groups": self.groups}) \
            .replace('\n', '').strip()
        # tidy up a bit of HTML whitespace to save space (do this after we complete testing)
        # return re.sub(r'\s{2,}', ' ', groups_html)

        warning_text = '<br>'.join(self.warnings())
        citations_html = self.citations_html()
        date_str = self.data.source.date_str
        combined_data = f'{warning_text}{groups_html}<p>Data as of {date_str} {self.variant_anchor_tag}</p>{citations_html}'

        return combined_data

    @export_column('report abstract', categories={"format": "tsv"})
    def report_abstract_tsv(self):
        return self._report_abstract()

    @export_column('reportAbstract', categories={"format": "json"})
    def report_abstract_json(self):
        return self._report_abstract()

    def _report_abstract(self):
        if self.formatter.is_shell:
            return 'This is a test'
        return self.variant_anchor_tag


@register_classification_exporter("mvl")
class ClassificationExportFormatterMVL(ClassificationExportFormatter):
    """
    Exports data in the format that Agilent's Alissa can import it
    """

    def __init__(self, classification_filter: ClassificationFilter, format_details: FormatDetailsMVL):
        self.format_details = format_details
        self.grouping_utils = ClassificationGroupUtils()
        super().__init__(classification_filter=classification_filter)

    @property
    def file_format(self) -> FormatDetailsMVLFileFormat:
        return self.format_details.format

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterMVL':
        classification_filter = ClassificationFilter.from_request(request)
        return ClassificationExportFormatterMVL(
            classification_filter=classification_filter,
            format_details=FormatDetailsMVL.from_request(request)
        )

    def header(self):
        # reset first row as we could be split over multiple files
        if self.file_format == FormatDetailsMVLFileFormat.TSV:
            return [delimited_row(MVLEntry.csv_header(categories={"format": "tsv"}), delimiter='\t')]
        elif self.file_format == FormatDetailsMVLFileFormat.JSON:

            def make_json_safe(obj: Any) -> str:
                if isinstance(obj, bool):
                    return "true" if obj else "false"
                return str(obj).replace("\n", "").replace("\"", "")

            # can't juse simple JSON because we don't want to close this off yet
            return ['{"molecularVariants":[']
        elif self.file_format == FormatDetailsMVLFileFormat.HTML:
            return ["<html><body>"]
        else:
            raise ValueError(f"Unexpected file format {self.format_details.format}")

    @property
    def delimiter_for_row(self):
        if self.file_format == FormatDetailsMVLFileFormat.JSON:
            return ","
        else:
            return ""

    def footer(self) -> list[str]:
        if self.file_format == FormatDetailsMVLFileFormat.JSON:
            return ["]}"]
        elif self.file_format == FormatDetailsMVLFileFormat.HTML:
            return ["</body></html>"]
        else:
            return []

    def row(self, allele_data: AlleleData) -> list[str]:
        c_datas = CHGVSData.split_into_c_hgvs(allele_data, use_full=True)

        if self.file_format == FormatDetailsMVLFileFormat.HTML:
            output = []
            for c_data in c_datas:
                mvl_data = MVLCHGVSData(
                    c_data,
                    self.format_details,
                    self.grouping_utils
                )
                mvl_entry = MVLEntry(mvl_data)

                html_data = render_to_string('classification/mvl_html_export.html', {"mvl_entry": mvl_entry})
                output.append(html_data)

            return ["".join(output)]

        elif self.file_format == FormatDetailsMVLFileFormat.TSV:
            return list(MVLEntry.csv_generator(
                (MVLCHGVSData(
                    c_data,
                    self.format_details,
                    self.grouping_utils
                ) for c_data in c_datas),
                delimiter='\t',
                include_header=False,
                categories={"format": "tsv"}
            ))
        elif self.file_format == FormatDetailsMVLFileFormat.JSON:
            output: list[str] = []
            for c_data in c_datas:
                export_row = MVLCHGVSData(
                    c_data,
                    self.format_details,
                    self.grouping_utils
                )
                row_json = json.dumps(MVLEntry(export_row).to_json(categories={"format": "json"}))
                output.append(row_json)
            return output

    def extension(self) -> str:
        if self.file_format == FormatDetailsMVLFileFormat.TSV:
            return "tsv"
        elif self.file_format == FormatDetailsMVLFileFormat.JSON:
            return "json"
        elif self.file_format == FormatDetailsMVLFileFormat.HTML:
            return "html"

    def content_type(self) -> str:
        if self.file_format == FormatDetailsMVLFileFormat.TSV:
            return "text/tab-separated-values"
        elif self.file_format == FormatDetailsMVLFileFormat.JSON:
            return "application/json"
        elif self.file_format == FormatDetailsMVLFileFormat.HTML:
            return "text/html"
