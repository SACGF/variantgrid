from dataclasses import dataclass
from typing import List
from django.http import HttpRequest
from django.template.loader import render_to_string
from django.urls import reverse
from lazy import lazy
from annotation.views import simple_citation_html
from classification.enums import SpecialEKeys
from classification.models import ClassificationModification, EvidenceKeyMap, ClassificationGroups
from classification.views.classification_export_mvl import CitationCounter
from classification.views.classification_export_utils import ConflictStrategy
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_formatter2 import ClassificationExportFormatter2
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter
from classification.views.exports.classification_export_utils import CHGVSData
from library.django_utils import get_url_from_view_path
from library.utils import delimited_row, export_column, ExportRow
from snpdb.models import Allele


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

        return format_details


@dataclass
class MVLCHGVSData:
    """
    Provides all the data needed to create a row in a MVL
    """
    data: CHGVSData

    @property
    def cms(self) -> [ClassificationModification]:
        return self.data.cms

    format: FormatDetailsMVL


@dataclass
class MVLClinicalSignificance:
    """
    This is the equivilent of the evidence key 'clinical_significance'
    :var alissa: The overall value to be reported to Alissa
    """
    alissa: str
    special: List[str]
    all: List[str]

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
                    (best_score < score and strategy == ConflictStrategy.MOST_BENIGN) or \
                    (best_score > score and strategy == ConflictStrategy.MOST_PATHOGENIC):
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

    @property
    def formatter(self):
        return self.mvl_data.format

    @property
    def data(self):
        return self.mvl_data.data

    @property
    def _cm(self) -> ClassificationModification:
        return self.data.cms[0]

    @lazy
    def classifications_values(self) -> MVLClinicalSignificance:
        return MVLClinicalSignificance.from_data(self.mvl_data)

    @export_column()
    def transcript(self):
        return self.data.chgvs.transcript

    @export_column()
    def c_nomen(self):
        return self.data.chgvs.raw_c

    @export_column()
    def classification(self):
        if self.formatter.is_shell:
            return 'VOUS'
        return self.classifications_values.alissa

    @lazy
    def variant_anchor_tag(self):
        # if we want to produce the same URLs as before for comparison, at the cost of a lot of speed
        if self.formatter.compatability_mode:
            v = Allele.objects.get(pk=self._cm.classification.allele_id).variant_for_build(self.mvl_data.data.source.genome_build)
            url = v.get_absolute_url()
            url = get_url_from_view_path(url) + f'?refer=mvl&seen={self.data.source.date_str}'
            return f'<a href="{url}" target="_blank">Click here for up-to-date classifications on this variant.</a>'
        else:
            # restore this after a comparison
            url = reverse('view_allele', kwargs={'pk': self._cm.classification.allele_id})
            url = get_url_from_view_path(url) + f'?refer=mvl&seen={self.data.source.date_str}'
            return f'<a href="{url}" target="_blank">Click here for up-to-date classifications on this variant.</a>'

    @lazy
    def groups(self) -> ClassificationGroups:
        return ClassificationGroups(classification_modifications=self.data.cms,
                                    genome_build=self.data.source.genome_build)

    def warnings(self) -> List[str]:
        warnings: List[str] = list()

        if self.data.different_chgvs:
            warnings.append('Warning <b>c.hgvs representations are different across transcript versions</b>')
        if self.classifications_values.special:
            warnings.append('Warning <b>Contains non-standard clinical significance</b>')

        discordant_count = 0
        for cms in self.data.cms:
            if self.data.source.is_discordant(cms):
                discordant_count += 1

        if discordant_count:
            if discordant_count == 1:
                warnings.append(f'Warning <b>1 record is in discordance</b>')
            else:
                warnings.append(f'Warning <b>{discordant_count} records are in discordance</b>')

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
        for citation, labs in citation_counter.ordered_references():
            has_citation = True
            references = ", ".join(labs)
            citations_html += f"<p>{simple_citation_html(citation)}<br><i>Referenced by</i>: {references}</p>"

        if not has_citation:
            citations_html += "No citations provided"
        return citations_html


    @export_column('variant information')
    def variant_information(self):
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


    @export_column('report abstract')
    def report_abstract(self):
        if self.formatter.is_shell:
            return 'This is a test'

        return self.variant_anchor_tag


@register_classification_exporter("mvl")
class ClassificationExportFormatter2MVL(ClassificationExportFormatter2):
    """
    Exports data in the format that Agilent's Alissa can import it
    """

    def __init__(self, classification_filter: ClassificationFilter, format_details: FormatDetailsMVL):
        self.format_details = format_details
        super().__init__(classification_filter=classification_filter)

    @staticmethod
    def from_request(request: HttpRequest) -> 'ClassificationExportFormatter2MVL':
        classification_filter = ClassificationFilter.from_request(request)
        classification_filter.rows_per_file = 9999

        return ClassificationExportFormatter2MVL(
            classification_filter=classification_filter,
            format_details=FormatDetailsMVL.from_request(request)
        )

    def header(self):
        return [delimited_row(MVLEntry.csv_header(), delimiter='\t')]

    def row(self, allele_data: AlleleData) -> List[str]:
        c_datas = CHGVSData.split_into_c_hgvs(allele_data, use_full=True)
        return list(MVLEntry.csv_generator(
            (MVLCHGVSData(c_data, self.format_details) for c_data in c_datas),
            delimiter='\t',
            include_header=False
        ))

    def extension(self) -> str:
        return "tsv"

    def content_type(self) -> str:
        return "text/tab-separated-values"
