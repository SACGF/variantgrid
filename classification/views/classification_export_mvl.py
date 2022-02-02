import csv
import io
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Set, Any, Tuple, Iterable

from django.conf import settings
from django.template.loader import render_to_string
from django.utils.timezone import now
from pytz import timezone

from annotation.citations import get_citations, CitationDetails
from annotation.models import CitationSource, Citation
from annotation.views import simple_citation_html
from classification.enums.classification_enums import SpecialEKeys
from classification.models import ClassificationGroups, ClassificationModification
from classification.models.evidence_key import EvidenceKeyMap
from classification.views.classification_export_utils import ExportFormatter, \
    AlleleGroup, ConflictStrategy
from library.django_utils import get_url_from_view_path


@dataclass(frozen=True, eq=True)
class CitationStub:
    source: str
    idx: str

    def __lt__(self, other):
        if self.source < other.source:
            return True
        if self.source == other.source:
            return self.idx.rjust(10, '0') < other.idx.rjust(10, '0')
        return False


class CitationCounter:

    def __init__(self):
        self.all_citations: Dict[CitationStub, Set[str]] = defaultdict(set)

    def reference_citations(self, cm: ClassificationModification):
        for db_ref in cm.db_refs:
            db = db_ref.get('db')
            if source := CitationSource.CODES.get(db):
                idx = db_ref.get('idx')
                stub = CitationStub(source=source, idx=idx)
                self.all_citations[stub].add(cm.classification.lab.name)

    def ordered_references(self) -> Iterable[Tuple[CitationDetails, List[Any]]]:
        citations: List[Citation] = list()

        by_source: Dict[str, List[str]] = defaultdict(list)
        for stub in list(self.all_citations.keys()):
            by_source[stub.source].append(stub.idx)

        # bulk select and select cachedcitation since we're going to be asking for that soon
        for source, keys in by_source.items():
            found = set()
            for cit in Citation.objects.select_related('cachedcitation').filter(citation_source=source, citation_id__in=keys):
                found.add(cit.citation_id)
                citations.append(cit)

            for check in keys:
                if check not in found:
                    citations.append(Citation.objects.create(citation_source=source, citation_id=check))

        details = get_citations(citations)

        for citation_detail in details:
            stub = CitationStub(CitationSource.CODES.get(citation_detail.source), citation_detail.citation_id)
            references = list(self.all_citations[stub])
            references.sort()
            yield citation_detail, references


class ExportFormatterMVL(ExportFormatter):
    """
    Formats classifications for Agilent 5.2 MVL usage
    """

    @property
    def version(self):
        return '3'

    @property
    def use_full_chgvs(self):
        return True

    def __init__(self, conflict_strategy: str, cs_override_labels: Dict[str, str], *args, **kwargs):
        super().__init__(*args, **kwargs)
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

    def row(self, group: AlleleGroup) -> str:
        out = io.StringIO()
        writer = csv.writer(out, delimiter='\t')

        date_str = now().astimezone(tz=timezone(settings.TIME_ZONE)).strftime("%Y-%m-%d")
        url = get_url_from_view_path(group.target_variant.get_absolute_url()) + f'?refer=mvl&seen={date_str}'
        variant_details = f'<a href="{url}" target="_blank">Click here for up-to-date classifications on this variant.</a>'

        for c_parts, vcms_w_chgvs in group.iter_c_hgvs_versionless_transcripts():

            transcript = c_parts.transcript
            c_hgvs = c_parts.raw_c

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

                if vcm_w_chgvs.chgvs.raw_c != c_hgvs:
                    has_diff_chgvs = True

            groups = ClassificationGroups(classification_modifications=[cnchgvs.vcm for cnchgvs in vcms_w_chgvs], genome_build=self.genome_build)
            groups_html = render_to_string('classification/classification_groups_mvl.html', {"groups": groups}).replace('\n', '').strip()
            citation_counter = CitationCounter()
            for group in groups:
                citation_counter.reference_citations(group.most_recent)

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

            citations_html = "<br><b>Citations Latest</b>:<br>"
            has_citation = False
            for citation, labs in citation_counter.ordered_references():
                has_citation = True
                references = ", ".join(labs)
                citations_html += f"<p>{simple_citation_html(citation)}<br><i>Referenced by</i>: {references}</p>"

            if not has_citation:
                citations_html += "No citations provided"

            combined_data = f'{warning_text}{groups_html}<p>Data as of {date_str} <a href="{url}" target="_blank">Click here for up-to-date classifications on this variant.</a></p>{citations_html}'

            self.row_count += 1
            writer.writerow([transcript, c_hgvs, classification, combined_data, variant_details])

        return out.getvalue()

    def filename(self) -> str:
        return self.generate_filename(suffix='mvl', extension='tsv')
