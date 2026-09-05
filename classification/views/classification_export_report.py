import re
from itertools import groupby
from typing import Optional

from bs4 import BeautifulSoup
from django.contrib.auth.models import User
from django.http import HttpResponse
from django.template import engines

from classification.enums import SpecialEKeys
from classification.models import ClassificationJsonParams, ClassificationReportTemplate
from classification.models.classification import Classification, ClassificationModification
from classification.models.evidence_key import EvidenceKeyMap
from snpdb.models import GenomeBuild


class ClassificationReport:
    """
    Formats using report for the corresponding lab.

    One record fills "record" the way it always has. A multi-variant report (the sample/patient page's
    Classify & Report tab) passes the rest in as well - "classifications" is every record in gene order and
    "gene_groups" the same records grouped by gene symbol - so single-variant templates keep working unchanged.
    """

    def __init__(self, classification: ClassificationModification, user: User,
                 classifications: Optional[list[ClassificationModification]] = None,
                 report_template: Optional[ClassificationReportTemplate] = None,
                 extra_context: Optional[dict] = None):
        self.classification = classification
        self.classifications = classifications if classifications is not None else [classification]
        self.report_template = report_template
        self.extra_context = extra_context or {}  # eg the case header a multi-variant report opens with
        self.user = user

    def get_template(self):
        report = self.report_template or ClassificationReportTemplate.preferred_template_for(self.classification)
        template_str = report.template or 'No report template has been configured'
        django_engine = engines['django']
        return django_engine.from_string(template_str)

    @staticmethod
    def _gene_symbol(record: ClassificationModification) -> str:
        return record.get(SpecialEKeys.GENE_SYMBOL) or ""

    def context(self) -> dict:
        records = sorted(self.classifications, key=self._gene_symbol)
        row_data_by_pk = {record.pk: self.row_data(record) for record in records}
        row_data = [row_data_by_pk[record.pk] for record in records]

        gene_groups = []
        for gene_symbol, group in groupby(records, key=self._gene_symbol):
            gene_groups.append({"gene_symbol": gene_symbol,
                                "classifications": [row_data_by_pk[record.pk] for record in group]})

        record = row_data_by_pk.get(self.classification.pk) or self.row_data(self.classification)
        return {
            **self.extra_context,
            'record': record,
            'classifications': row_data,
            'gene_groups': gene_groups,
        }

    def serve(self):
        template = self.get_template()
        content = template.render(self.context())
        response = HttpResponse(content=content, content_type='text/html')
        return response

    def row_data(self, record: ClassificationModification) -> dict:
        context = {}
        evidence = record.as_json(ClassificationJsonParams(self.user, include_data=True))['data']
        e_keys = EvidenceKeyMap.instance().with_overrides(record.classification.evidence_key_overrides)

        for e_key in e_keys.all_keys:
            blob = evidence.get(e_key.key) or {}

            report_blob = {
                'value': blob.get('value', None),
                'note': blob.get('note', None),
                'formatted': e_key.pretty_value(blob),
                'label': e_key.pretty_label
            }
            # Vue/JS can't handle ":" in names
            key = e_key.key.replace(':', '_')
            context[key] = report_blob

        for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
            c_hgvs = record.classification.get_c_hgvs(genome_build)
            key = "c_hgvs_" + genome_build.pk.lower()
            report_blob = {
                'value': c_hgvs,
                'note': None,
                'formatted': c_hgvs,
                'label': "c.HGVS"
            }
            context[key] = report_blob

        context['condition_resolved'] = record.classification.condition_resolution
        context['citations'] = record.loaded_citations().to_json()
        context['evidence_weights'] = Classification.summarize_evidence_weights(evidence)
        context['acmg_criteria'] = record.criteria_strength_summary(e_keys)
        context['editable'] = record.classification.can_write(self.user)
        return context

    def get_unknown_evidence(self) -> list[str]:
        row_data = self.row_data(self.classification)
        template = self.get_template()
        content = template.render({'record': row_data})
        return self._get_unknown_evidence(content, row_data)

    @staticmethod
    def _get_unknown_evidence(content, row_data) -> list[str]:
        soup = BeautifulSoup(content, "html.parser")
        unknown_evidence = set()
        for tag in soup.find_all(attrs={":evidence": True}):
            evidence = tag[":evidence"]
            if re.match(r"[{(.]", evidence):
                continue
            if evidence not in row_data:
                unknown_evidence.add(evidence)
        return sorted(unknown_evidence)
