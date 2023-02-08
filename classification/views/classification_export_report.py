from typing import Optional

from django.contrib.auth.models import User
from django.http import HttpResponse
from django.template import engines

from classification.models import ClassificationJsonParams, ClassificationReportTemplate
from classification.models.classification import ClassificationModification, \
    Classification
from classification.models.evidence_key import EvidenceKeyMap
from snpdb.models import GenomeBuild


class ClassificationReport:
    """
    Formats using report for the corresponding lab.
    Typically you'd only use for a single record
    """

    def __init__(self, classification: ClassificationModification, user: User):
        self.classification = classification
        self.user = user

    def serve(self):
        row_datas = []
        # only support 1 record for now
        vcm = self.classification
        row_data = self.row_data(vcm)
        row_datas.append(row_data)

        report = ClassificationReportTemplate.preferred_template_for(vcm)
        template_str = report.template or 'No report template has been configured'
        django_engine = engines['django']
        template = django_engine.from_string(template_str)
        content = template.render({'record': row_datas[0]})

        response = HttpResponse(content=content, content_type='text/html')
        # response['Content-Disposition'] = f'attachment; filename="{self.filename()}"'
        return response

    def row_data(self, record: ClassificationModification) -> dict:
        context = {}
        evidence = record.as_json(ClassificationJsonParams(self.user, include_data=True))['data']
        e_keys = EvidenceKeyMap.instance(lab=record.classification.lab)

        for e_key in e_keys.all_keys:
            blob = evidence.get(e_key.key) or {}

            report_blob = {}
            report_blob['value'] = blob.get('value', None)
            report_blob['note'] = blob.get('note', None)
            report_blob['formatted'] = e_key.pretty_value(blob)
            report_blob['label'] = e_key.pretty_label
            context[e_key.key] = report_blob

        for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
            c_hgvs = record.classification.get_c_hgvs(genome_build)
            key = "c_hgvs_" + genome_build.pk.lower()
            report_blob = {}
            report_blob['value'] = c_hgvs
            report_blob['note'] = None
            report_blob['formatted'] = c_hgvs
            report_blob['label'] = "c.HGVS"
            context[key] = report_blob

        context['condition_resolved'] = record.classification.condition_resolution
        context['citations'] = record.loaded_citations().to_json()
        context['evidence_weights'] = Classification.summarize_evidence_weights(evidence)
        context['acmg_criteria'] = record.criteria_strength_summary(e_keys)
        context['editable'] = record.classification.can_write(self.user)
        return context
