import csv
from enum import StrEnum
from functools import cached_property
from typing import Optional
from django.conf import settings
from django.http import HttpRequest
from django.template.loader import render_to_string
from django.urls import reverse
from classification.enums import SpecialEKeys, AlleleOriginBucket
from classification.models import EvidenceKeyMap
from classification.models.evidence_mixin import SomaticClinicalSignificanceValue
from classification.views.classification_export_view import InvalidExportParameter
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from classification.views.exports.classification_export_utils import CHGVSData
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column
from ontology.models import OntologyTerm
from snpdb.models import GenomeBuild, AlleleOriginFilterDefault


class TestingContextMode(StrEnum):
    GERMLINE = "GERMLINE"
    SOMATIC = "SOMATIC"
    ONCOGENIC = "ONCOGENIC"


class FranklinExportRow(ExportRow):

    def __init__(self, data: CHGVSData, mode: TestingContextMode):
        self.data = data
        self.mode = mode

    @cached_property
    def c_hgvs_obj(self):
        return self.data.chgvs

    @export_column("Gene")
    def gene(self):
        return self.c_hgvs_obj.gene

    @export_column("transcript")
    def transcript(self):
        return self.c_hgvs_obj.transcript

    @export_column("DNA")
    def dna(self):
        return self.c_hgvs_obj.raw_c

    # don't need to provide Chromosome, Position, Ref or Alt when gene, transcript and DNA are provided

    @export_column("Classification Date")
    def classification_date(self):
        max_date = max(cm.curated_date_check for cm in self.data.cms)
        return max_date.relevant_date.date_str

    GERMLINE_CLASSIFICATION_MAPPER = {
        "B": "BENIGN",
        "LB": "LIKELY PATHOGENIC",
        "VUS": "VUS",
        "VUS_A": "VUS",
        "VUS_B": "VUS",
        "VUS_C": "VUS",
        "LP": "LIKELY PATHOGENIC",
        "P": "PATHOGENIC"
    }
    ONCOGENIC_CLASSIFICATION_MAPPER = {
        "B": "BENIGN",
        "LB": "LIKELY PATHOGENIC",
        "VUS": "VUS",
        "VUS_A": "VUS",
        "VUS_B": "VUS",
        "VUS_C": "VUS",
        # TODO should we map LP, P to Likely Oncogenic, Oncogenic too?
        "LO": "LIKELY ONCOGENIC",
        "O": "ONCOGENIC"
    }
    SOMATIC_CLIN_SIG_MAPPER = {
        SomaticClinicalSignificanceValue("tier_1"): "TIER_1",
        SomaticClinicalSignificanceValue("tier_1", "A"): "TIER_1A",
        SomaticClinicalSignificanceValue("tier_1", "B"): "TIER_1B",
        SomaticClinicalSignificanceValue("tier_1_or_2"): "TIER_1",  # TODO make sure there is a warning about this
        SomaticClinicalSignificanceValue("tier_2"): "TIER_2",
        SomaticClinicalSignificanceValue("tier_2", "C"): "TIER_2C",
        SomaticClinicalSignificanceValue("tier_2", "D"): "TIER_2D",
        SomaticClinicalSignificanceValue("tier_3"): "TIER_3",
        SomaticClinicalSignificanceValue("tier_4"): "TIER_4",
    }

    @export_column("Classification")
    def classification(self) -> Optional[str]:
        # return most pathogenic/oncogenic/important tier value (based on mode)
        if self.mode in {TestingContextMode.GERMLINE, TestingContextMode.ONCOGENIC}:
            all_values = [cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE) for cm in self.data.cms]
            all_values = [v for v in all_values if v is not None]  # clear out unclassified
            for max_value in reversed(all_values):
                if self.mode == TestingContextMode.GERMLINE:
                    if translated := FranklinExportRow.GERMLINE_CLASSIFICATION_MAPPER.get(max_value):
                        return translated
                elif self.mode == TestingContextMode.ONCOGENIC:
                    if translated := FranklinExportRow.ONCOGENIC_CLASSIFICATION_MAPPER.get(max_value):
                        return translated
        elif self.mode == TestingContextMode.SOMATIC:
            all_values = [cm.somatic_clinical_significance_value for cm in self.data.cms]
            all_values = [v for v in all_values if v is not None]  # clear out unclassified
            for max_value in reversed(all_values):
                if translated := FranklinExportRow.SOMATIC_CLIN_SIG_MAPPER.get(max_value):
                    return translated
        return None

    @export_column("Classification System")
    def classification_system(self):
        return str(self.mode)

    @export_column("Conditions")
    def conditions_column(self):
        return f"{settings.SITE_NAME} {self.mode[0]}{self.mode[1:].lower()}"

    def conditions(self) -> list[str]:
        # Don't make condition a column, as if it changes it'll go into a new section in Franklin
        condition_set: set[OntologyTerm] = set()
        plain_texts: set[str] = set()
        for cm in self.data.cms:
            condition_obj = cm.condition_resolution_obj_fallback
            if terms := condition_obj.terms:
                condition_set.update(terms)
            elif plain_text := condition_obj.plain_text:
                plain_texts.add(plain_text)

        return [str(term) for term in sorted(condition_set)] + list(sorted(plain_texts))

    @export_column("Interpretation Text")
    def summary(self):
        partial_url = reverse('view_allele', kwargs={"allele_id": self.data.allele.allele_id})
        allele_url = f"{get_url_from_view_path(partial_url)}"

        all_classification_values = [cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE) for cm in self.data.cms]
        all_classification_values = [v for v in all_classification_values if v is not None]  # clear out unclassified
        classification_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        all_classification_values = list(classification_key.sort_values(set(all_classification_values)))
        formatted_classification_values = [classification_key.pretty_value(v) for v in all_classification_values]

        all_clinsig_values = [cm.somatic_clinical_significance_value for cm in self.data.cms]
        all_clinsig_values = [v for v in all_clinsig_values if v is not None]  # clear out unclassified
        all_clinsig_values = [cs.pretty_str for cs in sorted(set(all_clinsig_values))]

        latest_date = None
        latest_date_str = None

        all_labs = set()
        for cms in self.data.cms:
            if not latest_date or latest_date < cms.curated_date_check:
                latest_date = cms.curated_date_check
            all_labs.add(cms.classification.lab)

        if latest_date:
            latest_date_str = latest_date.relevant_date.date_str

        summary = render_to_string('classification/snippets/franklin_export_summary_cell.html', {
            "labs": list(sorted(all_labs)),
            "latest_date_str": latest_date_str,
            "conditions": self.conditions(),
            "allele_url": allele_url,
            "site_name": settings.SITE_NAME,
            "date_str": self.data.source.date_str,
            "all_classification_values": formatted_classification_values,
            "all_clinsig_values": all_clinsig_values
        })
        summary = summary.replace('\n', '').replace('\t', ' ')
        # remove new lines and tabs as to not break the tsv file (Franklin doesn't like tabs to be quoted)
        return summary

    @export_column("Submitter")
    def submitter(self):
        return settings.SITE_NAME

    @export_column("Genome Build")
    def genome_build(self):
        genome_build = self.data.source.genome_build
        if genome_build == GenomeBuild.grch38():
            return "hg38"
        elif genome_build == GenomeBuild.grch37():
            return "hg19"
        else:
            raise ValueError(f"Unsupported genome build for Franklin {genome_build}")

    # Inheritance is a bit messy, not including in export
    # it also might subdivide the rows


@register_classification_exporter("franklin")
class ClassificationExportFormatterFranklin(ClassificationExportFormatter):
    """
    Exports data in the format that Franklin TSV upload expects
    """

    def __init__(self, classification_filter: ClassificationFilter):
        super().__init__(classification_filter=classification_filter)

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterFranklin':
        classification_filter = ClassificationFilter.from_request(request)
        if classification_filter.allele_origin_filter != AlleleOriginFilterDefault.GERMLINE:
            raise InvalidExportParameter("Fraknlin export requires Allele Origin to be set to Germline.")
        return ClassificationExportFormatterFranklin(
            classification_filter=classification_filter
        )

    def header(self) -> list[str]:
        return ["\t".join(FranklinExportRow.csv_header())]

    def row(self, allele_data: AlleleData) -> list[str]:
        if not allele_data.cms:
            return []
        c_data = CHGVSData(
            allele=allele_data,
            chgvs=allele_data.cms[0].c_hgvs_best(self.genome_build),
            different_chgvs=False,
            cms=allele_data.cms
        )

        def data_iterator():
            for by_allele_origin in c_data.split_by_allele_origin():
                if by_allele_origin.allele_origin == AlleleOriginBucket.GERMLINE:
                    yield FranklinExportRow(by_allele_origin, TestingContextMode.GERMLINE)
                elif by_allele_origin.allele_origin == AlleleOriginBucket.SOMATIC:
                    for mode in (TestingContextMode.SOMATIC, TestingContextMode.ONCOGENIC):
                        row = FranklinExportRow(by_allele_origin, mode)
                        if row.classification():
                            yield row

        return list(FranklinExportRow.csv_generator(
            data_iterator(),
            delimiter='\t',
            include_header=False,
            quoting=csv.QUOTE_MINIMAL
            # export_tweak=ExportTweak(categories={"format": "tsv"})
        ))

    def extension(self) -> str:
        return "tsv"

    def content_type(self) -> str:
        return "text/tab-separated-values"
