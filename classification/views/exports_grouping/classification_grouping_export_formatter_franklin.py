import csv
from dataclasses import dataclass
from enum import StrEnum
from functools import cached_property
from typing import Iterator, Optional

from django.conf import settings
from django.template.loader import render_to_string
from django.urls import reverse

from classification.enums import ClassificationResultValue, SpecialEKeys, AlleleOriginBucket
from classification.models import EvidenceKeyMap
from classification.models.evidence_mixin import SomaticClinicalSignificanceValue
from classification.views.exports_grouping.classification_grouping_export_filter import \
    ClassificationGroupingExportFormat, ClassificationGroupingExportFormatProperties, \
    ClassificationGroupingExportFilter, ClassificationGroupingByAlleleAndOrigin
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, delimited_row, export_column
from ontology.models import OntologyTerm
from snpdb.models import GenomeBuild


class TestingContextMode(StrEnum):
    GERMLINE = "GERMLINE"
    SOMATIC = "SOMATIC"
    ONCOGENIC = "ONCOGENIC"


@dataclass(frozen=True)
class FranklinFormatDetails:
    genome_build: GenomeBuild


# Copy and pasted as the old export will eventually be deleted
class FranklinExportRow(ExportRow):

    def __init__(self, data: ClassificationGroupingByAlleleAndOrigin, mode: TestingContextMode, date_str: str):
        self.data = data
        self.mode = mode
        self.date_str = date_str

    @cached_property
    def c_hgvs_obj(self):
        return self.data.c_hgvs

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
        # hardcoded to 2000-01-01 so as to not change and cause duplicate records
        # as well as be old enough that it shouldn't be suggested to copy the details into a new record
        return "2000-01-01"

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
            all_values = [grouping.latest_cached_summary_obj.pathogenicity.classification for grouping in self.data.classification_groupings]
            all_values = [v for v in all_values if v is not None]  # clear out unclassified
            for max_value in reversed(all_values):
                if self.mode == TestingContextMode.GERMLINE:
                    if translated := FranklinExportRow.GERMLINE_CLASSIFICATION_MAPPER.get(max_value):
                        return translated
                elif self.mode == TestingContextMode.ONCOGENIC:
                    if translated := FranklinExportRow.ONCOGENIC_CLASSIFICATION_MAPPER.get(max_value):
                        return translated
        elif self.mode == TestingContextMode.SOMATIC:
            all_values = [grouping.latest_cached_summary_obj.somatic.somatic_clinical_significance_value for grouping in self.data.classification_groupings]
            all_values = [v for v in all_values if v is not None]  # clear out unclassified
            for max_value in reversed(all_values):
                if translated := FranklinExportRow.SOMATIC_CLIN_SIG_MAPPER.get(max_value):
                    return translated
        return None

    @export_column("Classification System")
    def classification_system(self):
        return str(self.mode)

    # Franklin (non-free version) only allows a subset of valid conditions, so just leave condition blank
    @export_column("Conditions")
    def conditions_column(self):
        return ""
        #return f"{settings.SITE_NAME} {self.mode[0]}{self.mode[1:].lower()}"

    def conditions(self) -> list[str]:
        # Don't make condition a column, as if it changes it'll go into a new section in Franklin
        condition_set: set[OntologyTerm] = set()
        plain_texts: set[str] = set()
        for cm in self.data.classification_groupings:
            condition_obj = cm.conditions_obj
            if terms := condition_obj.terms:
                condition_set.update(terms)
            elif plain_text := condition_obj.plain_text:
                plain_texts.add(plain_text)

        return [str(term) for term in sorted(condition_set)] + list(sorted(plain_texts))

    @export_column("Interpretation Text")
    def summary(self):
        partial_url = reverse('view_allele', kwargs={"allele_id": self.data.allele_id})
        allele_url = f"{get_url_from_view_path(partial_url)}"

        classification_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        all_classification_values = [cm.latest_cached_summary_obj.pathogenicity.classification for cm in self.data.classification_groupings]
        all_classification_values = [v for v in all_classification_values if v is not None]  # clear out unclassified
        all_classification_values = list(classification_key.sort_values(set(all_classification_values)))
        formatted_classification_values = [classification_key.pretty_value(v) for v in all_classification_values]

        all_clinsig_values = [cm.latest_cached_summary_obj.somatic.somatic_clinical_significance_value for cm in self.data.classification_groupings]
        all_clinsig_values = [v for v in all_clinsig_values if v is not None]  # clear out unclassified
        all_clinsig_values = [cs.pretty_value for cs in sorted(set(all_clinsig_values))]

        latest_date = None
        latest_date_str = None

        all_labs = set()
        for cms in self.data.classification_groupings:
            if effective_date := cms.latest_cached_summary_obj.date:
                if not latest_date or latest_date < effective_date:
                    latest_date = effective_date
            all_labs.add(cms.lab)

        if latest_date:
            latest_date_str = latest_date.date

        summary = render_to_string('classification/snippets/franklin_export_summary_cell.html', {
            "labs": list(sorted(all_labs)),
            "latest_date_str": latest_date_str,
            "conditions": self.conditions(),
            "allele_url": allele_url,
            "site_name": settings.SITE_NAME,
            "date_str": self.date_str,
            "all_classification_values": formatted_classification_values,
            "all_clinsig_values": all_clinsig_values
        })
        summary = summary.replace('\n', '').replace('\t', ' ')
        # remove new lines and tabs as to not break the tsv file (Franklin doesn't like tabs to be quoted)
        return summary

    @export_column("Submitter")
    def submitter(self):
        return settings.SITE_NAME

    @export_column("Classification Tags")
    def classification_tags(self):
        # TODO make this a setting
        return "Shariant"

    @export_column("Genome Build")
    def genome_build(self):
        genome_build = self.data.genome_build
        if genome_build == GenomeBuild.grch38():
            return "hg38"
        elif genome_build == GenomeBuild.grch37():
            return "hg19"
        else:
            raise ValueError(f"Unsupported genome build for Franklin {genome_build}")


class ClassificationGroupingExportFormatterFranklin(ClassificationGroupingExportFormat):

    def __init__(self,
                 classification_grouping_filter: ClassificationGroupingExportFilter,
                 franklin_formatter_details: FranklinFormatDetails
                 ):
        self.franklin_formatter_details = franklin_formatter_details
        super().__init__(classification_grouping_filter)

    @classmethod
    def format_properties(cls) -> ClassificationGroupingExportFormatProperties:
        return ClassificationGroupingExportFormatProperties(
            http_content_type="text/tsv",
            extension="tsv"
        )

    @property
    def genome_build(self) -> GenomeBuild:
        return self.franklin_formatter_details.genome_build

    def header(self) -> list[str]:
        return [delimited_row(FranklinExportRow.csv_header(), delimiter='\t', include_new_line=False)]

    def single_row_generator(self) -> Iterator[str]:

        def data_iterator():
            for allele_grouping in self.allele_group_iterator():
                for allele_origin_grouping in allele_grouping.sub_by_allele_origin():
                    if allele_origin_grouping.allele_origin_bucket == AlleleOriginBucket.GERMLINE:
                        row = FranklinExportRow(allele_origin_grouping, TestingContextMode.GERMLINE, date_str=self.classification_grouping_filter.date_str)
                        if row.c_hgvs_obj and row.classification():
                            yield row

                    elif allele_origin_grouping.allele_origin_bucket == AlleleOriginBucket.SOMATIC:
                        for mode in (TestingContextMode.SOMATIC, TestingContextMode.ONCOGENIC):
                            row = FranklinExportRow(allele_origin_grouping, mode, date_str=self.classification_grouping_filter.date_str)
                            if row.c_hgvs_obj and row.classification():
                                yield row

        for row in data_iterator():
            yield delimited_row(row.to_csv(), delimiter='\t', include_new_line=False, quoting=csv.QUOTE_MINIMAL)
        # for row in FranklinExportRow.csv_generator(
        #     data_iterator(),
        #     delimiter='\t',
        #     include_header=False,
        #     quoting=csv.QUOTE_MINIMAL
        #     # export_tweak=ExportTweak(categories={"format": "tsv"})
        # ):
        #     print(f"** {row}")
        #     # yield row
        #     yield "jojo"
