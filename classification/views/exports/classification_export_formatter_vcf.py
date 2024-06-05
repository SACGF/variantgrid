from dataclasses import dataclass
from enum import Enum
from functools import cached_property

from django.conf import settings
from django.http import HttpRequest
from django.urls import reverse
from more_itertools import first

from classification.enums import SpecialEKeys, AlleleOriginBucket
from classification.models import EvidenceKeyMap
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from classification.views.exports.vcf_export_utils import ExportVCF, export_vcf_info_cell, VCFHeaderType, VCFHeaderNumberSpecial, VCFHeader, VCFExportTweak
from library.django_utils import get_url_from_view_path
from snpdb.models import Variant, Allele


class VCFTargetSystem(str, Enum):
    EMEDGENE = "emedgene"


@dataclass(frozen=True)
class FormatDetailsVCF:
    target_system: VCFTargetSystem = VCFTargetSystem.EMEDGENE

    @staticmethod
    def from_request(request: HttpRequest) -> 'FormatDetailsVCF':
        return FormatDetailsVCF()


ALLELE_ORIGIN_BUCKET_TO_LABEL = {
    AlleleOriginBucket.SOMATIC: "somatic",
    AlleleOriginBucket.GERMLINE: "germline",
    AlleleOriginBucket.UNKNOWN: "unknown"
}


CLASSIFICATION_VALUE_TO_NUMBER = {
    "B": 1,
    "LB": 2,
    "LP": 4,
    "LO": 4,
    "P": 5,
    "O": 5
}


@dataclass(frozen=True)
class ClassificationVCF(ExportVCF):
    allele_data: AlleleData

    def get_variant(self) -> Variant:
        return self.allele_data.cached_variant

    def get_allele(self) -> Allele:
        return self.allele_data.cached_allele

    @export_vcf_info_cell(
        header_id="SVLEN",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.Integer,
        description="Difference in length between REF and ALT alleles",
        categories={"standard": True})
    def svlen(self):
        return self.allele_data.cached_variant.svlen

    @export_vcf_info_cell(
        header_id="link",
        number=1,
        header_type=VCFHeaderType.String,
        description="$site_name URL",
        categories={"standard": True, VCFTargetSystem.EMEDGENE.value: True})
    def link(self):
        return get_url_from_view_path(reverse("view_allele_compact", kwargs={"allele_id": self.allele_data.cached_allele.id}))

    @export_vcf_info_cell(
        header_id="significance",
        number=1,
        header_type=VCFHeaderType.Integer,
        description="Classification, 1 = Benign, 2 = Likely Benign, 3 = VUS/Other, 4 = Likely Pathogenic/Oncogenic, 5 = Pathogenic/Oncogenic",
        categories={VCFTargetSystem.EMEDGENE.value: True})
    def significance(self):
        all_values = []
        for cm in self.allele_data.cms:
            classification_value = CLASSIFICATION_VALUE_TO_NUMBER.get(cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE), 3)
            all_values.append(classification_value)
        if all_values:
            return max(all_values)
        else:
            return 3

    @export_vcf_info_cell(
        header_id="allele_origin",
        number=1,
        header_type=VCFHeaderType.String,
        description="Allele origin bucket, values will be: somatic, germline, mixed or unknown",
        categories={"standard": True})
    def allele_origin(self):
        allele_origins = self.allele_data.cms_allele_origins
        if all_origin_counts := len(allele_origins):
            if all_origin_counts > 1:
                return "mixed"
            else:
                return ALLELE_ORIGIN_BUCKET_TO_LABEL.get(first(allele_origins))
        else:
            # not sure how we got 0 allele origins, but just in case
            return "unknown"

    @export_vcf_info_cell(
        header_id="category",
        number=1,
        header_type=VCFHeaderType.String,
        description="Free text containing various information such as discordance status, number of labs, full classification",
        categories={VCFTargetSystem.EMEDGENE.value: True}
    )
    def category(self):
        """
        Jam everything into this field as it gets displayed with no processing in emedgene
        """
        parts = []
        for cms in self.allele_data.cms:
            if discordance_status := self.allele_data.source.is_discordant(cms):
                parts.append("discordance_" + discordance_status.value)
                break
        labs = set()
        cs_values = set()
        scs_values = set()
        clin_sig_e_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        somatic_clin_sig_e_key = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        for cms in self.allele_data.cms:
            labs.add(cms.lab)
            if cs_val := cms.get(SpecialEKeys.CLINICAL_SIGNIFICANCE):
                cs_values.add(cs_val)
            if scs_val := cms.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE):
                scs_values.add(scs_val)

        parts.append(f"labs:{len(labs)}")
        parts += [clin_sig_e_key.pretty_value(cs_val) for cs_val in clin_sig_e_key.sort_values(cs_values)]
        parts += [somatic_clin_sig_e_key.pretty_value(cs_val) for cs_val in somatic_clin_sig_e_key.sort_values(scs_values)]
        return "|".join(parts)


@register_classification_exporter("vcf")
class ClassificationExportFormatterVCF(ClassificationExportFormatter):

    def __init__(self, classification_filter: ClassificationFilter, format_details: FormatDetailsVCF):
        self.format_details = format_details
        self.report_keys = [EvidenceKeyMap.cached_key(key) for key in [SpecialEKeys.CLINICAL_SIGNIFICANCE]]
        super().__init__(classification_filter=classification_filter)

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterVCF':
        classification_filter = ClassificationFilter.from_request(request)
        return ClassificationExportFormatterVCF(
            classification_filter=classification_filter,
            format_details=FormatDetailsVCF.from_request(request)
        )

    @cached_property
    def export_tweak(self) -> VCFExportTweak:
        return VCFExportTweak(
            # TODO support other formats
            categories={"emedgene": True}
        )

    def header(self) -> list[str]:
        genome_build_str = '37'
        if self.genome_build.is_version(38):
            genome_build_str = '38'
        contig_field = f'classification__allele_info__grch{genome_build_str}__variant__locus__contig'
        contig_qs = self.classification_filter.cms_qs.values_list(contig_field, flat=True).order_by(contig_field).distinct()

        all_vcf_headers = ClassificationVCF.complete_header(
            genome_build=self.genome_build,
            contigs=contig_qs,
            extras=[VCFHeader(settings.SITE_NAME + "_dataFilters", value=self.classification_filter.description)],
            export_tweak=self.export_tweak
        )
        return [str(header) for header in all_vcf_headers]

    def row(self, allele_data: AlleleData) -> list[str]:
        if allele_data.cms:
            if row := ClassificationVCF(allele_data).vcf_row(export_tweak=self.export_tweak):
                return [row]

    def content_type(self) -> str:
        return "text/plain"

    def extension(self) -> str:
        return "vcf"
