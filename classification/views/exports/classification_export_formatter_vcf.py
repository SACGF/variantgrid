from dataclasses import dataclass
from enum import Enum
from functools import cached_property
from typing import Optional

from django.conf import settings
from django.http import HttpRequest
from django.urls import reverse

from classification.enums import SpecialEKeys, AlleleOriginBucket
from classification.models import EvidenceKeyMap
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from classification.views.exports.vcf_export_utils import ExportVCF, export_vcf_info_cell, VCFHeaderType, \
    VCFHeaderNumberSpecial, VCFHeader, VCFExportTweak
from library.django_utils import get_url_from_view_path
from library.utils import local_date_str_no_dash
from ontology.models import OntologyTerm
from snpdb.models import Variant, Allele


class VCFTargetSystem(str, Enum):
    EMEDGENE = "emedgene"
    GENERIC = "generic"
    VARSEQ = "varseq"

    CORE = "core"  # should be used for everyone


@dataclass(frozen=True)
class FormatDetailsVCF:
    target_system: VCFTargetSystem = VCFTargetSystem.GENERIC

    @staticmethod
    def from_request(request: HttpRequest) -> 'FormatDetailsVCF':
        target_system_str = request.GET.get("target_system")
        target_system: VCFTargetSystem = VCFTargetSystem.GENERIC
        try:
            target_system = VCFTargetSystem(target_system_str)
        except ValueError:
            pass
        return FormatDetailsVCF(target_system=target_system)


ALLELE_ORIGIN_BUCKET_TO_LABEL = {
    AlleleOriginBucket.SOMATIC: "somatic",
    AlleleOriginBucket.GERMLINE: "germline",
    AlleleOriginBucket.UNKNOWN: "allele_origin_unknown"
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
    include_id: bool = True
    link_extra: str = ""

    def get_variant(self) -> Variant:
        return self.allele_data.cached_variant

    def get_allele(self) -> Allele:
        return self.allele_data.cached_allele

    def get_variant_id(self) -> str:
        if self.include_id:
            return super().get_variant_id()
        else:
            # some systems don't know what to do with a VariantGrid variant ID
            return "."

    @export_vcf_info_cell(
        header_id="SVLEN",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.Integer,
        description="Difference in length between REF and ALT alleles",
        categories={"system": VCFTargetSystem.CORE})
    def svlen(self):
        return self.allele_data.cached_variant.svlen

    # @export_vcf_info_cell(
    #     header_id="id",
    #     number=1,
    #     header_type=VCFHeaderType.String,
    #     description="ID in info for Emedgene",
    #     categories={"system": VCFTargetSystem.EMEDGENE})
    # def info_id(self):
    #     return super().get_variant_id()

    @export_vcf_info_cell(
        header_id="link",
        number=1,
        header_type=VCFHeaderType.String,
        description="$site_name URL",
        categories={"system": VCFTargetSystem.CORE})
    def link(self):
        return get_url_from_view_path(reverse("view_allele_compact", kwargs={"allele_id": self.allele_data.cached_allele.id})) + self.link_extra

    @export_vcf_info_cell(
        header_id="labs",
        header_type=VCFHeaderType.String,
        description="Allele origin bucket, values will be 1 or more of: somatic, germline, allele_origin_unknown",
        categories={"system": {VCFTargetSystem.GENERIC, VCFTargetSystem.VARSEQ}})
    def lab_names(self):
        all_labs = sorted({cm.classification.lab for cm in self.allele_data.cms})
        return [f"{lab.organization.shortest_name}__{lab.name}" for lab in all_labs]

    @export_vcf_info_cell(
        header_id="significance",
        number=1,
        header_type=VCFHeaderType.Integer,
        description="Most Pathogenic Classification, 1 = Benign, 2 = Likely Benign, 3 = VUS/Other, 4 = Likely Pathogenic/Oncogenic, 5 = Pathogenic/Oncogenic",
        categories={"system": VCFTargetSystem.EMEDGENE})
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
        header_id="germline",
        header_type=VCFHeaderType.Flag,
        description="Allele origin bucket, values will be 1 or more of: somatic, germline, allele_origin_unknown",
        categories={"system": VCFTargetSystem.GENERIC})
    def allele_origin_germline(self):
        if AlleleOriginBucket.GERMLINE in self.allele_data.cms_allele_origins:
            return True

    @export_vcf_info_cell(
        header_id="somatic",
        header_type=VCFHeaderType.Flag,
        description="Allele origin bucket, values will be 1 or more of: somatic, germline, allele_origin_unknown",
        categories={"system": VCFTargetSystem.GENERIC})
    def allele_origin_somatic(self):
        if AlleleOriginBucket.SOMATIC in self.allele_data.cms_allele_origins:
            return True

    @export_vcf_info_cell(
        header_id="allele_origin_unknown",
        header_type=VCFHeaderType.Flag,
        description="Allele origin bucket, values will be  1 or more of: somatic, germline, allele_origin_unknown",
        categories={"system": VCFTargetSystem.GENERIC})
    def allele_origin_unknown(self):
        if AlleleOriginBucket.UNKNOWN in self.allele_data.cms_allele_origins:
            return True

    def unique_values(self, evidence_key: str, raw: bool = False):
        e_key = EvidenceKeyMap.cached_key(evidence_key)
        unique_values = set()
        for cms in self.allele_data.cms:
            value = cms.get(evidence_key)
            if value is not None:
                unique_values.add(value)
        sorted_values = list(e_key.sort_values(unique_values))
        if raw:
            return sorted_values
        else:
            return [e_key.pretty_value(val) for val in sorted_values]

    @export_vcf_info_cell(
        header_id="condition_terms",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.String,
        description="All unique conditions this variant is being curated against, only provides values that could be matched to exact ontology terms",
        categories={"system": {VCFTargetSystem.GENERIC, VCFTargetSystem.VARSEQ}}
    )
    def condition_terms(self):
        all_terms: set[OntologyTerm] = set()
        for cm in self.allele_data.cms:
            if cr := cm.classification.condition_resolution_obj:
                all_terms.update(cr.terms)
        if all_terms:
            sorted_terms = sorted(all_terms)
            return [f"{t.id} {t.name}" for t in sorted_terms]

    @export_vcf_info_cell(
        header_id="classification",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.String,
        description="All unique classification values for this allele, B, LB, VUS, VUS_A, VUS_B, VUS_C, LP, P, LO O, D, R, B=Benign, P=Pathogenic, O=Oncogenic, D=Drug Response, R=Risk Factor",
        categories={"system": VCFTargetSystem.GENERIC}
    )
    def classification_values(self):
        return self.unique_values(SpecialEKeys.CLINICAL_SIGNIFICANCE, raw=True)

    @export_vcf_info_cell(
        header_id="classification",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.String,
        description="All unique classification values for this allele, Benign, Likely Benign, VUS, VUS A, VUS B, VUS C, Likely Pathogenic, Pathogenic, Likely Oncogenic, Oncogenic, Drug Response, Risk Factor",
        categories={"system": VCFTargetSystem.VARSEQ}
    )
    def classification_labels(self):
        return self.unique_values(SpecialEKeys.CLINICAL_SIGNIFICANCE, raw=False)

    @export_vcf_info_cell(
        header_id="somatic_clinical_significance",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.String,
        description="All unique somatic clinical significance values for this allele, values include tier_1, tier_2, tier_1_or_2, tier_3, tier_4",
        categories={"system": VCFTargetSystem.GENERIC}
    )
    def somatic_clinical_significance_values(self):
        return self.unique_values(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE, raw=True)

    @export_vcf_info_cell(
        header_id="somatic_clinical_significance",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.String,
        description="All unique somatic clinical significance values for this allele, values include Tier I, Tier II, Tier I / II, Tier III, Tier IV",
        categories={"system": VCFTargetSystem.VARSEQ}
    )
    def somatic_clinical_significance_labels(self):
        return self.unique_values(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE, raw=False)

    @property
    def _discordance_value(self) -> Optional[str]:
        for cms in self.allele_data.cms:
            if discordance_status := self.allele_data.source.is_discordant(cms):
                return discordance_status.value

    @export_vcf_info_cell(
        header_id="lab_count",
        number=1,
        header_type=VCFHeaderType.Integer,
        description="Number of unique labs with classifications for this allele",
        categories={"system": VCFTargetSystem.GENERIC}
    )
    def lab_count(self):
        return len(set(cms.lab for cms in self.allele_data.cms))

    @export_vcf_info_cell(
        header_id="discordance_status",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.String,
        description="If present, describes the discordance status of allele, if absent, the allele is not in discordance",
        categories={"system": VCFTargetSystem.GENERIC}
    )
    def discordance_status(self):
        return self._discordance_value

    @export_vcf_info_cell(
        header_id="category",
        number=1,
        header_type=VCFHeaderType.String,
        description="Free text containing various information such as discordance status, number of labs, full classification",
        categories={"system": VCFTargetSystem.EMEDGENE}
    )
    def category(self):
        """
        Jam everything into this field as it gets displayed with no processing in emedgene
        """
        parts = []
        if discordance_value := self._discordance_value:
            parts.append(f"discordance_{discordance_value}")
        parts.append(f"labs:{self.lab_count()}")

        if allele_origins := self.allele_data.cms_allele_origins:
            for allele_origin in sorted(ALLELE_ORIGIN_BUCKET_TO_LABEL.get(ao) for ao in allele_origins):
                parts.append(allele_origin)
        parts += self.unique_values(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        parts += self.unique_values(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
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
            categories={
                "system": {self.target_system.value, VCFTargetSystem.CORE}
            }
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

    @property
    def target_system(self) -> VCFTargetSystem:
        return self.format_details.target_system

    @cached_property
    def link_extra(self) -> str:
        return f"?seen={local_date_str_no_dash()}"

    def row(self, allele_data: AlleleData) -> list[str]:
        if allele_data.cms:
            if row := ClassificationVCF(
                    allele_data,
                    include_id=True,  # self.target_system != VCFTargetSystem.EMEDGENE,
                    link_extra=self.link_extra
            ).vcf_row(export_tweak=self.export_tweak):
                return [row]

    def content_type(self) -> str:
        return "text/plain"

    def extension(self) -> str:
        return "vcf"
