from dataclasses import dataclass
from functools import cached_property
from typing import Iterable, Iterator

from django.urls import reverse

from classification.enums import SpecialEKeys
from classification.models import EvidenceKeyMap
from classification.views.exports.classification_export_formatter_vcf import VCFTargetSystem, \
    CLASSIFICATION_VALUE_TO_NUMBER
from classification.views.exports.vcf_export_utils import ExportVCF, export_vcf_info_cell, VCFHeaderNumberSpecial, \
    VCFHeaderType, VCFHeader
from classification.views.exports_grouping.classification_grouping_export_filter import \
    ClassificationGroupingExportFormat, ClassificationGroupingExportFormatProperties, ClassificationGroupingByAllele, \
    ClassificationGroupingExportFilter
from library.django_utils import get_url_from_view_path
from library.utils import ExportTweak
from snpdb.models import Variant, GenomeBuild, Contig, GenomeBuildContig


@dataclass(frozen=True)
class VCFFormatDetails:
    genome_build: GenomeBuild
    target_system: VCFTargetSystem


@dataclass(frozen=True)
class VCFRow(ExportVCF):
    entry: ClassificationGroupingByAllele
    date_str: str

    def get_variant(self) -> Variant:
        return self.entry.variant

    def unique_values(self, evidence_key: str, raw: bool = False):
        e_key = EvidenceKeyMap.cached_key(evidence_key)
        unique_values = set()
        for cgs in self.entry.classification_groupings:
            cms = cgs.latest_classification_modification
            value = cms.get(evidence_key)
            if value is not None:
                unique_values.add(value)
        sorted_values = list(e_key.sort_values(unique_values))
        if raw:
            return sorted_values
        else:
            return [e_key.pretty_value(val) for val in sorted_values]


    @export_vcf_info_cell(
        header_id="SVLEN",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.Integer,
        description="Difference in length between REF and ALT alleles",
        categories={"system": VCFTargetSystem.CORE})
    def svlen(self):
        return self.entry.variant.svlen

    @export_vcf_info_cell(
        header_id="significance",
        number=1,
        header_type=VCFHeaderType.Integer,
        description="Most Pathogenic Classification, 1 = Benign, 2 = Likely Benign, 3 = VUS/Other, 4 = Likely Pathogenic/Oncogenic, 5 = Pathogenic/Oncogenic",
        categories={"system": VCFTargetSystem.EMEDGENE})
    def significance(self):
        all_values = []
        for cg in self.entry.classification_groupings:
            classification_value = CLASSIFICATION_VALUE_TO_NUMBER.get(cg.latest_classification_modification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE), 3)
            all_values.append(classification_value)
        if all_values:
            return max(all_values)
        else:
            return 3

    @export_vcf_info_cell(
        header_id="link",
        number=1,
        header_type=VCFHeaderType.String,
        description="$site_name URL",
        categories={"system": VCFTargetSystem.GENERIC})
    def link(self):
        return get_url_from_view_path(reverse("view_allele_compact", kwargs={"allele_id": self.entry.allele_id})) + f"?seen={self.date_str}"


    @export_vcf_info_cell(
        header_id="classification",
        number=VCFHeaderNumberSpecial.UNBOUND,
        header_type=VCFHeaderType.String,
        description="All unique classification values for this allele, B, LB, VUS, VUS_A, VUS_B, VUS_C, LP, P, LO, O, D, R, B=Benign, P=Pathogenic, O=Oncogenic, D=Drug Response, R=Risk Factor",
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


class ClassificationGroupingExportFormatterVCF(ClassificationGroupingExportFormat):

    def __init__(self,
                 classification_grouping_filter: ClassificationGroupingExportFilter,
                 vcf_format_details: VCFFormatDetails):
        self.vcf_format_details = vcf_format_details
        super().__init__(classification_grouping_filter)

    @cached_property
    def export_tweaks(self):
        return ExportTweak(
            categories={"system": (VCFTargetSystem.CORE, self.vcf_format_details.target_system)},
        )

    @cached_property
    def genome_build(self):
        return self.vcf_format_details.genome_build

    @classmethod
    def format_properties(cls) -> ClassificationGroupingExportFormatProperties:
        return ClassificationGroupingExportFormatProperties(
            http_content_type="text/csv",
            extension="vcf"
        )

    def extra_filename_parts(self) -> list[str]:
        return [self.vcf_format_details.target_system.value, self.genome_build.name]

    def header(self) -> list[str]:

        # TODO reduce contigs to only those that are used
        # GenomeBuildContig.objects.filter(genome_build=genome_build, contig__in=contigs).order_by('order').select_related('contig')
        genome_build_str = '37'
        if self.genome_build.is_version(38):
            genome_build_str = '38'
        # self.queryset(genome_build=self.genome_build)

        # contig_field = f'latest_allele_info__grch{genome_build_str}__variant__locus__contig'
        # self.queryset(genome_build=self.genome_build).values_list(contig_field, flat=True).order_by(contig_field)

        contig_field = f'latest_allele_info__grch{genome_build_str}__variant__locus__contig'
        all_used_contig_ids = set(self.queryset(genome_build=self.genome_build).select_related(contig_field).values_list(contig_field, flat=True).order_by(contig_field).distinct())

        gb_contigs = GenomeBuildContig.objects.filter(genome_build=self.genome_build).order_by('order').select_related('contig')
        contigs = [gb_contig.contig for gb_contig in gb_contigs if gb_contig.contig_id in all_used_contig_ids]
        """
        contig_field = f'classification__allele_info__grch{genome_build_str}__variant__locus__contig'
        contig_qs = self.classification_filter.cms_qs.values_list(contig_field, flat=True).order_by(contig_field).distinct()

        if self.target_system == VCFTargetSystem.CLASSIFICATION_DUMP:
            all_vcf_headers = ClassificationSimpleVCF.complete_header(
                genome_build=self.genome_build,
                contigs=contig_qs,
                extras=[VCFHeader(settings.SITE_NAME + "_dataFilters", value=self.classification_filter.description)],
                export_tweak=self.export_tweak
            )

        :return:
        """
        return VCFRow.complete_header(
            genome_build=self.genome_build,
            contigs=contigs,
            export_tweak=self.export_tweaks,
            extras=[VCFHeader("Exported For", self.vcf_format_details.target_system.value)]
        )

    def single_row_generator(self) -> Iterator[str]:
        date_str = self.classification_grouping_filter.date_str
        for entry in self.allele_group_iterator(genome_build=self.genome_build):
            if row := VCFRow(entry=entry, date_str=date_str).vcf_row(export_tweak=self.export_tweaks):
                yield row
