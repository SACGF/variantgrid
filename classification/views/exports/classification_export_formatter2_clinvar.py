from typing import Set, Optional, Iterable, List

from django.http import HttpRequest
from django.urls import reverse
from lazy import lazy

from annotation.models import ClinVar, ClinVarVersion
from classification.enums import SpecialEKeys
from classification.views.classification_export_utils import ExportFormatter, AlleleGroup
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter2 import ClassificationExportFormatter2
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column
from snpdb.models import VariantAllele


class ClinVarCompareRow(ExportRow):

    CLINSIG_TO_VGSIG = {
        "Benign": "B",
        "Likely_benign": "LB",
        "Uncertain_significance": "VUS",
        "Likely_pathogenic": "LP",
        "Pathogenic": "P",
        "Pathogenic/Likely_pathogenic": ["P", "LP"],
        "Benign/Likely_benign": ["B", "LB"],
        "risk_factor": "R",
        "Conflicting_interpretations_of_pathogenicity": "Conflicting"
    }
    CLINSIG_BUCKETS = {
        "B": 1,
        "LB": 1,
        "VUS": 2,
        "LP": 3,
        "P": 3
    }

    def __init__(self, allele_group: AlleleData, clinvar_version: ClinVarVersion):
        self.allele_group = allele_group
        self.clinvar_version = clinvar_version

    @export_column()
    def allele_url(self) -> str:
        return get_url_from_view_path(reverse('view_allele', kwargs={"allele_id": self.allele_group.allele_id}))

    @export_column()
    def clinvar_url(self) -> str:
        clinvar: ClinVar
        if clinvar := self.clinvar:
            return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar.clinvar_variation_id}/"

    @export_column()
    def c_hgvs(self) -> str:
        for cm in self.allele_group.cms:
            if c_hgvs := cm.c_hgvs_best(self.allele_group.genome_build):
                if full_c_hgvs := c_hgvs.full_c_hgvs:
                    return full_c_hgvs

    @export_column()
    def gene_symbol(self) -> str:
        for cm in self.allele_group.cms:
            if c_hgvs := cm.c_hgvs_best(self.allele_group.genome_build):
                if gene_symbol := c_hgvs.gene_symbol:
                    return gene_symbol

    @lazy
    def server_clinical_significance_set(self) -> Set[str]:
        cs_set: Set[str] = set()
        for cm in self.allele_group.cms:
            if classified := cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE):
                cs_set.add(classified)
        return cs_set

    @lazy
    def clinvar_clinical_significance_set(self) -> Set[str]:
        clinvar: ClinVar
        cs_set: Set[str] = set()
        if clinvar := self.clinvar:
            if cs := clinvar.clinical_significance:
                for part in cs.split(","):
                    part = part.strip()
                    m_parts = ClinVarCompareRow.CLINSIG_TO_VGSIG.get(part, part)
                    if isinstance(m_parts, str):
                        cs_set.add(m_parts)
                    else:
                        for m_part in m_parts:
                            cs_set.add(m_part)
        return cs_set

    @export_column()
    def servers_clinical_significance(self) -> str:
        (cs_list := list(self.server_clinical_significance_set)).sort()
        return ";".join(cs_list)

    @export_column()
    def clinvar_clinical_significance(self) -> str:
        (cs_list := list(self.clinvar_clinical_significance_set)).sort()
        return ";".join(cs_list)

    def clinvar_clinical_significance_raw(self) -> str:
        if clinvar := self.clinvar:
            return clinvar.clinical_significance
        return ""

    @export_column()
    def clinvar_origin(self) -> str:
        clinvar: ClinVar
        if clinvar := self.clinvar:
            return clinvar.get_origin_display()

    @export_column()
    def clinvar_stars(self) -> str:
        clinvar: ClinVar
        if clinvar := self.clinvar:
            return clinvar.stars

    @export_column()
    def diff(self):
        server_clins: Set[str] = self.server_clinical_significance_set
        clinvar_clins: Set[str] = self.clinvar_clinical_significance_set
        if "not_provided" in clinvar_clins:
            clinvar_clins.remove("not_provided")

        for vus_x in ('VUS_A', 'VUS_B', 'VUS_C'):
            if vus_x in server_clins:
                server_clins.remove(vus_x)
                server_clins.add('VUS')

        if not clinvar_clins:
            return "6 novel"
        if "Conflicting" in clinvar_clins:
            return "3 unknown"
        if server_clins == clinvar_clins:
            return "5 agreement"
        if server_clins.intersection(clinvar_clins):
            return "4 overlap"

        def to_buckets(cs_set: Iterable[str]) -> Set[int]:
            bucket_set = set()
            for cs in cs_set:
                if bucket := ClinVarCompareRow.CLINSIG_BUCKETS.get(cs):
                    bucket_set.add(bucket)
            return bucket_set

        server_buckets = to_buckets(server_clins)
        clinsig_buckets = to_buckets(clinvar_clins)
        if server_clins.intersection(clinsig_buckets):
            return "2 confidence"

        return "1 discordant"

    @lazy
    def clinvar(self) -> Optional[ClinVar]:
        # do we want to try all clinvar versions?
        variant_ids = list(VariantAllele.objects.filter(allele_id=self.allele_group.allele_id).values_list('variant_id', flat=True))
        return ClinVar.objects.filter(variant__in=variant_ids, version=self.clinvar_version).first()


@register_classification_exporter("clinvar_compare")
class ClassificationExportFormatter2ClinVarCompare(ClassificationExportFormatter2):

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatter2CSV':
        return ClassificationExportFormatter2ClinVarCompare(
            classification_filter=ClassificationFilter.from_request(request)
        )

    @lazy
    def clinvar_version(self) -> ClinVarVersion:
        return ClinVarVersion.objects.filter(genome_build=self.classification_filter.genome_build).order_by('-created').first()

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> List[str]:
        return [ExportFormatter.write_single_row(ClinVarCompareRow.csv_header())]

    def row(self, allele_data: AlleleData) -> List[str]:
        return [ExportFormatter.write_single_row(ClinVarCompareRow(allele_data, self.clinvar_version).to_csv())]
