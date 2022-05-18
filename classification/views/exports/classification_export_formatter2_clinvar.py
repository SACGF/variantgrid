from enum import Enum
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


class ClinVarCompareValue(int, Enum):
    # the order is as such as discordant and confidence both show a distinct difference to the only value provided by ClinVar - so you're the one causing waves
    # unknown means that while you disagree with someone, you most likely agree with someone else
    # overlap, agreement means that you agree completely with someone and we know it
    UNCLASSIFIED = 0
    DISCORDANT = 1
    CONFIDENCE = 2
    UNKNOWN = 3
    OVERLAP = 4
    AGREEMENT = 5
    NOVEL = 6
    NOT_CALCULATED = 7

    @property
    def label(self) -> str:
        return {
            ClinVarCompareValue.UNCLASSIFIED: "unclassified",
            ClinVarCompareValue.DISCORDANT: "discordant",
            ClinVarCompareValue.CONFIDENCE: "confidence",
            ClinVarCompareValue.UNKNOWN: "unknown",
            ClinVarCompareValue.OVERLAP: "overlap",
            ClinVarCompareValue.AGREEMENT: "agreement",
            ClinVarCompareValue.NOVEL: "novel",
            ClinVarCompareValue.NOT_CALCULATED: "error"
        }[self.value]

    def __str__(self):
        return f"{self.value} {self.label}"


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
        "drug_response": "D",
        "Conflicting_interpretations_of_pathogenicity": "Conflicting",
        "Conflicting_interpretations_of_pathogenicity|_association": "Conflicting"
    }
    CLINSIG_BUCKETS = {
        "B": 1,
        "LB": 1,
        "VUS": 2,
        "LP": 3,
        "P": 3,
        "R": 4,
        "D": 5
    }

    def __init__(self, allele_group: AlleleData, clinvar_version: ClinVarVersion):
        self.allele_group = allele_group
        self.clinvar_version = clinvar_version

    @export_column("Allele URL")
    def allele_url(self) -> str:
        if allele_id := self.allele_group.allele_id:
            return get_url_from_view_path(reverse('view_allele', kwargs={"allele_id": allele_id}))

    @export_column("ClinVar URL")
    def clinvar_url(self) -> str:
        clinvar: ClinVar
        if clinvar := self.clinvar:
            return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar.clinvar_variation_id}/"

    @export_column("c.HGVS")
    def c_hgvs(self) -> str:
        for cm in self.allele_group.cms:
            if c_hgvs := cm.c_hgvs_best(self.allele_group.genome_build):
                if full_c_hgvs := c_hgvs.full_c_hgvs:
                    return full_c_hgvs

    @export_column("Gene Symbol")
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
                    if part.startswith("_"):
                        part = part[1:]
                    m_parts = ClinVarCompareRow.CLINSIG_TO_VGSIG.get(part, part)
                    if isinstance(m_parts, str):
                        cs_set.add(m_parts)
                    else:
                        for m_part in m_parts:
                            cs_set.add(m_part)
        return cs_set

    @export_column("Clinical Significances")
    def servers_clinical_significance(self) -> str:
        (cs_list := list(self.server_clinical_significance_set)).sort()
        return ";".join(cs_list)

    @export_column("ClinVar Clinical Significances")
    def clinvar_clinical_significance(self) -> str:
        (cs_list := list(self.clinvar_clinical_significance_set)).sort()
        return ";".join(cs_list)

    def clinvar_clinical_significance_raw(self) -> str:
        if clinvar := self.clinvar:
            return clinvar.clinical_significance
        return ""

    @export_column("ClinVar Origin")
    def clinvar_origin(self) -> str:
        clinvar: ClinVar
        if clinvar := self.clinvar:
            return clinvar.get_origin_display()

    @export_column("ClinVar Stars")
    def clinvar_stars(self) -> str:
        clinvar: ClinVar
        if clinvar := self.clinvar:
            return clinvar.stars

    @export_column("Comparison")
    def diff_value(self) -> ClinVarCompareValue:
        our_clins: Set[str] = self.server_clinical_significance_set
        clinvar_clins: Set[str] = self.clinvar_clinical_significance_set
        clinvar_clins.discard("not_provided")

        if not clinvar_clins:
            return ClinVarCompareValue.NOVEL
        if not our_clins:
            return ClinVarCompareValue.UNCLASSIFIED

        if "Conflicting" in clinvar_clins:
            return ClinVarCompareValue.UNKNOWN

        # simplify VUS_A/B/C to just VUS
        our_clins = {"VUS" if sc.startswith("VUS") else sc for sc in our_clins}

        lowest_value = ClinVarCompareValue.NOT_CALCULATED
        clingen_buckets = {ClinVarCompareRow.CLINSIG_BUCKETS.get(cs, -1) for cs in clinvar_clins}

        def compare_single(server_clin: str) -> ClinVarCompareValue:
            nonlocal clinvar_clins
            nonlocal clingen_buckets
            server_bucket = ClinVarCompareRow.CLINSIG_BUCKETS.get(server_clin, 0)
            if server_clin in clinvar_clins:
                return ClinVarCompareValue.AGREEMENT
            elif server_bucket in clingen_buckets:
                return ClinVarCompareValue.CONFIDENCE
            else:
                return ClinVarCompareValue.DISCORDANT

        compare_values = [compare_single(sc) for sc in our_clins]
        lowest_value = min(compare_values)

        if lowest_value == ClinVarCompareValue.AGREEMENT and len(clinvar_clins) > 1:
            return ClinVarCompareValue.OVERLAP
        else:
            return lowest_value

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
        if allele_data.cms and allele_data.allele_id:
            return [ExportFormatter.write_single_row(ClinVarCompareRow(allele_data, self.clinvar_version).to_csv())]
        else:
            return list()