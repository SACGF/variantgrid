from enum import Enum
from functools import cached_property
from typing import Set, Optional, List, Tuple, Callable

from django.http import HttpRequest
from django.urls import reverse

from annotation.models import ClinVar, ClinVarVersion
from classification.enums import SpecialEKeys
from classification.models import EvidenceKeyMap
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, delimited_row


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
            ClinVarCompareValue.OVERLAP: "agreement (partial)",
            ClinVarCompareValue.AGREEMENT: "agreement",
            ClinVarCompareValue.NOVEL: "novel",
            ClinVarCompareValue.NOT_CALCULATED: "error"
        }[self]

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
        "drug_response": "D"
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
        for cm in self.allele_group.cms_regardless_of_issues:
            if c_hgvs := cm.c_hgvs_best(self.allele_group.genome_build):
                if full_c_hgvs := c_hgvs.full_c_hgvs:
                    return full_c_hgvs

    @export_column("Gene Symbol")
    def gene_symbol(self) -> str:
        for cm in self.allele_group.cms_regardless_of_issues:
            if c_hgvs := cm.c_hgvs_best(self.allele_group.genome_build):
                if gene_symbol := c_hgvs.gene_symbol:
                    return gene_symbol

    @cached_property
    def server_clinical_significance_set(self) -> Set[str]:
        cs_set: Set[str] = set()
        for cm in self.allele_group.cms_regardless_of_issues:
            if classified := cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE):
                cs_set.add(classified)
        return cs_set

    @cached_property
    def clinvar_clinical_significance_set(self) -> Tuple[Set[str], Set[str]]:
        clinvar: ClinVar
        cs_set: Set[str] = set()
        unknown_set: Set[str] = set()
        if clinvar := self.clinvar:
            if cs := clinvar.clinical_significance:
                cs = cs.replace(",", "|")
                for part in cs.split("|"):
                    part = part.strip()
                    if part.startswith("_"):
                        part = part[1:]
                    if m_parts := ClinVarCompareRow.CLINSIG_TO_VGSIG.get(part):
                        if isinstance(m_parts, str):
                            cs_set.add(m_parts)
                        else:
                            for m_part in m_parts:
                                cs_set.add(m_part)
                    else:
                        unknown_set.add(part)
        return cs_set, unknown_set

    @export_column("Clinical Significances")
    def servers_clinical_significance(self) -> str:
        (cs_list := list(self.server_clinical_significance_set)).sort()
        return ";".join(cs_list)

    @export_column("ClinVar Clinical Significances")
    def clinvar_clinical_significance(self) -> str:
        known, unknown = self.clinvar_clinical_significance_set
        resolved_values = sorted(known)
        if unknown:
            resolved_values.append("?")

        if resolved_values:
            return ";".join(resolved_values)

    @export_column("ClinVar Clinical Significances Raw")
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
        clinvar_clins: Set[str]
        clinvar_unknown_clins: Set[str]
        clinvar_clins, clinvar_unknown_clins = self.clinvar_clinical_significance_set

        if (not clinvar_clins) and (not clinvar_unknown_clins):
            return ClinVarCompareValue.NOVEL
        if (not clinvar_clins) and clinvar_unknown_clins:
            return ClinVarCompareValue.UNKNOWN
        if not our_clins:
            return ClinVarCompareValue.UNCLASSIFIED

        has_unknown = any(clin_sig.startswith("?") for clin_sig in clinvar_clins)

        # simplify VUS_A/B/C to just VUS
        our_clins = {"VUS" if sc.startswith("VUS") else sc for sc in our_clins}

        lowest_value = ClinVarCompareValue.NOT_CALCULATED
        clingen_buckets = {EvidenceKeyMap.clinical_significance_to_bucket().get(cs, -1) for cs in clinvar_clins}

        def compare_single(server_clin: str) -> ClinVarCompareValue:
            nonlocal clinvar_clins
            nonlocal clingen_buckets
            server_bucket = EvidenceKeyMap.clinical_significance_to_bucket().get(server_clin, 0)
            if server_clin in clinvar_clins:
                return ClinVarCompareValue.AGREEMENT
            elif server_bucket in clingen_buckets:
                return ClinVarCompareValue.CONFIDENCE
            else:
                return ClinVarCompareValue.DISCORDANT

        compare_values = [compare_single(sc) for sc in our_clins]
        lowest_value = min(compare_values)

        if lowest_value == ClinVarCompareValue.AGREEMENT and (len(clinvar_clins) + len(clinvar_unknown_clins)) > 1:
            return ClinVarCompareValue.OVERLAP
        else:
            if clinvar_unknown_clins:
                return ClinVarCompareValue.UNKNOWN
            return lowest_value

    @export_column("$site_name Resolution Issues")
    def issues(self) -> str:
        if issues := self.allele_group.issues:
            return "\n".join(sorted(set(issue.message for issue in issues)))

    @cached_property
    def clinvar(self) -> Optional[ClinVar]:
        return self.allele_group["clinvar"]


@register_classification_exporter("clinvar_compare")
class ClassificationExportFormatterClinVarCompare(ClassificationExportFormatter):

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatter2CSV':
        return ClassificationExportFormatterClinVarCompare(
            classification_filter=ClassificationFilter.from_request(request)
        )

    def batch_pre_cache(self) -> Optional[Callable[[List[AlleleData]], None]]:
        # do we want to try all clinvar versions?
        def handle_batch(batch: List[AlleleData]):
            variant_to_batches = {}
            for ad in batch:
                if variant := ad.variant:
                    variant_to_batches[variant.pk] = ad
            clinvars = ClinVar.objects.filter(variant__in=variant_to_batches.keys(), version=self.clinvar_version)
            for clinvar in clinvars:
                if ad := variant_to_batches.get(clinvar.variant_id):
                    ad["clinvar"] = clinvar
        return handle_batch

    @cached_property
    def clinvar_version(self) -> ClinVarVersion:
        return ClinVarVersion.objects.filter(genome_build=self.classification_filter.genome_build).order_by('-created').first()

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> List[str]:
        return [delimited_row(ClinVarCompareRow.csv_header())]

    def row(self, allele_data: AlleleData) -> List[str]:
        if allele_data.allele_id:
            return [delimited_row(ClinVarCompareRow(allele_data, self.clinvar_version).to_csv())]
        else:
            return []
