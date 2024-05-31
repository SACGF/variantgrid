import logging
from datetime import datetime
from enum import Enum
from functools import cached_property
from itertools import groupby
from typing import Optional, Callable

from django.db.models import QuerySet
from django.http import HttpRequest
from django.urls import reverse

from annotation.clinvar_fetch_request import ClinVarFetchRequest
from annotation.models import ClinVar, ClinVarVersion, ClinVarRecord, ClinVarReviewStatus
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from classification.enums import SpecialEKeys, AlleleOriginBucket
from classification.models import EvidenceKeyMap, ClassificationModification
from classification.views.exports.classification_export_decorator import register_classification_exporter
from classification.views.exports.classification_export_filter import ClassificationFilter, AlleleData
from classification.views.exports.classification_export_formatter import ClassificationExportFormatter
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, delimited_row, first, ExportDataType
from snpdb.models import GenomeBuild, VariantAllele


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


class ClinVarCompareRowAbstract(ExportRow):

    def __init__(self,
                 allele_group: AlleleData,
                 clinvar: ClinVar):
        self.allele_group = allele_group
        self.clinvar = clinvar

    @cached_property
    def allele_url(self):
        if allele_id := self.allele_group.allele_id:
            return get_url_from_view_path(reverse('view_allele', kwargs={"allele_id": allele_id}))

    @cached_property
    def c_hgvs(self) -> str:
        for cm in self.allele_group.cms_regardless_of_issues:
            if c_hgvs := cm.c_hgvs_best(self.allele_group.genome_build):
                if full_c_hgvs := c_hgvs.full_c_hgvs:
                    return full_c_hgvs

    @cached_property
    def clinvar_url(self) -> str:
        clinvar: ClinVar
        if clinvar := self.clinvar:
            return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{clinvar.clinvar_variation_id}/"

    @cached_property
    def gene_symbol(self) -> str:
        for cm in self.allele_group.cms_regardless_of_issues:
            if c_hgvs := cm.c_hgvs_best(self.allele_group.genome_build):
                if gene_symbol := c_hgvs.gene_symbol:
                    return gene_symbol

    @cached_property
    def clinical_significance_set(self) -> set[str]:
        cs_set: set[str] = set()
        for cm in self.allele_group.cms_regardless_of_issues:
            if classified := cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE):
                cs_set.add(classified)
        return cs_set

    @cached_property
    def somatic_clinical_significance_set(self) -> set[str]:
        cs_set: set[str] = set()
        for cm in self.allele_group.cms_regardless_of_issues:
            if classified := cm.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE):
                cs_set.add(classified)
        return cs_set

    @cached_property
    def allele_origins(self):
        allele_origins = set()
        for cm in self.allele_group.cms_regardless_of_issues:
            allele_origins |= set(cm.get_value_list(SpecialEKeys.ALLELE_ORIGIN))
        e_key = EvidenceKeyMap.cached_key(SpecialEKeys.ALLELE_ORIGIN)

        allele_origin_values = ", ".join(e_key.pretty_value(val) for val in e_key.sort_values(allele_origins)) if allele_origins else "no-value"
        return f"{allele_origin_values} ({self.allele_group.allele_origin_bucket.label})"


class ClinVarCompareRow(ClinVarCompareRowAbstract):

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
    # not sure what other values are supported
    CLINSIG_TO_VGSSIG = {
        "Tier_II_-_Potential": "tier_2"
    }

    @export_column("Allele URL")
    def _allele_url(self) -> str:
        return self.allele_url

    @export_column("Clinvar URL")
    def _clinvar_url(self) -> str:
        return self.clinvar_url

    @export_column("c.HGVS")
    def _server_c_hgvs(self) -> str:
        return self.c_hgvs

    @export_column("Gene Symbol")
    def _gene_symbol(self) -> str:
        return self.gene_symbol

    @export_column("$site_name Resolution Issues")
    def issues(self) -> str:
        if issues := self.allele_group.issues:
            return "\n".join(sorted(set(issue.message for issue in issues)))

    @cached_property
    def clinvar_clinical_significance_set(self) -> tuple[set[str], set[str]]:
        clinvar: ClinVar
        cs_set: set[str] = set()
        unknown_set: set[str] = set()
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

    @export_column("$site_name Allele Origins")
    def _allele_origins(self):
        return self.allele_origins

    @export_column("ClinVar Allele Origins")
    def clinvar_origin(self) -> str:
        if clinvar := self.clinvar:
            return ", ".join(clinvar.allele_origins) + f" ({clinvar.allele_origin_bucket.label})"

    @export_column("$site_name Classification")
    def servers_classification(self) -> str:
        (cs_list := list(self.clinical_significance_set)).sort()
        return ";".join(cs_list)

    @export_column("ClinVar Classification")
    def clinvar_classification(self) -> str:
        known, unknown = self.clinvar_clinical_significance_set
        resolved_values = sorted(known)
        if unknown:
            resolved_values.append("?")

        if resolved_values:
            return ";".join(resolved_values)

    @export_column("$site_name Somatic Clinical Significances")
    def servers_somatic_clinical_significance(self) -> str:
        (cs_list := list(self.somatic_clinical_significance_set)).sort()
        return ";".join(cs_list)

    @export_column("ClinVar Somatic Clinical Significances")
    def clinvar_clinical_significance_raw(self) -> str:
        return "not-supported-yet"

    @export_column("ClinVar Stars")
    def clinvar_stars(self) -> str:
        clinvar: ClinVar
        if clinvar := self.clinvar:
            return clinvar.stars

    @export_column("Comparison")
    def diff_value(self) -> ClinVarCompareValue:
        our_clins: set[str] = self.clinical_significance_set
        clinvar_clins: set[str]
        clinvar_unknown_clins: set[str]
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


@register_classification_exporter("clinvar_compare")
class ClassificationExportFormatterClinVarCompare(ClassificationExportFormatter):

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterClinVarCompare':
        cf = ClassificationFilter.from_request(request)
        cf.allele_origin_split = True
        return ClassificationExportFormatterClinVarCompare(
            classification_filter=cf
        )

    def filter_clinvars(self, queryset: QuerySet[ClinVar]) -> QuerySet[ClinVar]:
        return queryset

    def batch_pre_cache(self) -> Optional[Callable[[list[AlleleData]], None]]:
        # do we want to try all clinvar versions?
        def handle_batch(batch: list[AlleleData]):
            allele_id_to_batch = {ad.allele_id: ad for ad in batch}
            variant_id_to_allele_id = {}
            all_allele_ids = set([ad.allele_id for ad in batch])
            for allele_id, variant_id in VariantAllele.objects.filter(allele_id__in=all_allele_ids).values_list(
                    'allele_id', 'variant_id'):
                variant_id_to_allele_id[variant_id] = allele_id

            expert_panel_clinvars = ClinVar.objects.filter(
                variant__in=variant_id_to_allele_id.keys(),
                version__in=self.clinvar_versions
            )
            expert_panel_clinvars = self.filter_clinvars(expert_panel_clinvars)

            for clinvar in expert_panel_clinvars:
                variant_id = clinvar.variant_id
                if allele_id := variant_id_to_allele_id.get(variant_id):
                    if ad := allele_id_to_batch.get(allele_id):
                        if clinvar.allele_origin_bucket == AlleleOriginBucket.UNKNOWN or clinvar.allele_origin_bucket == ad.allele_origin_bucket:
                            if (existing := ad["clinvar"]) and existing.stars > clinvar.stars:
                                # most likely these are just clinvar in different builds pointing to the same clinvar variation ID
                                # if one has stars and one doesn't, take the bigger stars
                                pass
                            else:
                                ad["clinvar"] = clinvar

        return handle_batch

    @cached_property
    def clinvar_versions(self) -> list[ClinVarVersion]:
        versions = []
        for genome_build in GenomeBuild.builds_with_annotation_cached():
            if clinvar_version := ClinVarVersion.objects.filter(genome_build=genome_build).order_by('-created').first():
                versions.append(clinvar_version)
        return versions

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> list[str]:
        return [delimited_row(ClinVarCompareRow.csv_header())]

    def row(self, allele_data: AlleleData) -> list[str]:
        if allele_data.allele_id:
            return [delimited_row(ClinVarCompareRow(
                allele_data,
                allele_data["clinvar"]
            ).to_csv())]
        else:
            return []

####################
### EXPERT PANEL ###
####################


class ClinVarExpertCompareRow(ClinVarCompareRowAbstract):

    def __init__(self,
                 allele_group: AlleleData,
                 cms: list[ClassificationModification],
                 clinvar: ClinVar,
                 clinvar_expert_record: ClinVarRecord):

        super().__init__(allele_group, clinvar)
        self.cms = cms
        self.clinvar_expert_record = clinvar_expert_record

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
    def _c_hgvs(self) -> str:
        return self.c_hgvs

    @export_column("Gene Symbol")
    def _gene_symbol(self) -> str:
        return self.gene_symbol

    @export_column("$site_name Resolution Issues")
    def issues(self) -> str:
        if issues := self.allele_group.issues:
            return "\n".join(sorted(set(issue.message for issue in issues)))

    @export_column("$site_name Clinical Significances")
    def _server_clinical_significance_set_labels(self):
        (cs_list := list(self.clinical_significance_set)).sort()
        return ";".join(cs_list)

    @export_column("ClinVar Expert Panel Classification")
    def clinvar_clinical_significance(self):
        return self.clinvar_expert_record.clinical_significance

    @export_column("ClinVar Expert Panel Somatic Significance")
    def clinvar_clinical_significance(self):
        return self.clinvar_expert_record.somatic_clinical_significance

    @export_column("Status")
    def status(self):
        bucket_map = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).option_dictionary_property("bucket")
        bucket_map["no known pathogenicity"] = 1

        server_cs_set = self.clinical_significance_set

        bucket_list = [bucket_map.get(cs) for cs in server_cs_set]
        buckets = {b for b in bucket_list if b is not None}  # Drug Response doesn't have a bucket for example

        if len(buckets) > 1:
            return "Internal Discordance"

        if bucket := bucket_map.get(self.clinvar_expert_record.clinical_significance):
            buckets.add(bucket)

        if len(buckets) > 1:
            return "Discordance"
        all_css = set().union(self.clinical_significance_set)
        all_css.add(self.clinvar_expert_record.clinical_significance)
        if len(all_css) > 1:
            return "Confidence"
        else:
            return "Agreement"

    @cached_property
    def server_last_evaluated_date(self):
        last_evaluateds = [cm.get(SpecialEKeys.CURATION_DATE) for cm in self.cms]
        last_evaluateds = [le for le in last_evaluateds if le]
        if dates := [datetime.strptime(le, "%Y-%m-%d") for le in last_evaluateds]:
            # do we want the minimum or the maximum date?
            earliest_date = min(dates)
            return earliest_date

        return None

    @export_column("$site_name Last Evaluated", data_type=ExportDataType.date)
    def server_last_evaluated(self):
        return self.server_last_evaluated_date

    @export_column("ClinVar Last Evaluated", data_type=ExportDataType.date)
    def clinvar_last_evaluated(self):
        return self.clinvar_expert_record.date_last_evaluated

    @export_column("ClinVar Allele Origins")
    def _clinvar_allele_origins(self) -> str:
        return (self.clinvar_expert_record.allele_origin or "no-value") + f" ({AlleleOriginBucket(self.clinvar_expert_record.allele_origin_bucket).label})"

    @export_column("$site_name Allele Origin")
    def _allele_origin(self):
        return self.allele_origins

    @export_column("Is ClinVar Newer")
    def is_clinvar_newer(self):
        if server_evaluated := self.server_last_evaluated_date:
            if clinvar_evaluated := self.clinvar_expert_record.date_last_evaluated:
                return datetime(year=clinvar_evaluated.year, month=clinvar_evaluated.month,
                                day=clinvar_evaluated.day) > server_evaluated


@register_classification_exporter("clinvar_compare_expert")
class ClassificationExportFormatterClinVarCompareExpert(ClassificationExportFormatterClinVarCompare):

    @classmethod
    def from_request(cls, request: HttpRequest) -> 'ClassificationExportFormatterClinVarCompareExpert':
        cf = ClassificationFilter.from_request(request)
        cf.allele_origin_split = True
        return ClassificationExportFormatterClinVarCompareExpert(
            classification_filter=cf,
        )

    def content_type(self) -> str:
        return "text/csv"

    def extension(self) -> str:
        return "csv"

    def header(self) -> list[str]:
        return [delimited_row(ClinVarExpertCompareRow.csv_header())]

    def filter_clinvars(self, queryset: QuerySet[ClinVar]) -> QuerySet[ClinVar]:
        return queryset.filter(clinvar_review_status__in=(ClinVarReviewStatus.REVIEWED_BY_EXPERT_PANEL[0], ClinVarReviewStatus.PRACTICE_GUIDELINE[0]))

    def row(self, allele_data: AlleleData) -> list[str]:
        rows = []
        clinvar_record: ClinVar
        if clinvar_record := allele_data["clinvar"]:
            records = ClinVarFetchRequest(
                clinvar_variation_id=clinvar_record.clinvar_variation_id,
                clinvar_versions=self.clinvar_versions
            ).fetch().records_with_min_stars(CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE)
            if records:
                if records := [rec for rec in records if rec.allele_origin_bucket in (AlleleOriginBucket.UNKNOWN, allele_data.allele_origin_bucket)]:
                    if len(records) > 1:
                        logging.warning(f"For allele {allele_data.allele_id}, {allele_data.allele_origin_bucket.label} we have {len(records)} expert panels")

                    def sort_key(cm):
                        return cm.classification.lab

                    all_cms = sorted([cm for cm in allele_data.all_cms if not cm.withdrawn], key=sort_key)
                    for lab, cms_by_lab in groupby(all_cms, key=sort_key):
                        new_rows = [delimited_row(
                            ClinVarExpertCompareRow(
                                allele_group=allele_data,
                                cms=list([cm.classification for cm in cms_by_lab]),
                                clinvar=clinvar_record,
                                clinvar_expert_record=first(records)
                            ).to_csv())
                        ]
                        rows += new_rows
                else:
                    # we have records for this allele, just not this allele origin
                    pass
            else:
                logging.warning(
                    f"Expected clinvar_variation_id {clinvar_record.clinvar_variation_id} to have an expert panel or higher")
        return rows
