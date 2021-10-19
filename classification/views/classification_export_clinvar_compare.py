from typing import Set, Optional, Iterable

from django.urls import reverse
from lazy import lazy

from annotation.models import ClinVar, AnnotationVersion, ClinVarVersion
from classification.enums import SpecialEKeys
from classification.views.classification_export_utils import ExportFormatter, AlleleGroup
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column


class ClinVarCompareRow(ExportRow):

    CLINSIG_TO_VGSIG = {
        "Benign": "B",
        "Likely_benign": "LB",
        "Uncertain_significance": "VUS",
        "Likely_pathogenic": "LP",
        "Pathogenic": "P",
        "Pathogenic/Likely_pathogenic": ["P", "LP"],
        "Benign/Likely_benign": ["B", "LB"]
    }
    CLINSIG_BUCKETS = {
        "B": 1,
        "LB": 1,
        "VUS": 2,
        "LP": 3,
        "P": 3
    }

    def __init__(self, allele_group: AlleleGroup, clinvar_version: ClinVarVersion):
        self.allele_group = allele_group
        self.clinvar_version = clinvar_version

    @export_column()
    def allele_url(self) -> str:
        return get_url_from_view_path(reverse('view_allele', kwargs={"pk": self.allele_group.allele_id}))

    @lazy
    def server_clinical_significance_set(self) -> Set[str]:
        cs_set: Set[str] = set()
        for cm in self.allele_group.data:
            cs_set.add(cm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))
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
    def clinvar_clinical_significance_raw(self) -> str:
        if clinvar := self.clinvar:
            return clinvar.clinical_significance
        return ""

    @export_column()
    def server_clinical_significance(self) -> str:
        (cs_list := list(self.server_clinical_significance_set)).sort()
        return ";".join(cs_list)

    @export_column()
    def clinvar_clinical_significance(self) -> str:
        (cs_list := list(self.clinvar_clinical_significance_set)).sort()
        return ";".join(cs_list)

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
            return "novel"
        if server_clins == clinvar_clins:
            return "same"
        if server_clins.intersection(clinvar_clins):
            return "overlap"

        def to_buckets(cs_set: Iterable[str]) -> Set[int]:
            bucket_set = set()
            for cs in cs_set:
                if bucket := ClinVarCompareRow.CLINSIG_BUCKETS.get(cs):
                    bucket_set.add(bucket)
            return bucket_set

        server_buckets = to_buckets(server_clins)
        clinsig_buckets = to_buckets(clinvar_clins)
        if server_clins.intersection(clinsig_buckets):
            return "confidence"

        return "discordant"

    @lazy
    def clinvar(self) -> Optional[ClinVar]:
        # do we want to try all clinvar versions?
        return ClinVar.objects.filter(variant__in=self.allele_group.variant_ids, version=self.clinvar_version).first()


class ExportFormatterClinVarCompare(ExportFormatter):

    @lazy
    def clinvar_version(self) -> ClinVarVersion:
        return ClinVarVersion.objects.filter(genome_build=self.genome_build).order_by('-created').first()

    def header(self) -> Optional[str]:
        return ExportFormatter.write_single_row(ClinVarCompareRow.csv_header())

    def row(self, group) -> Optional[str]:
        return ExportFormatter.write_single_row(ClinVarCompareRow(group, self.clinvar_version).to_csv())

    def filename(self) -> str:
        return self.generate_filename(suffix="clinvar_compare")
