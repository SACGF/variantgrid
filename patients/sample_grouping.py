""" Which Samples a grouping object (Extraction) reaches.

    One implementation, used by the analysis grouping node (SampleNode above sample level), its
    canvas badge and the samples preview endpoint - if the node and the editor could disagree about
    how many VCFs are in play the badge would stop being worth glancing at.

    Keyed on the object rather than a node, since the preview has to answer before a node is saved.
"""
from dataclasses import dataclass, field
from typing import Optional

from django.contrib.auth.models import User

from patients.models import Extraction
from snpdb.archive import DataArchivedError
from snpdb.models import (
    CohortGenotypeCollection,
    CohortGenotypeStats,
    GenomeBuild,
    ImportStatus,
    Sample,
)


@dataclass(frozen=True)
class ExcludedSample:
    """ A sample the user can see, that the group leaves out """
    sample: Sample
    reason: str


@dataclass
class SampleGroup:
    """ Samples resolved for a grouping object, and what was left out.

        Exclusions are reported rather than silently dropped - the user didn't hand pick these
        samples, so narrowing they can't see is the failure mode that bites here. """
    samples: list[Sample] = field(default_factory=list)
    excluded: list[ExcludedSample] = field(default_factory=list)
    # Samples on the object the user has no permission for - counted, not named
    hidden_count: int = 0

    @property
    def vcfs(self) -> list:
        vcfs = []
        for sample in self.samples:
            if sample.vcf not in vcfs:
                vcfs.append(sample.vcf)
        return vcfs

    @property
    def warnings(self) -> list[str]:
        warnings = []
        for excluded in self.excluded:
            warnings.append(f"Sample '{excluded.sample.name}' not included: {excluded.reason}")
        if self.hidden_count:
            warnings.append(f"{self.hidden_count} sample(s) not included: no permission to view")
        return warnings


def get_extraction_sample_group(user: User, extraction: Extraction,
                                genome_build: Optional[GenomeBuild] = None) -> SampleGroup:
    """ genome_build restricts to an analysis's build - an analysis has exactly one, so a sample from
        another build is a message rather than a decision """
    group = SampleGroup()
    if extraction is None:
        return group

    all_sample_ids = set(Sample.objects.filter(extraction=extraction).values_list("pk", flat=True))
    sample_qs = Sample.filter_for_user(user).filter(extraction=extraction).select_related("vcf")
    for sample in sample_qs.order_by("pk"):
        all_sample_ids.discard(sample.pk)
        if sample.import_status != ImportStatus.SUCCESS:
            group.excluded.append(ExcludedSample(sample, "import has not completed"))
        elif sample.data_archived:
            group.excluded.append(ExcludedSample(sample, f"VCF '{sample.vcf}' data is archived"))
        elif genome_build and sample.genome_build != genome_build:
            group.excluded.append(ExcludedSample(sample, f"genome build {sample.genome_build} (need {genome_build})"))
        else:
            group.samples.append(sample)

    group.hidden_count = len(all_sample_ids)
    return group


def get_sample_variant_count(sample: Sample) -> Optional[int]:
    """ From the stats rows, so this is a join with no variant table access.
        None where stats haven't been calculated yet - the UI says so rather than showing a zero """
    try:
        cgc = sample.cohort_genotype_collection
    except (CohortGenotypeCollection.DoesNotExist, DataArchivedError):
        return None
    stats = CohortGenotypeStats.objects.filter(cohort_genotype_collection=cgc, sample=sample,
                                               filter_key__isnull=True, passing_filter=False).first()
    return stats.variant_count if stats else None


def sample_group_as_json(group: SampleGroup) -> dict:
    samples = []
    for sample in group.samples:
        samples.append({
            "sample_id": sample.pk,
            "sample": sample.name,
            "vcf_id": sample.vcf_id,
            "vcf": str(sample.vcf),
            "genome_build": str(sample.genome_build),
            "has_genotype": sample.has_genotype,
            "variant_count": get_sample_variant_count(sample),
        })

    return {
        "samples": samples,
        "excluded": [{"sample": e.sample.name, "reason": e.reason} for e in group.excluded],
        "hidden_count": group.hidden_count,
        "warnings": group.warnings,
        "totals": {
            "samples": len(group.samples),
            "vcfs": len(group.vcfs),
            "variant_count": sum(s["variant_count"] or 0 for s in samples) or None,
        },
    }
