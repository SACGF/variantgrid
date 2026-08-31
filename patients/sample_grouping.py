""" Which Samples a level of the Patient -> Specimen -> Extraction hierarchy reaches.

    One implementation, used by the analysis grouping node (SampleNode above sample level), its
    canvas chips and the node editor's tree - if the node and the editor could disagree about how
    many VCFs are in play the card would stop being worth glancing at.

    Keyed on the object rather than a node, since the editor has to answer before a node is saved.
"""
from dataclasses import dataclass, field
from typing import Optional

from django.contrib.auth.models import User
from django.db.models.query_utils import Q

from patients.models import Extraction, Patient, Specimen
from patients.models_enums import SampleSourceLevel, TissueStatus
from snpdb.archive import DataArchivedError
from snpdb.models import (
    CohortGenotypeCollection,
    CohortGenotypeStats,
    GenomeBuild,
    ImportStatus,
    Sample,
)

SOURCE_LEVEL_MODELS = {
    SampleSourceLevel.SAMPLE: Sample,
    SampleSourceLevel.EXTRACTION: Extraction,
    SampleSourceLevel.SPECIMEN: Specimen,
    SampleSourceLevel.PATIENT: Patient,
}


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


def _get_source_level_q(level: str, source) -> Q:
    """ The samples a source object reaches.

        Patient is a union because the two links are set independently - the VCF import carries
        extraction down without setting sample.patient, while the patient CSV sets patient and may
        leave extraction null. """
    match level:
        case SampleSourceLevel.SAMPLE:
            return Q(pk=source.pk)
        case SampleSourceLevel.EXTRACTION:
            return Q(extraction=source)
        case SampleSourceLevel.SPECIMEN:
            return Q(extraction__specimen=source)
        case SampleSourceLevel.PATIENT:
            return Q(patient=source) | Q(extraction__specimen__patient=source)
    raise ValueError(f"Unknown sample source level: '{level}'")


def get_exclusion_reason(sample: Sample, genome_build: Optional[GenomeBuild]) -> Optional[str]:
    """ genome_build restricts to an analysis's build - an analysis has exactly one, so a sample from
        another build is a message rather than a decision """
    if sample.import_status != ImportStatus.SUCCESS:
        return "import has not completed"
    if sample.data_archived:
        return f"VCF '{sample.vcf}' data is archived"
    if genome_build and sample.genome_build != genome_build:
        return f"genome build {sample.genome_build} (need {genome_build})"
    return None


def get_sample_group(user: User, level: str, source, genome_build: Optional[GenomeBuild] = None) -> SampleGroup:
    group = SampleGroup()
    if source is None:
        return group

    q = _get_source_level_q(level, source)
    all_sample_ids = set(Sample.objects.filter(q).values_list("pk", flat=True))
    sample_qs = Sample.filter_for_user(user).filter(q).select_related("vcf").distinct()
    for sample in sample_qs.order_by("pk"):
        all_sample_ids.discard(sample.pk)
        if reason := get_exclusion_reason(sample, genome_build):
            group.excluded.append(ExcludedSample(sample, reason))
        else:
            group.samples.append(sample)

    group.hidden_count = len(all_sample_ids)
    return group


def get_patient_for_source(level: str, source) -> Optional[Patient]:
    """ Every level resolves to a patient - it's what the node's pedigree badge and the editor tree
        are drawn from """
    if source is None:
        return None
    match level:
        case SampleSourceLevel.SAMPLE:
            if source.patient:
                return source.patient
            return source.extraction.specimen.patient if source.extraction else None
        case SampleSourceLevel.EXTRACTION:
            return source.specimen.patient
        case SampleSourceLevel.SPECIMEN:
            return source.patient
        case SampleSourceLevel.PATIENT:
            return source
    raise ValueError(f"Unknown sample source level: '{level}'")


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


def sample_as_json(sample: Sample, genome_build: Optional[GenomeBuild] = None) -> dict:
    reason = get_exclusion_reason(sample, genome_build)
    return {
        "sample_id": sample.pk,
        "sample": sample.name,
        "vcf_id": sample.vcf_id,
        "vcf": str(sample.vcf),
        "genome_build": str(sample.genome_build),
        "has_genotype": sample.has_genotype,
        "variant_count": get_sample_variant_count(sample),
        "included": reason is None,
        "reason": reason,
    }


def sample_group_as_json(group: SampleGroup) -> dict:
    samples = [sample_as_json(sample) for sample in group.samples]

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


def get_specimen_label(specimen: Specimen) -> str:
    """ Short enough for a chip - the full str() is the hover """
    if specimen.tissue:
        return specimen.tissue.name
    if specimen.tissue_status != TissueStatus.UNKNOWN:
        return specimen.get_tissue_status_display()
    return specimen.reference_id or specimen.external_pk or f"({specimen.pk})"


def get_extraction_label(extraction: Extraction) -> str:
    if extraction.nucleic_acid_source:
        return extraction.get_nucleic_acid_source_display()
    return extraction.reference_id or extraction.external_pk or f"({extraction.pk})"


def _vcf_filters_as_json(vcf) -> list[dict]:
    """ A VCF's own FILTER codes - only PASS means the same thing in every VCF, so the editor draws
        these under each sample row rather than pooling them """
    return [{"filter_id": vf.filter_id, "description": vf.description}
            for vf in vcf.vcffilter_set.order_by("filter_id")]


def _tree_samples(sample_qs, genome_build, selected, in_selection: bool) -> list[dict]:
    samples = []
    for sample in sample_qs:
        data = sample_as_json(sample, genome_build)
        data["vcf_filters"] = _vcf_filters_as_json(sample.vcf)
        data["level"] = SampleSourceLevel.SAMPLE
        data["id"] = sample.pk
        data["in_selection"] = in_selection or selected == (SampleSourceLevel.SAMPLE, sample.pk)
        samples.append(data)
    return samples


def get_patient_sample_tree(user: User, level: str, source, genome_build: Optional[GenomeBuild] = None) -> dict:
    """ The whole patient a source object belongs to, with the picked row's subtree flagged.

        The node editor draws this rather than just the resolved samples, so moving up or down a
        level - the RNA arm, the other specimen - is a click rather than another search. """
    patient = get_patient_for_source(level, source)
    selected = (level, source.pk) if source is not None else None
    visible_samples = Sample.filter_for_user(user).select_related("vcf")

    tree = {
        "selected": {"level": level, "id": source.pk, "label": str(source)} if source else None,
        "patient": None,
        "specimens": [],
        "unlinked": None,
    }

    if patient is None:
        # A sample linked to neither a patient nor an extraction still gets its own row to configure
        if level == SampleSourceLevel.SAMPLE:
            tree["unlinked"] = {
                "level": None,
                "id": None,
                "label": "Unlinked samples",
                "in_selection": False,
                "hidden_count": 0,
                "sample_count": 1,
                "samples": _tree_samples(visible_samples.filter(pk=source.pk), genome_build, selected, False),
            }
        return tree

    patient_selected = selected == (SampleSourceLevel.PATIENT, patient.pk)
    tree["patient"] = {
        "level": SampleSourceLevel.PATIENT,
        "id": patient.pk,
        "label": str(patient),
        "detail": patient.get_sex_display() if patient.sex else "",
        "in_selection": patient_selected,
    }

    for specimen in patient.specimen_set.select_related("tissue").order_by("pk"):
        specimen_selected = patient_selected or selected == (SampleSourceLevel.SPECIMEN, specimen.pk)
        specimen_data = {
            "level": SampleSourceLevel.SPECIMEN,
            "id": specimen.pk,
            "label": get_specimen_label(specimen),
            "detail": str(specimen),
            "in_selection": specimen_selected,
            "extractions": [],
        }
        for extraction in specimen.extraction_set.order_by("pk"):
            extraction_selected = specimen_selected or selected == (SampleSourceLevel.EXTRACTION, extraction.pk)
            sample_qs = visible_samples.filter(extraction=extraction).order_by("pk")
            all_count = Sample.objects.filter(extraction=extraction).count()
            samples = _tree_samples(sample_qs, genome_build, selected, extraction_selected)
            specimen_data["extractions"].append({
                "level": SampleSourceLevel.EXTRACTION,
                "id": extraction.pk,
                "label": get_extraction_label(extraction),
                "detail": str(extraction),
                "in_selection": extraction_selected,
                "hidden_count": all_count - len(samples),
                "sample_count": len(samples),
                "samples": samples,
            })
        specimen_data["sample_count"] = sum(e["sample_count"] for e in specimen_data["extractions"])
        tree["specimens"].append(specimen_data)

    # Samples the patient CSV attached straight to the patient, with no extraction to sit under
    unlinked_qs = visible_samples.filter(patient=patient, extraction__isnull=True).order_by("pk")
    if unlinked_samples := _tree_samples(unlinked_qs, genome_build, selected, patient_selected):
        tree["unlinked"] = {
            "level": None,
            "id": None,
            "label": "Unlinked samples",
            "in_selection": patient_selected,
            "hidden_count": 0,
            "sample_count": len(unlinked_samples),
            "samples": unlinked_samples,
        }
    tree["patient"]["sample_count"] = (sum(sp["sample_count"] for sp in tree["specimens"])
                                       + (tree["unlinked"]["sample_count"] if tree["unlinked"] else 0))
    return tree
