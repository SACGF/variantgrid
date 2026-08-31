""" Which Samples a level of the Patient -> Specimen -> Extraction hierarchy reaches.

    One implementation, used by the analysis grouping node (SampleNode above sample level), its
    canvas chips and the node editor's tree - if the node and the editor could disagree about how
    many VCFs are in play the card would stop being worth glancing at.

    Keyed on the object rather than a node, since the editor has to answer before a node is saved.
"""
import operator
from collections import defaultdict
from dataclasses import dataclass, field
from functools import reduce
from typing import Optional

from django.contrib.auth.models import User
from django.db.models import Count, F
from django.db.models.query_utils import Q

from patients.models import Extraction, Patient, Specimen
from patients.models_enums import SampleSourceLevel, TissueStatus
from snpdb.archive import DataArchivedError
from snpdb.models import (
    CohortGenotypeCollection,
    CohortGenotypeCollectionType,
    CohortGenotypeStats,
    GenomeBuild,
    ImportStatus,
    Sample,
    VCFFilter,
)

@dataclass(frozen=True)
class SourceLevel:
    """ What one level of the hierarchy is, so the levels differ in a table rather than in branches """
    model: type
    # Sample lookups reaching this object's samples, OR'd together. Patient takes two because the
    # links are set independently - the VCF import carries extraction down without setting
    # sample.patient, while the patient CSV sets patient and may leave extraction null
    sample_paths: tuple[str, ...]
    # Attribute paths from the object to its patient, first one that resolves wins.
    # Empty means the object is already the patient
    patient_paths: tuple[str, ...] = ()


SOURCE_LEVELS = {
    SampleSourceLevel.SAMPLE: SourceLevel(
        model=Sample, sample_paths=("pk",),
        patient_paths=("patient", "extraction__specimen__patient")),
    SampleSourceLevel.EXTRACTION: SourceLevel(
        model=Extraction, sample_paths=("extraction",), patient_paths=("specimen__patient",)),
    SampleSourceLevel.SPECIMEN: SourceLevel(
        model=Specimen, sample_paths=("extraction__specimen",), patient_paths=("patient",)),
    # Keep in step with Patient.get_samples(), which expresses the same union model side
    SampleSourceLevel.PATIENT: SourceLevel(
        model=Patient, sample_paths=("patient", "extraction__specimen__patient")),
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


def get_source_level_q(level: str, source) -> Q:
    """ The samples a source object reaches. Every path takes the pk, so one lookup shape covers
        'pk' and the FK chains alike """
    paths = SOURCE_LEVELS[level].sample_paths
    return reduce(operator.or_, [Q(**{path: source.pk}) for path in paths])


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

    q = get_source_level_q(level, source)
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


def _follow_path(obj, path: str):
    """ Walk an attribute path, giving up at the first null link. A name that isn't there still
        raises - that's a typo in the table above, not a sample without a patient """
    for attribute in path.split("__"):
        obj = getattr(obj, attribute)
        if obj is None:
            return None
    return obj


def get_patient_for_source(level: str, source) -> Optional[Patient]:
    """ Every level resolves to a patient - it's what the node's pedigree badge and the editor tree
        are drawn from """
    if source is None:
        return None
    patient_paths = SOURCE_LEVELS[level].patient_paths
    if not patient_paths:
        return source  # A Patient is its own
    for path in patient_paths:
        if patient := _follow_path(source, path):
            return patient
    return None


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


def sample_as_json(sample: Sample, genome_build: Optional[GenomeBuild] = None,
                   variant_counts: Optional[dict] = None) -> dict:
    """ Pass variant_counts (sample pk -> count) to skip the per sample stats lookup """
    reason = get_exclusion_reason(sample, genome_build)
    if variant_counts is None:
        variant_count = get_sample_variant_count(sample)
    else:
        variant_count = variant_counts.get(sample.pk)
    return {
        "sample_id": sample.pk,
        "sample": sample.name,
        "vcf_id": sample.vcf_id,
        "vcf": str(sample.vcf),
        "genome_build": str(sample.genome_build),
        "has_genotype": sample.has_genotype,
        "variant_count": variant_count,
        "included": reason is None,
        "reason": reason,
    }


def get_sample_variant_counts(samples: list[Sample]) -> dict[int, int]:
    """ Batched get_sample_variant_count - the editor tree asks for a whole patient at once.

        Mirrors Cohort.cohort_genotype_collection: the current version's UNCOMMON collection, and
        nothing for an archived VCF (whose partition has been dropped). """
    live = [s for s in samples if not s.data_archived]
    if not live:
        return {}
    cgc_qs = CohortGenotypeCollection.objects.filter(
        cohort__in={s.vcf.cohort.pk for s in live},  # Reverse one-to-one, so no cohort_id
        cohort_version=F("cohort__version"),
        collection_type=CohortGenotypeCollectionType.UNCOMMON)
    cgc_ids = set(cgc_qs.values_list("pk", flat=True))
    stats_qs = CohortGenotypeStats.objects.filter(cohort_genotype_collection__in=cgc_ids,
                                                  sample__in=live,
                                                  filter_key__isnull=True, passing_filter=False)
    return dict(stats_qs.values_list("sample_id", "variant_count"))


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


def _get_vcf_filters(samples: list[Sample]) -> dict[int, list[dict]]:
    """ Each VCF's own FILTER codes, keyed by VCF - only PASS means the same thing in every VCF, so
        the editor draws these under each sample row rather than pooling them. One query for the
        patient rather than one per sample, since a patient's samples share few VCFs """
    vcf_filters = defaultdict(list)
    vcf_filter_qs = VCFFilter.objects.filter(vcf__in={s.vcf_id for s in samples}).order_by("filter_id")
    for vcf_id, filter_id, description in vcf_filter_qs.values_list("vcf_id", "filter_id", "description"):
        vcf_filters[vcf_id].append({"filter_id": filter_id, "description": description})
    return vcf_filters


def _tree_samples(samples: list[Sample], genome_build, variant_counts, vcf_filters) -> list[dict]:
    rows = []
    for sample in samples:
        data = sample_as_json(sample, genome_build, variant_counts=variant_counts)
        data["vcf_filters"] = vcf_filters.get(sample.vcf_id, [])
        data["level"] = SampleSourceLevel.SAMPLE
        data["id"] = sample.pk
        rows.append(data)
    return rows


def get_patient_sample_tree(user: User, level: str, source, genome_build: Optional[GenomeBuild] = None) -> dict:
    """ The whole patient a source object belongs to.

        The node editor draws this rather than just the resolved samples, so moving up or down a
        level - the RNA arm, the other specimen - is a click rather than another search. Which rows
        the pick actually reads is worked out client side, since the tree is redrawn as the user
        moves between them without fetching again.

        `samples` are the rows with no extraction above them: the patient's own where there is a
        patient, and just the picked sample where there is not - a deployment that hasn't set up
        specimens and extractions has no hierarchy to draw, so the editor says so and shows the one
        row rather than inventing containers around it. """
    patient = get_patient_for_source(level, source)
    tree = {
        "selected": {"level": level, "id": source.pk, "label": str(source)} if source else None,
        "patient": None,
        "specimens": [],
        "samples": [],
    }

    # Everything the rows need, up front - this is per editor open, so a query per sample shows
    visible_samples = Sample.filter_for_user(user).select_related("vcf__cohort")
    if patient is None:
        samples = list(visible_samples.filter(pk=source.pk)) if level == SampleSourceLevel.SAMPLE else []
    else:
        samples = list(visible_samples.filter(get_source_level_q(SampleSourceLevel.PATIENT, patient))
                       .distinct().order_by("pk"))
    row_kwargs = {
        "genome_build": genome_build,
        "variant_counts": get_sample_variant_counts(samples),
        "vcf_filters": _get_vcf_filters(samples),
    }

    if patient is None:
        tree["samples"] = _tree_samples(samples, **row_kwargs)
        return tree

    by_extraction = defaultdict(list)
    unlinked = []
    for sample in samples:
        if sample.extraction_id:
            by_extraction[sample.extraction_id].append(sample)
        else:
            unlinked.append(sample)
    # Samples on the extraction the user can't view are a count - naming them would leak them
    all_counts = dict(Sample.objects.filter(extraction__specimen__patient=patient)
                      .values_list("extraction_id").annotate(n=Count("pk")))

    tree["patient"] = {
        "level": SampleSourceLevel.PATIENT,
        "id": patient.pk,
        "label": str(patient),
        "detail": patient.get_sex_display() if patient.sex else "",
    }
    specimen_qs = patient.specimen_set.select_related("tissue").prefetch_related("extraction_set")
    for specimen in specimen_qs.order_by("pk"):
        specimen_data = {
            "level": SampleSourceLevel.SPECIMEN,
            "id": specimen.pk,
            "label": get_specimen_label(specimen),
            "detail": str(specimen),
            "extractions": [],
        }
        for extraction in sorted(specimen.extraction_set.all(), key=lambda e: e.pk):
            extraction_samples = by_extraction.get(extraction.pk, [])
            specimen_data["extractions"].append({
                "level": SampleSourceLevel.EXTRACTION,
                "id": extraction.pk,
                "label": get_extraction_label(extraction),
                "detail": str(extraction),
                "hidden_count": all_counts.get(extraction.pk, 0) - len(extraction_samples),
                "sample_count": len(extraction_samples),
                "samples": _tree_samples(extraction_samples, **row_kwargs),
            })
        specimen_data["sample_count"] = sum(e["sample_count"] for e in specimen_data["extractions"])
        tree["specimens"].append(specimen_data)

    # Samples the patient CSV attached straight to the patient, with no extraction to sit under
    tree["samples"] = _tree_samples(unlinked, **row_kwargs)
    tree["patient"]["sample_count"] = (sum(sp["sample_count"] for sp in tree["specimens"])
                                       + len(tree["samples"]))
    return tree
