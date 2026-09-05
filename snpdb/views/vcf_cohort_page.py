"""  Shared context for the unified VCF / Cohort page (issue #933).

     view_vcf and view_cohort render the same template - a VCF's samples are fixed by the file, a custom
     or sub cohort's membership is editable, and everything else (tabs, samples table, wizards, related
     data) is the same page.
"""
from collections import defaultdict

from django.db.models import ForeignKey
from django.db.models.expressions import F
from django.urls.base import reverse

from annotation.models import AnnotationVersion
from annotation.tasks.calculate_sample_stats import enqueue_cohort_stats_recompute
from patients.models_enums import Sex
from snpdb.archive import DataArchivedError
from snpdb.forms import SampleChoiceForm, VCFChoiceForm
from snpdb.models import (
    VCF,
    Cohort,
    CohortGenotypeCollection,
    CohortGenotypeStats,
    CohortSample,
    Duo,
    Quad,
    Sample,
    Trio,
)

FAMILY_GROUP_MODELS = (Trio, Duo, Quad)


def sample_membership_rows(samples: list[Sample]) -> list[dict]:
    """ Sample rows for the membership editor - stats come from each sample's own VCF cohort, so a
        custom cohort drawing on several VCFs shows what the VCF page shows """
    stats_by_sample_id = {}
    if samples:
        stats_qs = CohortGenotypeStats.objects.filter(
            sample__in=samples, filter_key__isnull=True, passing_filter=False,
            cohort_genotype_collection__cohort__vcf=F("sample__vcf"),
            cohort_genotype_collection__cohort_version=F("cohort_genotype_collection__cohort__version"))
        stats_by_sample_id = {stats.sample_id: stats for stats in stats_qs}

    rows = []
    for sample in samples:
        stats = stats_by_sample_id.get(sample.pk)
        rows.append({
            "id": sample.pk,
            "name": sample.name,
            "url": sample.get_absolute_url(),
            "vcf_id": sample.vcf_id,
            "vcf_name": sample.vcf.name,
            "het_hom": f"{stats.het_count:,} / {stats.hom_count:,}" if stats else None,
            "sex": (stats.chrx_sex_guess if stats else Sex.UNKNOWN).label,
        })
    return rows


def cohort_zygosity_stats(cohort: Cohort, cohort_samples: list[Sample]) -> tuple[list[str], tuple, bool]:
    """ Per-sample zygosity counts shaped for showStackedBar, in the cohort's display order.

        A sub cohort's CGC is its parent's, so the rows are filtered to this cohort's own samples.
        Returns (sample_names, sample_zygosities, stats_pending) - stats_pending when a sample has no
        row yet, which also queues a recompute against the cohort that owns the CGC. """
    try:
        cgc = cohort.cohort_genotype_collection
    except (CohortGenotypeCollection.DoesNotExist, DataArchivedError):
        return [], (), False

    ss_fields = ("sample_id", "sample__name", "ref_count", "het_count", "hom_count", "unk_count")
    stats_by_sample_id = {}
    for value_dict in CohortGenotypeStats.objects.filter(cohort_genotype_collection=cgc,
                                                         sample__in=cohort_samples,
                                                         filter_key__isnull=True,
                                                         passing_filter=False).values(*ss_fields):
        stats_by_sample_id[value_dict.pop("sample_id")] = value_dict

    sample_names = []
    sample_zygosities = defaultdict(list)
    for sample in cohort_samples:
        if value_dict := stats_by_sample_id.get(sample.pk):
            sample_names.append(value_dict.pop("sample__name"))
            for zygosity, count in value_dict.items():
                sample_zygosities[zygosity].append(count)

    stats_pending = len(sample_names) < len(cohort_samples)
    if stats_pending:
        # Stats are per (CGC, annotation version), so a version bump or a fresh rebuild leaves gaps
        base_cohort = cohort.get_base_cohort()
        enqueue_cohort_stats_recompute(base_cohort, AnnotationVersion.latest(cohort.genome_build))
    return sample_names, tuple(sample_zygosities.items()), stats_pending


def _family_groups_by_sample_id(cohort: Cohort) -> dict[int, list[str]]:
    """ Duo/Trio/Quad point at CohortSample with on_delete=CASCADE, so dropping a sample takes the
        family group built on it - the save bar names what a pending removal would destroy """
    sample_id_by_cohort_sample_id = dict(cohort.cohortsample_set.values_list("pk", "sample_id"))
    family_groups = defaultdict(list)
    for model in FAMILY_GROUP_MODELS:
        member_fields = [f.name for f in model._meta.get_fields()
                         if isinstance(f, ForeignKey) and f.related_model is CohortSample]
        for family_group in model.objects.filter(cohort=cohort):
            label = f"{model.__name__} '{family_group.name}'" if family_group.name \
                else f"{model.__name__} {family_group.pk}"
            for field_name in member_fields:
                cohort_sample_id = getattr(family_group, f"{field_name}_id")
                if sample_id := sample_id_by_cohort_sample_id.get(cohort_sample_id):
                    if label not in family_groups[sample_id]:
                        family_groups[sample_id].append(label)
    return family_groups


def _membership_editor_context(cohort: Cohort, cohort_samples: list[Sample],
                               cohort_genotype_collection, has_write_permission: bool) -> dict:
    """ Config + rows the membership editor is built from. A sub cohort picks from its parent VCF's
        samples so every edit stays on the instant path; a custom cohort searches all samples. """
    config = {
        "cohort_id": cohort.pk,
        "version": cohort.version,
        "import_status": cohort.import_status,
        "has_genotype_data": cohort_genotype_collection is not None,
        "has_write_permission": has_write_permission,
        "running_celery_task": cohort_genotype_collection.celery_task if cohort_genotype_collection else None,
        "sample_rows_url": reverse("cohort_sample_rows", kwargs={"cohort_id": cohort.pk}),
        "save_url": reverse("cohort_sample_edit", kwargs={"cohort_id": cohort.pk}),
        "create_genotype_url": reverse("create_cohort_genotype", kwargs={"cohort_id": cohort.pk}),
        "candidate_samples": None,
        "family_groups": _family_groups_by_sample_id(cohort),
    }

    if parent_vcf := (cohort.parent_cohort.vcf if cohort.parent_cohort else None):
        member_ids = [sample.pk for sample in cohort_samples]
        candidates = parent_vcf.sample_set.exclude(pk__in=member_ids).select_related("vcf").order_by("pk")
        config["candidate_samples"] = sample_membership_rows(list(candidates))

    return {
        "membership_config": config,
        "cohort_sample_rows": sample_membership_rows(cohort_samples),
        "sample_form": SampleChoiceForm(genome_build=cohort.genome_build),
        "vcf_form": VCFChoiceForm(genome_build=cohort.genome_build),
    }


def vcf_cohort_page_context(cohort: Cohort, has_write_permission: bool, vcf: VCF = None) -> dict:
    """ Context shared by view_vcf and view_cohort. Pass vcf for the VCF-backed page - its samples are
        fixed by the file, so it gets the sample formset rather than the membership editor. """
    cohort_genotype_collection = None
    cohort_samples = []
    if cohort:
        try:
            cohort_genotype_collection = cohort.cohort_genotype_collection
        except (CohortGenotypeCollection.DoesNotExist, DataArchivedError):
            cohort_genotype_collection = None
        cohort_samples = [cs.sample for cs in cohort.get_cohort_samples()]

    context = {
        "base_template": "snpdb/menu/menu_data_base.html" if vcf else "snpdb/menu/menu_patients_base.html",
        "vcf": vcf,
        "cohort": cohort,
        "cohort_genotype_collection": cohort_genotype_collection,
        "has_write_permission": has_write_permission,
        "source_vcf_count": len({sample.vcf_id for sample in cohort_samples}),
        "permission_class": "snpdb.models.VCF" if vcf else "snpdb.models.Cohort",
        "permission_pk": vcf.pk if vcf else (cohort.pk if cohort else None),
        "show_stats_tab": bool(vcf) or cohort_genotype_collection is not None,
    }

    if vcf:
        context["page_title"] = vcf.name
    else:
        context["page_title"] = "Cohort"
        context.update(_membership_editor_context(cohort, cohort_samples, cohort_genotype_collection,
                                                  has_write_permission))
        sample_names, sample_zygosities, stats_pending = cohort_zygosity_stats(cohort, cohort_samples)
        context.update({
            "sample_names": sample_names,
            "sample_zygosities": sample_zygosities,
            "stats_pending": stats_pending,
        })
    return context
