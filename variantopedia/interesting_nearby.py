from django.contrib.auth.models import User
from django.db.models import Count, Sum, Q

from annotation.annotation_version_querysets import get_variant_queryset_for_latest_annotation_version
from snpdb.models import Variant


def interesting_counts(qs, user, genome_build, clinical_significance=False):
    from classification.models import Classification

    agg_kwargs = {
        "total": Count("id"),
        "total_het": Sum("global_variant_zygosity__het_count"),
        "total_hom": Sum("global_variant_zygosity__hom_count"),
    }

    classifications = {
        "clinvar": (None, "highest_pathogenicity"),
        "classification": (Classification.get_variant_q(user, genome_build), "clinical_significance")
    }

    for classification, (classification_q, clinical_significance_path) in classifications.items():
        agg_kwargs[f"{classification}_count"] = Count("id", filter=classification_q)
        if clinical_significance:
            for cs in range(1, 6):
                q_clinical_significance = Q(**{f"{classification}__{clinical_significance_path}": cs})
                if classification_q:
                    q_clinical_significance &= classification_q
                agg_kwargs[f"{classification}_{cs}"] = Count("id", filter=q_clinical_significance)

    return qs.aggregate(**agg_kwargs)


def variant_interesting_summary(variant: Variant, user: User) -> str:
    qs = get_variant_queryset_for_latest_annotation_version(variant.genome_build)
    qs = qs.filter(pk=variant.pk)
    counts = interesting_counts(qs)

    FIELD_LABELS = {
        "clinvar_count": "ClinVar",
        "clinvar_path": ""
    }

    interesting_summary = ""

    return interesting_summary