"""
Summaries shared by the `vg inspect` kinds: one compact dict per domain noun (variant, allele, sample,
classification, lab, user, node), so the same object looks the same whichever inspector reached it.
Each returns a `ref`-shaped dict (kind, id, label) plus the few scalars that identify the thing, never
a nested graph - the inspector for that kind is where the graph lives.
"""
from typing import Any

from django.contrib.auth.models import User

from classification.enums.classification_enums import SpecialEKeys
from classification.models.classification import Classification
from flags.models.models import FlagStatus
from library.vg.inspect import ref
from snpdb.models import Allele, Lab, Sample, Variant


def build_names(variant: Variant) -> list[str]:
    return sorted(build.name for build in variant.genome_builds)


def variant_summary(variant: Variant) -> dict[str, Any]:
    return {**ref("variant", variant, variant.full_string if hasattr(variant, "full_string") else str(variant)),
            "builds": "/".join(build_names(variant))}


def allele_summary(allele: Allele) -> dict[str, Any]:
    summary = ref("allele", allele, str(allele))
    if allele.clingen_allele_id:
        summary["clingen"] = str(allele.clingen_allele)
    return summary


def sample_summary(sample: Sample) -> dict[str, Any]:
    return {**ref("sample", sample, sample.name), "vcf": sample.vcf_id, "status": sample.get_import_status_display()}


def classification_summary(classification: Classification) -> dict[str, Any]:
    last_published = classification.last_published_version
    return {**ref("classification", classification, classification.lab_record_id or ""),
            "lab": classification.lab.group_name,
            "significance": classification.get(SpecialEKeys.CLINICAL_SIGNIFICANCE),
            "share_level": classification.share_level,
            "published": last_published is not None,
            "withdrawn": classification.withdrawn}


def lab_summary(lab: Lab) -> dict[str, Any]:
    return ref("lab", lab, lab.group_name or lab.name)


def user_summary(user: User) -> dict[str, Any]:
    return ref("user", user, user.username)


def flags_summary(obj) -> dict[str, Any] | None:
    """ Open flags on a FlagsMixin object, by type; None when it has no collection """
    collection = getattr(obj, "flag_collection", None)
    if collection is None:
        return None
    open_flags = collection.flag_set.filter(resolution__status=FlagStatus.OPEN).select_related("flag_type")
    by_type: dict[str, int] = {}
    for flag in open_flags:
        by_type[flag.flag_type_id] = by_type.get(flag.flag_type_id, 0) + 1
    return {"open": sum(by_type.values()), **by_type} if by_type else {"open": 0}
