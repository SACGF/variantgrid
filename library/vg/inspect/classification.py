"""
`vg inspect classification <pk>`: a lab's Classification - lab and owner, the allele it resolved to
and the c.HGVS per build, the key evidence values, share level and publish state, modification
history, clinical context and discordance, groupings, ClinVar export and open flags.
"""
from typing import Any

from classification.enums.classification_enums import SpecialEKeys
from classification.models.classification import Classification, ClassificationModification
from classification.models.classification_grouping import ClassificationGroupingEntry
from classification.models.clinvar_export_models import ClinVarExport
from classification.models.discordance_models import DiscordanceReport
from library.vg.inspect import capped
from library.vg.inspect.common import (
    allele_summary,
    flags_summary,
    lab_summary,
    user_summary,
    variant_summary,
)

EVIDENCE_KEYS = (SpecialEKeys.CLINICAL_SIGNIFICANCE, SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE, SpecialEKeys.ALLELE_ORIGIN,
                 SpecialEKeys.CONDITION, SpecialEKeys.GENE_SYMBOL, SpecialEKeys.GENOME_BUILD, SpecialEKeys.C_HGVS,
                 SpecialEKeys.G_HGVS, SpecialEKeys.CURATION_DATE)


def load(key: str) -> Classification:
    try:
        return Classification.objects.select_related("lab__organization", "user", "allele", "clinical_context", "allele_info").get(pk=int(key))
    except (ValueError, Classification.DoesNotExist) as e:
        raise LookupError(f"No Classification with pk {key!r}") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    classification = load(key)
    last_published = classification.last_published_version
    data: dict[str, Any] = {
        "id": classification.pk,
        "lab_record_id": classification.lab_record_id,
        "lab": lab_summary(classification.lab),
        "user": user_summary(classification.user),
        "created": classification.created.date().isoformat(),
        "modified": classification.modified.date().isoformat(),
        "share_level": classification.share_level,
        "published": {"last": last_published.created.date().isoformat(), "share_level": last_published.share_level} if last_published else None,
        "withdrawn": classification.withdrawn,
        "evidence": {key: _value(classification, key) for key in EVIDENCE_KEYS},
        "allele": _allele_info(classification),
        "flags": flags_summary(classification),
        "url": classification.get_absolute_url(),
    }
    if classification.clinical_context_id:
        context = classification.clinical_context
        data["clinical_context"] = {"id": context.pk, "name": context.name, "origin": context.allele_origin_bucket, "status": context.status,
                                    "discordance_reports": capped(DiscordanceReport.objects.filter(clinical_context=context).order_by("-pk"),
                                                                  lambda r: {"id": r.pk, "resolution": r.get_resolution_display(), "created": r.created.date().isoformat()})}
    if depth >= 2:
        data["modifications"] = capped(ClassificationModification.objects.filter(classification=classification).select_related("user").order_by("-created"), _modification)
        data["groupings"] = capped(ClassificationGroupingEntry.objects.filter(classification=classification).select_related("grouping__lab").order_by("pk"),
                                   lambda e: {"grouping": e.grouping_id, "lab": e.grouping.lab.group_name, "dirty": e.grouping.dirty})
        data["clinvar_exports"] = capped(ClinVarExport.objects.filter(classification_based_on__classification=classification).order_by("-pk"),
                                         lambda x: {"id": x.pk, "status": x.status, "scv": x.scv or None})
    return data


def _value(classification: Classification, key: str):
    value = classification.get(key)
    if isinstance(value, list):
        return ", ".join(str(v) for v in value)
    return value


def _allele_info(classification: Classification) -> dict[str, Any] | None:
    info: dict[str, Any] = {}
    if classification.allele_id:
        info["allele"] = allele_summary(classification.allele)
    if classification.variant_id:
        info["variant"] = variant_summary(classification.variant)
    allele_info = classification.allele_info
    if allele_info:
        info["imported"] = {"build": allele_info.imported_genome_build_patch_version_id, "c_hgvs": allele_info.imported_c_hgvs,
                            "g_hgvs": allele_info.imported_g_hgvs, "transcript": allele_info.imported_transcript}
        for build in ("grch37", "grch38"):
            resolved = getattr(allele_info, build)
            if resolved:
                info[build] = {"variant": resolved.variant_id, "c_hgvs": resolved.c_hgvs, "gene": resolved.gene_symbol_id, "transcript": resolved.transcript_version_id and str(resolved.transcript_version)}
        validation = allele_info.latest_validation
        if validation:
            info["validation"] = {"include": validation.include, "confirmed": validation.confirmed, "tags": validation.validation_tags_list if hasattr(validation, "validation_tags_list") else None}
    return info or None


def _modification(modification: ClassificationModification) -> dict[str, Any]:
    return {"id": modification.pk, "created": modification.created.strftime("%Y-%m-%d %H:%M"), "user": modification.user.username,
            "source": modification.source, "published": modification.published, "share_level": modification.share_level,
            "is_last_published": modification.is_last_published, "changed_keys": len(modification.delta or {})}
