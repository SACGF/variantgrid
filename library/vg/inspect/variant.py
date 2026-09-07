"""
`vg inspect variant <pk>`: a build-specific Variant and everything hanging off it - locus and alt, its
Allele and sibling Variants per build, the ACTIVE VariantAnnotation per build, classifications,
tags, samples carrying it (capped) and ClinGen / liftover state through the Allele.
"""
from typing import Any

from analysis.models.models_variant_tag import VariantTag
from annotation.models.models import VariantAnnotation, VariantAnnotationVersion
from classification.models.classification import Classification
from library.vg.inspect import capped, ref
from library.vg.inspect.common import (
    allele_summary,
    build_names,
    classification_summary,
    flags_summary,
)
from snpdb.models import AlleleLiftover, CohortGenotype, Variant


def load(key: str) -> Variant:
    try:
        return Variant.objects.select_related("locus__contig", "locus__ref", "alt").get(pk=int(key))
    except (ValueError, Variant.DoesNotExist) as e:
        raise LookupError(f"No Variant with pk {key!r}") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    variant = load(key)
    builds = build_names(variant)
    data: dict[str, Any] = {
        "id": variant.pk,
        "variant": str(variant),
        "builds": builds,
        "locus": {"contig": variant.locus.contig.name, "position": variant.locus.position, "ref": variant.locus.ref.seq},
        "alt": variant.alt.seq,
        "svlen": variant.svlen,
        "is_reference": variant.is_reference,
        "is_symbolic": variant.is_symbolic,
        "url": variant.get_absolute_url(),
    }
    allele = variant.allele
    data["allele"] = _allele(allele, variant) if allele else None
    data["annotation"] = _annotation(variant)
    if depth >= 2:
        data["classifications"] = capped(Classification.objects.filter(variant=variant).order_by("pk"), classification_summary)
        data["tags"] = capped(VariantTag.objects.filter(variant=variant).select_related("tag", "user").order_by("pk"), _tag)
        data["samples"] = _samples(variant)
        if allele:
            data["liftover"] = capped(AlleleLiftover.objects.filter(allele=allele).select_related("liftover").order_by("-pk"), _liftover)
    return data


def _allele(allele, variant) -> dict[str, Any]:
    summary = allele_summary(allele)
    summary["variants"] = {va.genome_build_id: va.variant_id for va in allele.variant_alleles().select_related("genome_build")}
    summary["flags"] = flags_summary(allele)
    return summary


def _annotation(variant: Variant) -> dict[str, Any]:
    per_build: dict[str, Any] = {}
    for genome_build in sorted(variant.genome_builds, key=lambda b: b.name):
        vav = VariantAnnotationVersion.latest(genome_build)
        if vav is None:
            per_build[genome_build.name] = {"vav": None}
            continue
        annotation = VariantAnnotation.objects.filter(version=vav, variant=variant).first()
        if annotation is None:
            per_build[genome_build.name] = {"vav": vav.pk, "annotated": False}
            continue
        per_build[genome_build.name] = {
            "vav": vav.pk, "gene": annotation.symbol, "transcript": annotation.transcript_version_id and str(annotation.transcript_version),
            "consequence": annotation.consequence, "impact": annotation.get_impact_display(), "hgvs_c": annotation.hgvs_c, "hgvs_p": annotation.hgvs_p,
            "gnomad_af": annotation.gnomad_af, "dbsnp": annotation.dbsnp_rs_id,
        }
    return per_build


def _tag(tag: VariantTag) -> dict[str, Any]:
    return {"tag": tag.tag_id, "analysis": tag.analysis_id, "user": tag.user.username if tag.user_id else None,
            "created": tag.created.date().isoformat()}


def _samples(variant: Variant) -> dict[str, Any]:
    """ Samples with a call at this variant (ref/het/hom; '.' missing and 'U' unknown are skipped), read from the
        packed CohortGenotype row of every collection holding it """
    genotypes = []
    for cohort_genotype in CohortGenotype.objects.filter(variant=variant).select_related("collection__cohort"):
        for sample_genotype in cohort_genotype.get_sample_genotypes():
            if sample_genotype.zygosity not in (".", "U", None):
                genotypes.append({**ref("sample", sample_genotype.sample, sample_genotype.sample.name),
                                  "zygosity": sample_genotype.zygosity, "cohort": cohort_genotype.collection.cohort_id})
    unique = {g["id"]: g for g in genotypes}
    return capped(list(unique.values()), lambda g: g)


def _liftover(allele_liftover: AlleleLiftover) -> dict[str, Any]:
    run = allele_liftover.liftover
    return {"run": run.pk, "to": run.genome_build_id, "tool": run.get_conversion_tool_display(),
            "status": allele_liftover.get_status_display(), "error": (allele_liftover.error_tidy() or None) if allele_liftover.error else None}
