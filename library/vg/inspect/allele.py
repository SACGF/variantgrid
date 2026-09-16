"""
`vg inspect allele <pk|CA…>`: the build-independent Allele - its Variant per build, ClinGen record,
classifications, clinical contexts (discordance), liftover attempts and open flags.
"""
from typing import Any

from classification.models.classification import Classification
from classification.models.clinical_context_models import ClinicalContext
from library.vg.inspect import capped
from library.vg.inspect.common import classification_summary, flags_summary, variant_summary
from snpdb.models import Allele, AlleleLiftover, ClinGenAllele


def load(key: str) -> Allele:
    key = key.strip()
    if key.upper().startswith(ClinGenAllele.CLINGEN_ALLELE_CODE_PREFIX if hasattr(ClinGenAllele, "CLINGEN_ALLELE_CODE_PREFIX") else "CA"):
        clingen_id = ClinGenAllele.get_id_from_code(key) if hasattr(ClinGenAllele, "get_id_from_code") else int(key[2:])
        allele = Allele.objects.filter(clingen_allele_id=clingen_id).first()
        if allele is None:
            raise LookupError(f"No Allele linked to ClinGen {key}")
        return allele
    try:
        return Allele.objects.get(pk=int(key))
    except (ValueError, Allele.DoesNotExist) as e:
        raise LookupError(f"No Allele with pk {key!r} (a CA… ClinGen id is also accepted)") from e


def inspect(key: str, depth: int) -> dict[str, Any]:
    allele = load(key)
    data: dict[str, Any] = {
        "id": allele.pk,
        "allele": str(allele),
        "clingen": str(allele.clingen_allele) if allele.clingen_allele_id else None,
        "url": allele.get_absolute_url(),
        "variants": [{**variant_summary(va.variant), "build": va.genome_build_id, "origin": va.get_origin_display()}
                     for va in allele.variant_alleles().select_related("variant__locus__contig", "variant__locus__ref", "variant__alt", "genome_build")],
        "flags": flags_summary(allele),
    }
    if depth >= 2:
        data["classifications"] = capped(Classification.objects.filter(allele=allele).order_by("pk"), classification_summary)
        data["clinical_contexts"] = capped(ClinicalContext.objects.filter(allele=allele).order_by("pk"), _clinical_context)
        data["liftover"] = capped(AlleleLiftover.objects.filter(allele=allele).select_related("liftover").order_by("-pk"), _liftover)
    return data


def _clinical_context(clinical_context: ClinicalContext) -> dict[str, Any]:
    return {"id": clinical_context.pk, "name": clinical_context.name, "origin": clinical_context.allele_origin_bucket,
            "status": clinical_context.status, "classifications": clinical_context.classification_set.count()}


def _liftover(allele_liftover: AlleleLiftover) -> dict[str, Any]:
    run = allele_liftover.liftover
    return {"run": run.pk, "to": run.genome_build_id, "tool": run.get_conversion_tool_display(),
            "status": allele_liftover.get_status_display()}
