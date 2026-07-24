""" Variant / Allele <-> Beacon `g_variant` coordinate translation.

One coordinate-mapping module, both directions (§5.1, §9.2):
- inbound: parse a Beacon g_variants query -> a permission-agnostic `Variant` queryset,
  and render a resolved `Variant` -> a Beacon `g_variant` record.
- outbound: render a `Variant` -> Beacon query params, and parse a remote Beacon's
  response -> (exists, count).

Beacon variant coordinates are 0-based (`start`); VariantGrid loci are 1-based, so we
convert on the boundary (start_0based == position - 1).
"""
from typing import Optional

from django.db.models import QuerySet

from snpdb.models import Variant, GenomeBuild


class BeaconQueryError(ValueError):
    """ A malformed Beacon g_variants request (bad params / unknown assembly). """


def genome_build_from_assembly(assembly_id: str) -> GenomeBuild:
    """ Resolve a Beacon `assemblyId` (e.g. "GRCh37", "GRCh38") to a GenomeBuild. """
    try:
        return GenomeBuild.get_name_or_alias(assembly_id)
    except (GenomeBuild.DoesNotExist, ValueError) as e:
        raise BeaconQueryError(f"Unknown assemblyId '{assembly_id}'") from e


def variant_qs_for_beacon_query(reference_name: str, start_0based: int, reference_bases: str,
                                alternate_bases: str, genome_build: GenomeBuild) -> QuerySet[Variant]:
    """ Build a `Variant` queryset for an exact Beacon g_variants coordinate query,
        restricted to the contigs of `genome_build`. Permission scoping happens later
        (per dataset); this is the raw coordinate match. """
    position = start_0based + 1  # Beacon 0-based -> VariantGrid 1-based
    qs = Variant.objects.filter(Variant.get_contigs_q(genome_build))
    qs = qs.filter(Variant.get_chrom_q(reference_name))
    qs = qs.filter(locus__position=position,
                   locus__ref__seq=reference_bases.upper(),
                   alt__seq=alternate_bases.upper())
    return qs


def variant_to_g_variant(variant: Variant, genome_build: GenomeBuild) -> dict:
    """ Render a resolved `Variant` as a Beacon `g_variant` record (record-tier detail). """
    vc = variant.coordinate  # VariantCoordinate (1-based)
    g_variant = {
        "variantInternalId": str(variant.pk),
        "variation": {
            "referenceName": vc.chrom.replace("chr", ""),
            "assemblyId": genome_build.name,
            "start": vc.position - 1,  # 1-based -> 0-based
            "referenceBases": vc.ref,
            "alternateBases": vc.alt,
        },
    }
    if variant.svlen is not None:
        g_variant["variation"]["variantType"] = variant.alt.seq  # e.g. <DEL> / <DUP>
        g_variant["variation"]["svlen"] = variant.svlen

    identifiers = []
    allele = getattr(variant, "allele", None)
    if allele and allele.clingen_allele:
        identifiers.append({"clinGenAlleleId": str(allele.clingen_allele)})
    if identifiers:
        g_variant["identifiers"] = identifiers
    return g_variant


# ----------------------------- outbound (§9.2) -----------------------------

def variant_to_beacon_query_params(variant: Variant, genome_build: GenomeBuild) -> Optional[dict]:
    """ Build the g_variants query params we POST/GET to an external Beacon. Only a public
        coordinate leaves. Returns None for symbolic variants (no explicit ref/alt to send). """
    vc = variant.coordinate
    if vc.is_symbolic or not (vc.ref and vc.alt):
        return None
    return {
        "referenceName": vc.chrom.replace("chr", ""),
        "assemblyId": genome_build.name,
        "start": vc.position - 1,  # 1-based -> 0-based
        "referenceBases": vc.ref,
        "alternateBases": vc.alt,
    }


def parse_beacon_response(data: dict) -> tuple[Optional[bool], Optional[int]]:
    """ Parse a remote Beacon g_variants response into (exists, count).
        `exists` is None when the node did not answer: Beacon v2 uses `exists: null` for
        "cannot answer", and some live servers omit responseSummary entirely. That is not
        the same as absence, so it must not collapse to False - a node that declines would
        otherwise render as a confident "not found".
        `count` is None when the node reports presence but no numeric count. """
    summary = (data or {}).get("responseSummary") or {}
    exists = summary.get("exists")
    if exists is not None:
        exists = bool(exists)
    count = summary.get("numTotalResults")
    if count is not None:
        try:
            count = int(count)
        except (TypeError, ValueError):
            count = None
    return exists, count
