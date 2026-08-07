""" Outbound Beacon query targets: a per-server gate + query-builder (§9).

Before we fan a variant out to an external Beacon we test whether that server is even relevant
to it - a germline SNV against a cancer copy-number Beacon is wasted latency and a guaranteed
miss. Each target pairs a gate (`accepts`) with a coordinate query (`build_params`); a node in
settings.BEACON_QUERY_NODES selects its target via `type` (and may restrict `assemblies`).
`eligible_queries()` loops the configured nodes and returns params only for the servers whose
gate the variant passes, so a variant page shows only the Beacons that apply to that variant.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

from beacon.variant_mapping import variant_to_beacon_query_params
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models import Variant, GenomeBuild


class BeaconQueryTarget(ABC):
    """ A gate + query-builder for one class of external Beacon. """

    # Short noun phrase for the variants this Beacon class handles - shown on the variant page
    # to explain why a node was skipped (e.g. "only queried for a copy-number DEL/DUP variant").
    domain: str = ""

    @abstractmethod
    def accepts(self, variant: Variant, genome_build: GenomeBuild) -> bool:
        """ The gate: is this variant in this server's domain? """

    @abstractmethod
    def build_params(self, variant: Variant, genome_build: GenomeBuild) -> dict:
        """ The g_variants query params to send (only a public coordinate leaves). """


class SnvBeaconTarget(BeaconQueryTarget):
    """ Sequence Beacons: an exact referenceBases/alternateBases coordinate query. """

    domain = "sequence variant (SNV/indel)"

    def accepts(self, variant, genome_build) -> bool:
        vc = variant.coordinate
        return not vc.is_symbolic and bool(vc.ref and vc.alt)

    def build_params(self, variant, genome_build) -> dict:
        return variant_to_beacon_query_params(variant, genome_build)


class CnvBeaconTarget(BeaconQueryTarget):
    """ Copy-number Beacons (bycon/progenetix et al.): a symbolic <DEL>/<DUP> becomes a range
        query for any copy-number event of that class overlapping the variant's span. The
        per-sample copy number (on CohortGenotype) stays internal - only the public span and
        copy-number class leave. """

    domain = "copy-number DEL/DUP variant"

    # VG symbolic alt -> GA4GH/EFO copy-number class (validated live against progenetix).
    VARIANT_TYPE_EFO = {
        VCFSymbolicAllele.DEL: "EFO:0030067",  # copy number loss
        VCFSymbolicAllele.DUP: "EFO:0030070",  # copy number gain
    }

    def _efo(self, vc) -> Optional[str]:
        return self.VARIANT_TYPE_EFO.get(vc.alt)

    def accepts(self, variant, genome_build) -> bool:
        vc = variant.coordinate
        return bool(vc.is_symbolic and vc.svlen is not None and self._efo(vc))

    def build_params(self, variant, genome_build) -> dict:
        vc = variant.coordinate
        return {
            "referenceName": vc.chrom.replace("chr", ""),
            "assemblyId": genome_build.name,
            "start": vc.position - 1,  # 1-based -> 0-based
            "end": vc.end,             # position + abs(svlen)
            "variantType": self._efo(vc),
        }


TARGETS: dict[str, BeaconQueryTarget] = {
    "snv": SnvBeaconTarget(),
    "cnv": CnvBeaconTarget(),
}


@dataclass
class NodeEligibility:
    """ Whether one configured node is queried for a variant. When `eligible`, `params` carries
        the query to send; otherwise `reason` is a human-readable explanation for the page. """
    node_id: str
    node_type: Optional[str]
    eligible: bool
    params: Optional[dict] = None
    reason: Optional[str] = None


def evaluate_queries(variant: Variant, genome_build: GenomeBuild, node_configs: dict) -> list[NodeEligibility]:
    """ Per-node eligibility for this variant. A node is skipped (with a `reason`) when its
        `type` is unknown/omitted, its `assemblies` (if set) do not include this build, or its
        target's gate rejects the variant. The variant page uses the reasons to show which
        Beacons apply and why the rest were not queried; eligible_queries() is the params-only
        projection used by the headless outbound path. """
    results = []
    for node_id, node in node_configs.items():
        node_type = node.get("type")
        target = TARGETS.get(node_type)
        if target is None:
            results.append(NodeEligibility(node_id, node_type, False,
                                           reason=f"Beacon type {node_type!r} is not supported"))
            continue
        assemblies = node.get("assemblies")
        if assemblies and genome_build.name not in assemblies:
            results.append(NodeEligibility(node_id, node_type, False,
                                           reason=f"only queried for {', '.join(assemblies)} "
                                                  f"(this variant is {genome_build.name})"))
            continue
        if not target.accepts(variant, genome_build):
            results.append(NodeEligibility(node_id, node_type, False,
                                           reason=f"only queried for a {target.domain}"))
            continue
        results.append(NodeEligibility(node_id, node_type, True,
                                       params=target.build_params(variant, genome_build)))
    return results


def eligible_queries(variant: Variant, genome_build: GenomeBuild, node_configs: dict) -> dict:
    """ {node_id: query_params} for every configured node whose target accepts this variant. """
    return {e.node_id: e.params
            for e in evaluate_queries(variant, genome_build, node_configs) if e.eligible}
