""" Outbound Beacon query targets: a per-server gate + query-builder (§9).

Before we fan a variant out to an external Beacon we test whether that server is even relevant
to it - a germline SNV against a cancer copy-number Beacon is wasted latency and a guaranteed
miss. Each target pairs a gate (`accepts`) with a coordinate query (`build_params`); a node in
settings.BEACON_QUERY_NODES selects its target via `type` (and may restrict `assemblies`).
`eligible_queries()` loops the configured nodes and returns params only for the servers whose
gate the variant passes, so a variant page shows only the Beacons that apply to that variant.
"""
import re
from abc import ABC, abstractmethod
from typing import Optional

from beacon.variant_mapping import variant_to_beacon_query_params
from snpdb.models import Variant, GenomeBuild


class BeaconQueryTarget(ABC):
    """ A gate + query-builder for one class of external Beacon. """

    @abstractmethod
    def accepts(self, variant: Variant, genome_build: GenomeBuild) -> bool:
        """ The gate: is this variant in this server's domain? """

    @abstractmethod
    def build_params(self, variant: Variant, genome_build: GenomeBuild) -> dict:
        """ The g_variants query params to send (only a public coordinate leaves). """


class SnvBeaconTarget(BeaconQueryTarget):
    """ Sequence Beacons: an exact referenceBases/alternateBases coordinate query. """

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

    # VG symbolic alt -> GA4GH/EFO copy-number class (validated live against progenetix).
    VARIANT_TYPE_EFO = {
        "DEL": "EFO:0030067",  # copy number loss
        "DUP": "EFO:0030070",  # copy number gain
    }
    _SYMBOLIC_ALT = re.compile(r"<(DEL|DUP)>", re.IGNORECASE)

    def _efo(self, vc) -> Optional[str]:
        match = self._SYMBOLIC_ALT.search(vc.alt or "")
        return self.VARIANT_TYPE_EFO.get(match.group(1).upper()) if match else None

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


def eligible_queries(variant: Variant, genome_build: GenomeBuild, node_configs: dict) -> dict:
    """ Return {node_id: query_params} for every configured node whose target accepts this
        variant. A node is skipped when its `type` is unknown/omitted, its `assemblies` (if set)
        do not include this build, or its target's gate rejects the variant. """
    eligible = {}
    for node_id, node in node_configs.items():
        target = TARGETS.get(node.get("type"))
        if target is None:
            continue
        assemblies = node.get("assemblies")
        if assemblies and genome_build.name not in assemblies:
            continue
        if target.accepts(variant, genome_build):
            eligible[node_id] = target.build_params(variant, genome_build)
    return eligible
