import logging
from typing import Tuple

from more_itertools import all_equal

from genes.hgvs import HGVSVariant
from genes.hgvs.hgvs_converter import HGVSConverter, HgvsMatchRefAllele
from library.utils import all_equal
from snpdb.models import VariantCoordinate, GenomeBuild


class ComboCheckerHGVSConverter(HGVSConverter):
    def __init__(self, genome_build: GenomeBuild, converters, die_on_error=True):
        super().__init__(genome_build)
        self._converters = converters
        self.die_on_error = die_on_error

    def _call_converters(self, method_str, *args, **kwargs):
        results_and_methods = []
        for converter in self._converters:
            method = getattr(converter, method_str)
            results_and_methods.append((method(*args, **kwargs), converter.description()))

        results = [rm[0] for rm in results_and_methods]
        if not all_equal(results):
            msg = f"HGVS converters calling {method_str}({args}) returned different results! {results_and_methods})"
            if self.die_on_error:
                raise ValueError(msg)
            else:
                logging.error(msg)

        result = results[0]  # All same so any is fine
        # logging.debug("HGVS Combo: %s", result)
        return result

    def variant_coords_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        return self._call_converters("variant_coords_to_g_hgvs", vc)

    def variant_coords_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSVariant:
        return self._call_converters("variant_coords_to_c_hgvs", vc, transcript_version)

    def hgvs_to_variant_coords_and_reference_match(self, hgvs_string: str, transcript_version) -> Tuple[VariantCoordinate, HgvsMatchRefAllele]:
        return self._call_converters("hgvs_to_variant_coords_and_reference_match", hgvs_string, transcript_version)

    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        return self._call_converters("c_hgvs_remove_gene_symbol", hgvs_string)

    def get_transcript_accession(self, hgvs_string: str) -> str:
        return self._call_converters("get_transcript_accession", hgvs_string)

    def description(self) -> str:
        # This will be different so just return the 1st one
        return self._converters[0].description()
