import logging

from genes.hgvs import HGVSVariant, HGVSException
from genes.hgvs.hgvs_converter import HGVSConverter, HgvsMatchRefAllele
from library.utils import all_equal
from snpdb.models import VariantCoordinate, GenomeBuild


class ComboCheckerHGVSConverter(HGVSConverter):
    def get_hgvs_converter_type(self) -> 'HGVSConverterType':
        from genes.hgvs.hgvs_converter import HGVSConverterType
        return HGVSConverterType.COMBO

    def get_version(self) -> str:
        return "1.0"

    def __init__(self, genome_build: GenomeBuild, converters, die_on_error=True):
        super().__init__(genome_build)
        self._converters = converters
        self.die_on_error = die_on_error

    def _call_converters(self, method_str, *args, **kwargs):
        results_and_methods = []
        for converter in self._converters:
            method = getattr(converter, method_str)
            try:
                result = method(*args, **kwargs)
            except HGVSException as e:
                result = e
            results_and_methods.append((result, converter.description()))

        results = [rm[0] for rm in results_and_methods]
        # If all exceptions, they won't be equal (as they are raised from different ones)
        has_exception = False
        has_non_exception = False
        for r in results:
            if isinstance(r, Exception):
                has_exception = True
            else:
                has_non_exception = True

        result = results[0]  # First result
        if has_exception and not has_non_exception:
            raise result

        if not all_equal(results):
            msg = f"HGVS converters calling {method_str}({args}) returned different results! {results_and_methods})"
            if self.die_on_error:
                raise ValueError(msg)
            else:
                logging.error(msg)

        if has_exception:
            raise result
        return result

    def create_hgvs_variant(self, hgvs_string: str) -> HGVSVariant:
        return self._call_converters("create_hgvs_variant", hgvs_string)

    def _variant_coordinate_to_g_hgvs(self, vc: VariantCoordinate) -> HGVSVariant:
        return self._call_converters("_variant_coordinate_to_g_hgvs", vc)

    def variant_coordinate_to_c_hgvs(self, vc: VariantCoordinate, transcript_version) -> HGVSVariant:
        return self._call_converters("variant_coordinate_to_c_hgvs", vc, transcript_version)

    def hgvs_to_variant_coordinate_and_reference_match(self, hgvs_string: str, transcript_version) -> tuple[VariantCoordinate, HgvsMatchRefAllele]:
        return self._call_converters("hgvs_to_variant_coordinate_and_reference_match", hgvs_string, transcript_version)

    def c_hgvs_remove_gene_symbol(self, hgvs_string: str) -> str:
        return self._call_converters("c_hgvs_remove_gene_symbol", hgvs_string)

    def get_transcript_accession(self, hgvs_string: str) -> str:
        return self._call_converters("get_transcript_accession", hgvs_string)

    def description(self) -> str:
        # This will be different so just return the 1st one
        return self._converters[0].description()
