"""
Fake ClassificationModifications for the case report tests.

Ordering, amp_tier and kind detection are pure functions of a record's evidence and its Variant, so
these tests build the record rather than a database full of classifications, samples and VCFs.
"""
from datetime import datetime
from types import SimpleNamespace
from typing import Optional

from django.utils import timezone

from classification.enums import SpecialEKeys
from classification.report.case_report_context import ReportVariant


class FakeModification:
    """ Enough ClassificationModification for the context builder: get(), a pk and a classification
        carrying the resolved Variant """

    def __init__(self, pk: int = 1, values: Optional[dict] = None, variant=None, sample=None,
                 modified: Optional[datetime] = None, lab_record_id: str = "fake"):
        self.pk = pk
        self.classification_id = pk
        self._values = values or {}
        self.modified = modified or timezone.now()
        self.classification = SimpleNamespace(pk=pk, variant=variant, sample=sample,
                                              cr_lab_id=lab_record_id)

    def get(self, key, default=None):
        return self._values.get(key, default)


def fake_gene_level_variant(alt: str, gene_symbols: Optional[list[str]] = None):
    """ A gene-level Variant - the alt is what says whether it is a fusion or a copy number call """
    gene_fusion = None
    if gene_symbols:
        gene_fusion = SimpleNamespace(gene_level_ids=[SimpleNamespace(symbol_str=symbol)
                                                      for symbol in gene_symbols])
    return SimpleNamespace(is_gene_level=True, alt=SimpleNamespace(seq=alt), genefusion=gene_fusion)


def fake_report_variant(gene_symbol: str, tier: Optional[str] = None,
                        amp_levels: Optional[list[str]] = None, vaf: Optional[float] = None,
                        copy_number: Optional[int] = None, c_hgvs: Optional[str] = None,
                        reported: bool = True, variant=None, pk: int = 1,
                        modified: Optional[datetime] = None,
                        gene_summary: Optional[str] = None) -> ReportVariant:
    values = {SpecialEKeys.GENE_SYMBOL: gene_symbol}
    if tier:
        values[SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE] = tier
    for level in amp_levels or []:
        values[f"amp:level_{level.lower()}"] = "met"
    if vaf is not None:
        values[SpecialEKeys.ALLELE_FREQUENCY] = vaf
    if copy_number is not None:
        values[SpecialEKeys.COPY_NUMBER] = copy_number

    evidence = {"gene_symbol": {"value": gene_symbol, "note": None, "formatted": gene_symbol,
                                "label": "Gene"},
                "c_hgvs": {"value": c_hgvs, "note": None, "formatted": c_hgvs, "label": "c.HGVS"},
                "p_hgvs": {"value": None, "note": None, "formatted": "", "label": "p.HGVS"},
                "h_summary": {"value": gene_summary, "note": None, "formatted": gene_summary or "",
                              "label": "Gene summary"},
                "somatic_summary_interpretation": {"value": f"{gene_symbol} narrative.", "note": None,
                                                   "formatted": "", "label": "Interpretation"}}
    record = FakeModification(pk=pk, values=values, variant=variant, modified=modified)
    return ReportVariant.build(record, user=None, reported=reported, evidence=evidence)
