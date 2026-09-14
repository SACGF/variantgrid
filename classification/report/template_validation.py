"""
The render a case template has to survive before it can be saved.

A report template is lab maintained config edited in the admin, and a template that only fails when
someone builds a real case fails in front of a patient's report. So both templates are rendered over
a fixture case here: the HTML has to convert through xhtml2pdf, and the JSON has to json.loads.

Deliberately free of any classification.models import - classification_report_models imports this,
so anything model-shaped here would be a cycle. FIXTURE_CONTEXT is therefore a hand written stand-in
for classification/report/case_report_context.py:context_as_dict, and a test asserts they match.
"""
import io
import json
from datetime import UTC, date, datetime
from typing import Optional

from django.template import TemplateSyntaxError, engines
from xhtml2pdf import pisa


def _fixture_evidence(**values) -> dict:
    """ Enough of an evidence blob for a template to read .value / .note / .formatted off it """
    return {key: {"value": value, "note": None, "formatted": str(value), "label": key.title()}
            for key, value in values.items()}


def _fixture_variant(gene_symbol: str, reported: bool = True) -> dict:
    return {
        "classification_modification_id": 1,
        "classification_id": 1,
        "kind": "small_variant",
        "alteration": "var",
        "gene_symbol": gene_symbol,
        "gene_symbols": [gene_symbol],
        "gene_label": gene_symbol,
        "tier": "tier_1",
        "amp_tier": "IA",
        "tier_rank": 10,
        "vaf": 0.42,
        "vaf_percent": 42.0,
        "copy_number": None,
        "fold_change": None,
        "reported": reported,
        "sample": {"pk": 1, "str": "fixture sample", "name": "fixture sample"},
        "evidence": _fixture_evidence(c_hgvs="c.123A>G", p_hgvs="p.Lys41Arg",
                                      refseq_transcript_id="NM_000000.1",
                                      gene_symbol=gene_symbol, exon="5 of 11",
                                      somatic_summary_interpretation="Fixture narrative.",
                                      h_summary="Fixture gene paragraph."),
        "warnings": [],
    }


def _fixture_gene_group(gene_symbol: str, variants: list[dict]) -> dict:
    return {
        "gene_symbol": gene_symbol,
        "variants": variants,
        "classifications": [v["evidence"] for v in variants],
        "gene_summary": "Fixture gene paragraph.",
        "gene_summary_source": 1,
        "warnings": [],
    }


def _build_fixture_context() -> dict:
    reported = _fixture_variant("GENE1")
    unreported = _fixture_variant("GENE2", reported=False)
    gene_groups = [_fixture_gene_group("GENE1", [reported]), _fixture_gene_group("GENE2", [unreported])]
    return {
        "source_level": "P",
        "case_report": {"pk": 1, "str": "fixture report", "status": "D",
                        "external_report_id": "FIXTURE-1", "report_date": date(2026, 1, 1)},
        "patient": {"pk": 1, "str": "fixture patient", "patient_code": "C00000",
                    "date_of_birth": date(1970, 1, 1), "sex": "U"},
        "specimen": {"pk": 1, "str": "fixture specimen", "reference_id": "BLOCK-1",
                     "collection_date": date(2026, 1, 1), "received_date": date(2026, 1, 2)},
        "extractions": [{"pk": 1, "str": "fixture extraction", "reference_id": "EXT-1",
                         "nucleic_acid_source": "D"}],
        "samples": [{"pk": 1, "str": "fixture sample", "name": "fixture sample"}],
        "sequencing_runs": ["FIXTURE_RUN"],
        "measures": {"tmb": {"pk": 1, "str": "TMB 7.0mut/Mb (Low)", "value": 7.0, "unit": "mut/Mb",
                             "call": "Low", "threshold": "10", "method": "fixture"},
                     "msi": {"pk": 2, "str": "MSI 1.0% (Stable)", "value": 1.0, "unit": "%",
                             "call": "Stable", "threshold": "20", "method": "fixture"}},
        "variants": [reported, unreported],
        # Every kind, so a template's "none detected" branch is rendered before it can be saved
        "kind_groups": [{"kind": "small_variant", "label": "Somatic Variants", "variants": [reported]},
                        {"kind": "copy_number", "label": "Copy Number Changes", "variants": []},
                        {"kind": "fusion", "label": "Gene Fusions", "variants": []}],
        "tier_groups": [{"tier": "tier_1", "label": "Tier I - Variants of Strong Clinical Significance",
                         "genes": [gene_groups[0]], "unreported_count": 0},
                        {"tier": "tier_2", "label": "Tier II - Variants of Potential Clinical Significance",
                         "genes": [], "unreported_count": 1},
                        {"tier": "tier_3", "label": "Tier III - Variants of Unknown Clinical Significance",
                         "genes": [], "unreported_count": 0}],
        "gene_groups": gene_groups,
        "summary": "Fixture summary interpretation.",
        "case_values": {"panel": "FIXTURE500", "mutations_comment": "Fixture comment"},
        "lab": {"pk": 1, "str": "fixture lab", "name": "fixture lab"},
        "user": {"pk": 1, "str": "fixture user", "username": "fixture user"},
        "generated": datetime(2026, 1, 1, tzinfo=UTC),
        "versions": {"variantgrid": "fixture", "annotation": {"GRCh38": "fixture"}, "callers": []},
    }


FIXTURE_CONTEXT = _build_fixture_context()


def render_html(template_str: str, context: dict) -> str:
    return engines['django'].from_string(template_str).render(context)


def html_to_pdf(html: str) -> bytes:
    """ Raises ValueError naming what xhtml2pdf choked on - it reports errors rather than raising """
    buffer = io.BytesIO()
    status = pisa.CreatePDF(html, dest=buffer, encoding='utf-8')
    if status.err:
        raise ValueError(f"xhtml2pdf could not convert this HTML ({status.err} error(s))")
    return buffer.getvalue()


def validate_case_template(template_str: str) -> Optional[str]:
    """ The message to show the editor, or None if the template is fine """
    try:
        html = render_html(template_str, FIXTURE_CONTEXT)
    except TemplateSyntaxError as e:
        return f"Template error: {e}"
    try:
        html_to_pdf(html)
    except Exception as e:  # pylint: disable=broad-except
        return f"This template does not convert to PDF: {e}"
    return None


def validate_json_template(template_str: str) -> Optional[str]:
    if not template_str.strip():
        return None  # Blank is the canonical context dump, not a broken template
    try:
        rendered = render_html(template_str, FIXTURE_CONTEXT)
    except TemplateSyntaxError as e:
        return f"Template error: {e}"
    try:
        json.loads(rendered)
    except json.JSONDecodeError as e:
        return f"This template does not render valid JSON: {e}"
    return None
