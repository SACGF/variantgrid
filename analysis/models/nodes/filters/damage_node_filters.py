""" The Effect node's filters as people read them (#2836): one EffectFilter per row of the editor - the node fields it
    sets, the group it sits in and how it reads in words. The node's query (DamageNode._get_node_q) doesn't use this:
    it is every "All of" filter (the _required ones, plus the filters that always narrow) and at least one "Any of".

    Entry points: get_effect_filters(node) (the rows this node's annotation version offers, in editor order),
    get_effect_filter_terms(node) (each set filter in words) and describe_effect_filters(node) (the sentence the
    editor shows and the node's summary and hover card use). The wording lives here only: the editor gets each
    filter's `words` and SENTENCE and fills in the live values (damagenode_editor.js), so the two can't drift. """
import html
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from functools import cache
from typing import Any, Optional

from annotation.enums import EVEClass, Pathogenicity, PathogenicityImpact
from annotation.pathogenicity_predictions import TOOLS
from snpdb.variant_filters import get_variant_type_labels, restricts_variant_types

ANY_OF = "any"
ALL_OF = "all"

# The sentence's frame, around the filters' own words
SENTENCE = {
    "none": "No filters. Every variant from the node above passes.",
    "lead": "Keeps variants with ",
    "any_of": "at least one of",
    "all_of": "all of",
    "and": "and",
    "separator": " · ",
}

# A prediction dropdown starts at the value the damage count calls pathogenic - EVE's class isn't counted
_PRED_DEFAULTS = {"eve_class": EVEClass.PATHOGENIC}

# Score filters outside the coding sequence, which the editor's Add filter palette lists apart from the missense ones
_NON_CODING_RAW_FIELDS = frozenset({"promoter_ai_score", "avi_phred"})


@dataclass(frozen=True)
class EffectFilter:
    label: str  # the row's label, the tool's own name
    category: str  # heading in the Add filter palette
    value_fields: tuple[str, ...]  # DamageNode fields holding the value - the row is drawn while any is set
    words: tuple[str, ...]  # each value field in words, "{}" standing for its value
    op: str = ""  # comparison after the label ("≥")
    required_field: Optional[str] = None  # set: in All of, clear: in Any of. None: always in All of
    allow_null_field: Optional[str] = None  # only has an effect in All of
    null_label: str = "keep if not scored"
    null_words: str = "or not scored"
    help: str = ""
    default: Any = None  # a dropdown's value when the filter is added (None: its first option)

    @property
    def key(self) -> str:
        return self.value_fields[0]

    @property
    def fixed(self) -> bool:
        """ Always narrows, so it can only be in All of (variant type, repeat and SV artefact regions) """
        return self.required_field is None

    def group(self, node) -> str:
        if self.fixed or getattr(node, self.required_field):
            return ALL_OF
        return ANY_OF

    def term(self, node) -> Optional[str]:
        """ This filter in words, or None when the node has no value for it """
        parts = []
        for field, words in zip(self.value_fields, self.words):
            if (value := _value_words(node, field)) is not None:
                parts.append(words.replace("{}", value))
        if not parts:
            return None
        text = " or ".join(parts)
        if self.allow_null_field and self.group(node) == ALL_OF and getattr(node, self.allow_null_field):
            text += f" ({self.null_words})"
        return text


def _value_words(node, field: str) -> Optional[str]:
    value = getattr(node, field)
    if field == "variant_types":
        if restricts_variant_types(value):
            return " or ".join(get_variant_type_labels(value))
        return None
    if value is None or value is False or value == "":
        return None
    if value is True:
        return ""
    if node._meta.get_field(field).choices:
        return getattr(node, f"get_{field}_display")()
    return f"{value:g}"


def _slider(label: str, category: str, field: str, stem: str, op: str = "≥", allow_null: bool = True,
            **kwargs) -> EffectFilter:
    """ A threshold on one score, with its own <stem>_required and <stem>_allow_null """
    return EffectFilter(label, category, (field,), (f"{label} {op} {{}}",), op=op, required_field=f"{stem}_required",
                        allow_null_field=f"{stem}_allow_null" if allow_null else None, **kwargs)


def _flag(label: str, category: str, field: str, words: str, fixed: bool = False, **kwargs) -> EffectFilter:
    return EffectFilter(label, category, (field,), (words,), required_field=None if fixed else f"{field}_required",
                        **kwargs)


def _always() -> bool:
    return True


def _all_effect_filters(node, editor: bool) -> Iterator[tuple[EffectFilter, Callable[[], bool]]]:
    """ The filters in editor order, each with whether the editor offers it. For the editor, those the node's columns
        version has. Otherwise (describing the node) every filter, reading nothing off the annotation version that
        the words don't need - the node card asks for its chips on every analysis load """
    columns_version = node.columns_version if editor else None
    always = _always

    def from_version(version: int) -> bool:
        return not editor or columns_version >= version

    yield EffectFilter("Impact", "Consequence", ("impact_min",), ("impact ≥ {}",), op="≥",
                       required_field="impact_required", default=PathogenicityImpact.HIGH,
                       help="VEP's impact. MODERATE* does not include missense variants, so use MODERATE* together "
                            "with one or more pathogenicity predictors to include damaging missense variants"), always

    lof = "Loss of function"
    if from_version(2):
        yield EffectFilter("PTC NMD", lof, ("ptc_nmd",), ("PTC NMD: {}",), required_field="ptc_nmd_required",
                           help="Nonsense Mediated Decay predicted from where a frameshift's premature termination "
                                "codon lands. 'Escapes NMD' keeps frameshifts making a truncated protein, 'Escapes NMD "
                                "or NMD predicted' every frameshift whose new stop codon was located"), \
            lambda: columns_version >= 5 or bool(node.ptc_nmd)
        yield _flag("NMD escaping variant (LOFTEE)", lof, "nmd_escaping_variant", "NMD escaping variant (LOFTEE)",
                    help="Nonsense Mediated Decay escaping variant"), node.shows_nmd_escaping_variant
        yield EffectFilter("ALoFT", lof, ("aloft",), ("ALoFT: {}",), required_field="aloft_required",
                           allow_null_field="aloft_allow_null",
                           help="Annotation of Loss-of-Function Transcripts - predicts whether a loss-of-function "
                                "variant is tolerated, or causes recessive or dominant disease. Only high-confidence "
                                "predictions match"), node.shows_aloft

    splicing = "Splicing"
    splice_fields = ("splice_min", "alphagenome_splicing_min")
    splice_words = ("splicing ≥ {}", "AlphaGenome splicing ≥ {}")
    splice_filter = EffectFilter("Splicing (SpliceAI, dbscSNV)", splicing, splice_fields, splice_words, op="≥",
                                 required_field="splice_required", allow_null_field="splice_allow_null",
                                 help="SpliceAI or dbscSNV score, or an indel in a splice region")
    yield splice_filter, always
    if from_version(2):
        yield EffectFilter("SpliceVault", splicing, ("splicevault_out_of_frame",), ("SpliceVault: {}",),
                           required_field="splicevault_out_of_frame_required",
                           allow_null_field="splicevault_out_of_frame_allow_null",
                           help="SpliceVault's top mis-splicing events at a splice site SpliceAI predicts lost. 'All "
                                "top events' keeps a variant whose top events are all out of frame (Dawes 2023's "
                                "consistent-with-null case), 'Top event' one whose most frequent event is"), \
            node.has_splicevault

    missense = "Missense predictors"
    num_prediction_fields = None
    if editor or node.damage_predictions_min is not None:
        num_prediction_fields = node.num_prediction_fields
    of_total = f" of {num_prediction_fields}" if num_prediction_fields is not None else ""
    yield EffectFilter("Damage predictions", missense, ("damage_predictions_min",),
                       (f"at least {{}}{of_total} damage predictions",), op="≥",
                       required_field="damage_predictions_required", allow_null_field="damage_predictions_allow_null",
                       null_label="don't count unscored tools as benign",
                       null_words="unscored tools not counted as benign",
                       help=node.damage_predictions_description() if editor else ""), node.shows_damage_predictions
    if not editor or columns_version == 1:
        yield _slider("CADD (phred)", missense, "cadd_score_min", stem="cadd_score"), always
    if editor and columns_version == 1:  # otherwise the REVEL raw score below, the same field
        yield _slider("REVEL", missense, "revel_score_min", stem="revel_score"), always

    raw_score_filters = []
    rankscore_filters = []
    if from_version(4):
        editor_raw_tools = cache(node.get_editor_raw_score_tools)
        for tool in TOOLS:
            if tool.raw_field:
                category = "Non-coding" if tool.raw_field in _NON_CODING_RAW_FIELDS else missense
                raw_score_filters.append((
                    _slider(tool.node_label or tool.name, category, tool.node_threshold_field, tool.raw_field,
                            op=tool.raw_direction.comparison, help=tool.note),
                    lambda tool=tool: tool in editor_raw_tools()))
        yield from (f for f in raw_score_filters if f[0].category == missense)

        editor_pred_tools = cache(node.get_editor_pred_tools)
        for tool in TOOLS:
            if tool.pred_field:
                label = f"{tool.name} class" if tool.pred_field.endswith("class") else f"{tool.name} prediction"
                yield EffectFilter(label, missense, (tool.pred_field,), (f"{label}: {{}}",),
                                   required_field=f"{tool.pred_field}_required",
                                   allow_null_field=f"{tool.pred_field}_allow_null",
                                   default=_PRED_DEFAULTS.get(tool.pred_field,
                                                              next(iter(tool.pred_pathogenic_values), None))), \
                    lambda tool=tool: tool in editor_pred_tools()

    if from_version(2):
        # Columns v2/v3 have only rankscores. From v4 they are legacy, offered while a node still has one set
        editor_rankscore_tools = cache(node.get_editor_rankscore_tools)
        category = "Legacy rankscores" if from_version(4) else missense
        rankscore_filters = [(_slider(f"{tool.name} rankscore", category, f"{tool.rankscore_field}_min",
                                      tool.rankscore_field),
                              lambda tool=tool: tool in editor_rankscore_tools())
                             for tool in TOOLS if tool.rankscore_field]
        if not from_version(4):
            yield from rankscore_filters

    non_coding = "Non-coding"
    if from_version(2):
        yield EffectFilter("5' UTR (UTRAnnotator)", non_coding, ("utr_high_impact",), ("{}",),
                           required_field="utr_high_impact_required", allow_null_field="utr_high_impact_allow_null",
                           help="UTRAnnotator: uAUG gain or uFrameshift creating an oORF, or uSTOP loss with no "
                                "alternative stop, each with strong/moderate Kozak; or uAUG loss with strong Kozak "
                                "(Martin-Geary 2025). Out-of-frame only drops an oORF in frame with the CDS - an "
                                "N-terminal extension, not necessarily loss of function (Ellingford 2022)"), \
            node.has_utr_annotator
    yield from (f for f in raw_score_filters if f[0].category == non_coding)

    gene = "Gene"
    gene_null = {"null_label": "keep genes with no score", "null_words": "or gene has no score"}
    yield _slider("gnomAD LOEUF", gene, "gene_loeuf_max", op="≤", stem="gene_loeuf",
                  help="gnomAD v4 LOEUF of the gene's MANE Select transcript - the upper bound of the loss-of-function "
                       "observed/expected 90% CI. Lower is more intolerant of loss of function (PVS1's 'LoF is a "
                       "known mechanism'). Starts at gnomAD v4's constrained-gene cutoff, 0.6", **gene_null), \
        lambda: node.has_gene_constraint or node.gene_loeuf_max is not None
    yield _slider("gnomAD missense Z", gene, "gene_mis_z_min", stem="gene_mis_z",
                  help="gnomAD v4 missense Z of the gene's MANE Select transcript - higher means fewer missense "
                       "variants than expected (PP2 / BP1). Starts at 3.09, ExAC's missense-constrained gene cutoff",
                  **gene_null), \
        lambda: node.has_gene_constraint or node.gene_mis_z_min is not None
    yield _flag("Gene has an OMIM phenotype", gene, "omim_phenotype", "gene has an OMIM phenotype"), always

    known = "Known variant"
    yield _flag("In COSMIC", known, "in_cosmic", "in COSMIC",
                help="Has a COSMIC id (from VEP's Existing_variation) - unlike COSMIC count, which needs the "
                     "licensed COSMIC data"), always
    yield _slider("COSMIC count", known, "cosmic_count_min", stem="cosmic_count", allow_null=False), \
        node.shows_cosmic_count
    yield _flag("Published (PubMed, Mastermind)", known, "published", "published (PubMed, Mastermind)"), always
    yield _flag("In a protein domain", known, "protein_domain", "in a protein domain"), always

    region = "Type and region"
    yield EffectFilter("Variant type", region, ("variant_types",), ("variant type is {}",)), always
    yield _flag("Outside repeat regions", region, "outside_repeat_region", "outside repeat regions", fixed=True,
                help="Keep only variants with no RepeatMasker overlap - PM4's in-frame indel or stop loss 'in a "
                     "non-repeat region'"), \
        lambda: node.has_repeat_masker() or node.outside_repeat_region

    sv = "Structural variants (AnnotSV)"
    sv_null = {"null_label": "keep small variants", "null_words": "or not a structural variant"}
    yield EffectFilter("AnnotSV ACMG class", sv, ("annotsv_acmg_class_min",), ("AnnotSV ACMG class ≥ {}",), op="≥",
                       required_field="annotsv_acmg_class_required", allow_null_field="annotsv_acmg_class_allow_null",
                       default=Pathogenicity.LIKELY_PATHOGENIC,
                       help="AnnotSV's ACMG class for the SV, from its ranking score on the ACMG/ClinGen CNV points "
                            "system. Every small variant has none", **sv_null), \
        lambda: node.has_annotsv or node.annotsv_acmg_class_min is not None
    yield _slider("AnnotSV ranking score", sv, "annotsv_acmg_score_min", stem="annotsv_acmg_score",
                  help="AnnotSV ranking score: 0.90 and above is likely pathogenic, 0.99 and above pathogenic. Every "
                       "small variant has none", **sv_null), \
        lambda: node.has_annotsv or node.annotsv_acmg_score_min is not None
    yield _flag("No SV breakpoint in an artefact region", sv, "exclude_sv_artefact_regions",
                "no SV breakpoint in an artefact region", fixed=True,
                help="No SV breakpoint in an ENCODE blacklist region, segmental duplication or repeat (AnnotSV) - "
                     "where SV callers most often report false positives. Small variants pass"), \
        lambda: node.has_annotsv or node.exclude_sv_artefact_regions

    if from_version(4):
        yield from rankscore_filters


def get_effect_filters(node) -> list[EffectFilter]:
    """ The editor's rows: the filters this annotation version populates, plus any the node has a value for
        (DamageNode.shows_filter) so it can be seen and removed """
    return [effect_filter for effect_filter, offered in _all_effect_filters(node, editor=True) if offered()]


def get_effect_filter_terms(node) -> dict[str, tuple[str, str]]:
    """ key -> (group, words) for each filter the node has a value for, in editor order """
    terms = {}
    for effect_filter, _ in _all_effect_filters(node, editor=False):
        if (term := effect_filter.term(node)) is not None:
            terms[effect_filter.key] = (effect_filter.group(node), term)
    return terms


def effect_filter_sentence(any_of: list[str], all_of: list[str], as_html: bool = False) -> str:
    """ "Keeps variants with at least one of: <any of>, and all of: <all of>." - damagenode_editor.js builds the
        same from SENTENCE """
    if not (any_of or all_of):
        return SENTENCE["none"]
    escape = html.escape if as_html else (lambda s: s)
    bold = (lambda s: f"<b>{s}</b>") if as_html else (lambda s: s)

    def group(terms: list[str], name: str) -> str:
        listed = SENTENCE["separator"].join(escape(t) for t in terms)
        return listed if len(terms) == 1 else f"{bold(SENTENCE[name])}: {listed}"

    text = SENTENCE["lead"]
    if any_of:
        text += group(any_of, "any_of")
    if all_of:
        if any_of:
            text += f", {bold(SENTENCE['and'])} "
        text += group(all_of, "all_of")
    return text + "."


def describe_effect_filters(node, as_html: bool = False) -> str:
    terms = get_effect_filter_terms(node).values()
    any_of = [term for group, term in terms if group == ANY_OF]
    all_of = [term for group, term in terms if group == ALL_OF]
    return effect_filter_sentence(any_of, all_of, as_html=as_html)
