"""
What a report template is handed, and the order it is handed it in.

ReportVariant is one classification's row: the evidence blob every template reads, plus the facts
the document is organised by - what kind of event it is, its AMP sub-tier, its gene(s) and its VAF.
Ordering, grouping and tier derivation live here rather than in the templates, which only loop, so
the HTML preview, the PDF, the Word file and the JSON cannot disagree about what the report says.

Entry points:
- evidence_row_data(modification, user) - the per-record evidence dict (the single record report's
  `record`, and every ReportVariant.evidence)
- build_report_variants(modifications, user, reported_by_pk) - ReportVariants in report order
- build_report_context(...) - the whole ReportContext for a case
- context_as_dict(report_context) - what the templates and CaseReport.context_snapshot see

Building, rebuilding and finalising the CaseReport itself is classification/report/case_report_builder.py.
"""
import re
from dataclasses import dataclass, field
from datetime import datetime
from decimal import ROUND_HALF_UP, Decimal
from typing import Any, Optional

from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist
from django.db.models import F, Q
from django.utils import timezone

from annotation.models import VariantAnnotationVersion
from classification.enums import SpecialEKeys
from classification.enums.classification_enums import SomaticClinicalSignificance
from classification.models.classification import Classification, ClassificationModification
from classification.models.classification_json import ClassificationJsonParams
from classification.models.classification_report_models import CaseReportStatus, Measure
from classification.models.evidence_key import EvidenceKeyMap
from library.django_utils import get_cached_project_git_hash
from library.genomics.vcf_enums import GeneLevelSymbolicAlt, VariantClass
from patients.models import Extraction, Patient, Specimen, SpecimenMeasure
from patients.models_enums import (
    MEASURE_CONTEXT_KEYS,
    SampleSourceLevel,
)
from patients.sample_grouping import get_patient_for_source, get_sample_group
from seqauto.models import DragenTSO500CombinedVariantOutput, LibraryQC
from seqauto.models.models_enums import CVO_MEASURE_CONTEXT_KEYS, LIBRARY_QC_CONTEXT_KEYS
from snpdb.models import GenomeBuild, Lab, Sample


class ReportVariantKind:
    """ What sort of event the report is printing. The Results Summary is one table per kind """
    SMALL_VARIANT = 'small_variant'
    COPY_NUMBER = 'copy_number'
    COPY_NUMBER_LOSS = 'copy_number_loss'
    FUSION = 'fusion'
    SPLICE = 'splice'

    ORDER = [SMALL_VARIANT, COPY_NUMBER, COPY_NUMBER_LOSS, FUSION, SPLICE]
    LABELS = {
        SMALL_VARIANT: "Somatic Variants",
        COPY_NUMBER: "Copy Number Changes",
        COPY_NUMBER_LOSS: "Copy Number Losses",
        FUSION: "Gene Fusions",
        SPLICE: "Splicing Variants",
    }


class Alteration:
    """ The short form the JSON record uses for what changed - the four values the TSO 500 reports
        have ever carried, so a downstream consumer never meets one it has not seen. A gene level
        loss goes out as an amplification with its real copy number below two, rather than as a
        fifth value no consumer has a parser for """
    VARIANT = 'var'
    AMPLIFICATION = 'amp'
    FUSION = 'fusion'
    SPLICE = 'spl'


# amp_tier, ranked. The bare tiers are the fallbacks §Tier gives a classification whose AMP levels
# do not say which sub-tier it is - they sort with their family so the report order still reads right
AMP_TIER_RANK = {
    'IA': 10, 'IB': 11, 'I': 12,
    'I/II': 13,
    'IIC': 20, 'IID': 21, 'II': 22,
    'III': 30,
    'IV': 40,
}
UNTIERED_RANK = 99
# Behind every real copy number, counting either way about - @see ReportVariant.copies_rank
NO_COPY_NUMBER_RANK = 10 ** 9

TIER_GROUP_LABELS = {
    SomaticClinicalSignificance.TIER_1: "Tier I - Variants of Strong Clinical Significance",
    SomaticClinicalSignificance.TIER_1_OR_2: "Tier I/II - Variants of Strong or Potential Clinical Significance",
    SomaticClinicalSignificance.TIER_2: "Tier II - Variants of Potential Clinical Significance",
    SomaticClinicalSignificance.TIER_3: "Tier III - Variants of Unknown Clinical Significance",
    SomaticClinicalSignificance.TIER_4: "Tier IV - Benign or Likely Benign Variants",
}
# Printed whether or not the case has any - a report that skips Tier II does not say "none found"
ALWAYS_PRINTED_TIERS = [SomaticClinicalSignificance.TIER_1,
                        SomaticClinicalSignificance.TIER_2,
                        SomaticClinicalSignificance.TIER_3]
TIER_GROUP_ORDER = [SomaticClinicalSignificance.TIER_1,
                    SomaticClinicalSignificance.TIER_1_OR_2,
                    SomaticClinicalSignificance.TIER_2,
                    SomaticClinicalSignificance.TIER_3,
                    SomaticClinicalSignificance.TIER_4]
UNTIERED = ""
UNTIERED_LABEL = "Not tiered"

# The value of variant_reported that means "seen, not on the report"
NOT_REPORTED = "not_included"

GENE_SUMMARY_KEY = SpecialEKeys.H_SUMMARY

# The protein designation a curator gives a variant that changes no codon. A small variant carrying
# it is an intronic call the report has to say something about, or the reader is left with a blank
# protein column. TERT is left out: its promoter variants are non-coding for a different reason
SPLICE_NOTE_PROTEIN = "p.(?)"
SPLICE_NOTE_EXCLUDED_GENES = {"TERT"}


def evidence_row_data(record: ClassificationModification, user: User) -> dict:
    """ Every evidence key as {value, note, formatted, label}, plus the derived blocks a report
        template reads. Keys have ':' replaced with '_' because Vue/JS can't handle it in a name """
    context = {}
    evidence = record.as_json(ClassificationJsonParams(user, include_data=True))['data']
    e_keys = EvidenceKeyMap.instance().with_overrides(record.classification.evidence_key_overrides)

    for e_key in e_keys.all_keys:
        blob = evidence.get(e_key.key) or {}
        context[e_key.key.replace(':', '_')] = {
            'value': blob.get('value', None),
            'note': blob.get('note', None),
            'formatted': e_key.pretty_value(blob),
            'label': e_key.pretty_label,
        }

    for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
        c_hgvs = record.classification.get_resolved_hgvs(genome_build)
        context["c_hgvs_" + genome_build.pk.lower()] = {
            'value': c_hgvs,
            'note': None,
            'formatted': c_hgvs,
            'label': "c.HGVS",
        }

    context['condition_resolved'] = record.classification.condition_resolution
    context['citations'] = record.loaded_citations().to_json()
    context['evidence_weights'] = Classification.summarize_evidence_weights(evidence)
    context['acmg_criteria'] = record.criteria_strength_summary(e_keys)
    context['editable'] = record.classification.can_write(user)
    return context


def amp_tier(record: ClassificationModification) -> tuple[str, list[str]]:
    """ The sub-tier the printed report needs ('IIC'), which is the tier and the AMP evidence level
        together. Returns (amp_tier, warnings) - a tier whose levels don't say which sub-tier it is
        gets the bare tier and a warning, so nothing is ever printed as a silently wrong sub-tier """
    tier = record.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
    if not tier:
        return "", []

    levels = {level for key, level in SpecialEKeys.AMP_LEVELS_TO_LEVEL.items() if record.get(key)}
    short = SomaticClinicalSignificance.SHORT_LABELS.get(tier)
    if short is None:
        return "", [f"Unknown somatic clinical significance '{tier}'"]

    if tier == SomaticClinicalSignificance.TIER_1:
        for level in ("A", "B"):
            if level in levels:
                return f"I{level}", []
        return short, ["Tier I with no AMP level A or B - the report cannot say whether it is IA or IB"]
    if tier == SomaticClinicalSignificance.TIER_2:
        for level in ("C", "D"):
            if level in levels:
                return f"II{level}", []
        return short, ["Tier II with no AMP level C or D - the report cannot say whether it is IIC or IID"]
    if tier == SomaticClinicalSignificance.TIER_1_OR_2:
        return short, ["Tier I/II has not been resolved to a tier, so it has no AMP sub-tier"]
    return short, []


def _kind_and_alteration(record: ClassificationModification) -> tuple[str, str]:
    """ A gene-level alt says outright what the event is - a fusion, a copy number call or a splice
        junction; otherwise the variant_class evidence key is what a classification with no resolved
        variant has to go on """

    variant = record.classification.variant
    if variant is not None and variant.is_gene_level:
        if parsed := GeneLevelSymbolicAlt.parse(variant.alt.seq):
            kind_alt = parsed[0]
            if kind_alt in (GeneLevelSymbolicAlt.FUSION, GeneLevelSymbolicAlt.FUSION_UNORDERED):
                return ReportVariantKind.FUSION, Alteration.FUSION
            if kind_alt == GeneLevelSymbolicAlt.GAIN:
                return ReportVariantKind.COPY_NUMBER, Alteration.AMPLIFICATION
            if kind_alt == GeneLevelSymbolicAlt.LOSS:
                return ReportVariantKind.COPY_NUMBER_LOSS, Alteration.AMPLIFICATION
            if kind_alt == GeneLevelSymbolicAlt.SPLICE:
                return ReportVariantKind.SPLICE, Alteration.SPLICE

    variant_class = record.get(SpecialEKeys.VARIANT_CLASS)
    if variant_class == VariantClass.COPY_NUMBER_GAIN.label:
        return ReportVariantKind.COPY_NUMBER, Alteration.AMPLIFICATION
    if variant_class == VariantClass.COPY_NUMBER_LOSS.label:
        return ReportVariantKind.COPY_NUMBER_LOSS, Alteration.AMPLIFICATION
    return ReportVariantKind.SMALL_VARIANT, Alteration.VARIANT


def _gene_symbols(record: ClassificationModification, kind: str) -> list[str]:
    """ Both partners for a fusion, so a template can print GENE1::GENE2. The 5' partner is first -
        it is the anchor the Variant is filed under, and what the report sorts on """
    symbols = []
    if kind == ReportVariantKind.FUSION:
        variant = record.classification.variant
        gene_fusion = getattr(variant, "genefusion", None) if variant else None
        if gene_fusion:
            symbols = [gene_level_id.symbol_str for gene_level_id in gene_fusion.gene_level_ids]
    if not symbols:
        if symbol := record.get(SpecialEKeys.GENE_SYMBOL):
            symbols = [symbol]
    return symbols


def _as_float(value) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_int(value) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


@dataclass
class ReportVariant:
    """ One classification as the report prints it """
    modification: ClassificationModification
    kind: str
    alteration: str
    gene_symbol: str  # sort key; a fusion uses the 5' partner
    gene_symbols: list[str]
    splice_label: Optional[str]  # what the lab calls a splice event, eg "MET exon 14 skipping"
    tier: Optional[str]  # the raw somatic:clinical_significance value
    amp_tier: str
    tier_rank: int
    vaf: Optional[float]  # allele_frequency, as a fraction
    copy_number: Optional[int]
    fold_change: Optional[float]  # fold_change key, as the caller wrote it
    reported: bool
    sample: Optional[Sample]
    evidence: dict
    warnings: list[str] = field(default_factory=list)

    @staticmethod
    def build(record: ClassificationModification, user: User, reported: Optional[bool] = None,
              evidence: Optional[dict] = None) -> 'ReportVariant':
        kind, alteration = _kind_and_alteration(record)
        gene_symbols = _gene_symbols(record, kind)
        tier = record.get(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
        tier_label, warnings = amp_tier(record)
        if reported is None:
            reported = record.get(SpecialEKeys.VARIANT_REPORTED) != NOT_REPORTED
        return ReportVariant(
            modification=record,
            kind=kind,
            alteration=alteration,
            gene_symbol=gene_symbols[0] if gene_symbols else "",
            gene_symbols=gene_symbols,
            splice_label=record.get(SpecialEKeys.SPLICE_LABEL),
            tier=tier,
            amp_tier=tier_label,
            tier_rank=AMP_TIER_RANK.get(tier_label, UNTIERED_RANK),
            vaf=_as_float(record.get(SpecialEKeys.ALLELE_FREQUENCY)),
            copy_number=_as_int(record.get(SpecialEKeys.COPY_NUMBER)),
            fold_change=_as_float(record.get(SpecialEKeys.FOLD_CHANGE)),
            reported=reported,
            sample=record.classification.sample,
            evidence=evidence if evidence is not None else evidence_row_data(record, user),
            warnings=warnings,
        )

    @property
    def gene_label(self) -> str:
        """ 'BCR::ABL1' for a fusion, the symbol otherwise """
        if self.kind == ReportVariantKind.FUSION and len(self.gene_symbols) > 1:
            return "::".join(self.gene_symbols)
        return self.gene_symbol

    @property
    def vaf_percent(self) -> Optional[float]:
        if self.vaf is None:
            return None
        return round(self.vaf * 100, 1)

    @property
    def vaf_percent_whole(self) -> Optional[int]:
        """ The whole percent the document prints, rounded half up off the fraction the caller
            wrote - a 12.5% VAF is 13%. Rounding here rather than in the template is what stops a
            value being rounded twice, which is where half up and Python's half to even disagree """
        if self.vaf is None:
            return None
        return int((Decimal(str(self.vaf)) * 100).quantize(Decimal("1"), rounding=ROUND_HALF_UP))

    @property
    def c_hgvs(self) -> Optional[str]:
        return (self.evidence.get("c_hgvs") or {}).get("value")

    @property
    def copies_rank(self) -> int:
        """ Most copies first, except a loss, where the fewest copies is the biggest finding. A
            call with no count sorts behind the ones that have one, whichever way the kind sorts """
        if self.copy_number is None:
            return NO_COPY_NUMBER_RANK
        if self.kind == ReportVariantKind.COPY_NUMBER_LOSS:
            return self.copy_number
        return -self.copy_number

    @property
    def sort_key(self) -> tuple:
        """ Tier, then gene, then the biggest finding first - two variants in one gene print highest
            VAF first. Negated rather than reverse sorted, so the gene stays ascending """
        return (self.tier_rank, self.gene_symbol or "￿",
                -(self.vaf if self.vaf is not None else -1),
                self.copies_rank,
                -(self.fold_change if self.fold_change is not None else -1),
                self.c_hgvs or "")


@dataclass
class GeneGroup:
    """ One gene's variants and the paragraph that follows them """
    gene_symbol: str
    variants: list[ReportVariant]
    gene_summary: str = ""
    gene_summary_source: Optional[int] = None
    warnings: list[str] = field(default_factory=list)

    @property
    def classifications(self) -> list[dict]:
        """ The evidence dicts, the shape templates written against the sapath#246 context read """
        return [v.evidence for v in self.variants]


@dataclass
class TierGroup:
    """ One AMP tier of the Variant Interpretation section. `genes` holds only what is being
        reported; unreported_count is how many of the tier's variants were seen but left off """
    tier: str
    label: str
    genes: list[GeneGroup]
    unreported_count: int = 0


@dataclass
class KindGroup:
    """ One table of the Results Summary """
    kind: str
    label: str
    variants: list[ReportVariant]


@dataclass
class ReportContext:
    source_level: str
    variants: list[ReportVariant]
    kind_groups: list[KindGroup]
    tier_groups: list[TierGroup]
    gene_groups: list[GeneGroup]
    lab: Optional[Lab] = None
    user: Optional[User] = None
    case_report: Optional[Any] = None
    patient: Optional[Patient] = None
    specimen: Optional[Specimen] = None
    extractions: list[Extraction] = field(default_factory=list)
    samples: list[Sample] = field(default_factory=list)
    sequencing_runs: list[str] = field(default_factory=list)
    combined_variant_output: Optional[DragenTSO500CombinedVariantOutput] = None
    measures: dict[str, Measure] = field(default_factory=dict)
    summary: str = ""
    splice_note: Optional[str] = None
    case_values: dict = field(default_factory=dict)
    # Anything that is not a finalised report is a draft, so a preview carries the watermark too
    draft: bool = True
    generated: Optional[datetime] = None
    versions: dict = field(default_factory=dict)


def sort_report_variants(variants: list[ReportVariant]) -> list[ReportVariant]:
    """ Report order: kind, then sub-tier, then gene, then the biggest finding first """
    return sorted(variants, key=lambda v: (ReportVariantKind.ORDER.index(v.kind), *v.sort_key))


def build_report_variants(modifications: list[ClassificationModification], user: User,
                          reported_by_pk: Optional[dict[int, bool]] = None) -> list[ReportVariant]:
    """ Every classification as a ReportVariant, in report order. reported_by_pk overrides what
        the variant_reported evidence key says """
    reported_by_pk = reported_by_pk or {}
    return sort_report_variants([ReportVariant.build(record, user, reported=reported_by_pk.get(record.pk))
                                 for record in modifications])


def _gene_summary(variants: list[ReportVariant]) -> tuple[str, Optional[int], list[str]]:
    """ The gene's one paragraph. Gene level copy consensus (#1419) means these normally agree;
        where they don't the most recently modified record wins and the rest are named, so the
        scientist can bring them into line before the report is finalised """
    with_summary = [v for v in variants if (v.evidence.get(GENE_SUMMARY_KEY) or {}).get("value")]
    if not with_summary:
        return "", None, []

    by_recency = sorted(with_summary, key=lambda v: v.modification.modified, reverse=True)
    winner = by_recency[0]
    summary = winner.evidence[GENE_SUMMARY_KEY]["value"]
    disagree = [v for v in with_summary if v.evidence[GENE_SUMMARY_KEY]["value"] != summary]
    warnings = []
    if disagree:
        names = ", ".join(str(v.modification.classification.cr_lab_id) for v in disagree)
        warnings.append(f"Gene level text differs from {names} - the most recently curated was used")
    return summary, winner.modification.pk, warnings


def build_gene_groups(variants: list[ReportVariant]) -> list[GeneGroup]:
    """ The variants grouped by gene symbol, alphabetically, each with its gene level paragraph """
    by_symbol: dict[str, list[ReportVariant]] = {}
    for variant in variants:
        by_symbol.setdefault(variant.gene_symbol, []).append(variant)

    gene_groups = []
    for gene_symbol in sorted(by_symbol, key=lambda s: s or "￿"):
        in_gene = sorted(by_symbol[gene_symbol], key=lambda v: v.sort_key)
        summary, source, warnings = _gene_summary(in_gene)
        gene_groups.append(GeneGroup(gene_symbol=gene_symbol, variants=in_gene, gene_summary=summary,
                                     gene_summary_source=source, warnings=warnings))
    return gene_groups


def build_kind_groups(variants: list[ReportVariant]) -> list[KindGroup]:
    """ The Results Summary - one table per kind, in report order, of what is being reported. Every
        kind is printed whether or not the case has any, the way the tiers are: a summary that skips
        Gene Fusions leaves the reader to wonder whether they were looked for """
    return [KindGroup(kind=kind, label=ReportVariantKind.LABELS[kind],
                      variants=[v for v in variants if v.kind == kind and v.reported])
            for kind in ReportVariantKind.ORDER]


def build_splice_note(variants: list[ReportVariant]) -> Optional[str]:
    """ The sentence the document puts above the tier sections when a reported small variant has
        no protein change to print. Generated rather than typed, so the document and the JSON's
        mutations comment name the same genes; a template's own comment field overrides it """
    genes = []
    for variant in variants:
        if not variant.reported or variant.kind != ReportVariantKind.SMALL_VARIANT:
            continue
        if (variant.evidence.get(SpecialEKeys.P_HGVS) or {}).get("value") != SPLICE_NOTE_PROTEIN:
            continue
        symbol = variant.gene_symbol
        if symbol and symbol not in SPLICE_NOTE_EXCLUDED_GENES and symbol not in genes:
            genes.append(symbol)

    if not genes:
        return None
    if len(genes) == 1:
        return (f"Note that the {genes[0]} variant with protein designation {SPLICE_NOTE_PROTEIN} "
                "is an intronic variant that is predicted to disrupt splicing, see below.")
    named = f"{', '.join(genes[:-1])} and {genes[-1]}"
    return (f"Note that the {named} variants with protein designation {SPLICE_NOTE_PROTEIN} are "
            "intronic variants that are predicted to disrupt splicing, see below.")


def build_tier_groups(variants: list[ReportVariant]) -> list[TierGroup]:
    """ The Variant Interpretation - tier, then gene, with the kinds interleaved so a Tier IIC
        amplification prints under Tier II beside the Tier IIC small variants.

        A tier outside ALWAYS_PRINTED_TIERS with nothing being reported from it is left out
        entirely rather than printed empty. That is what keeps Tier IV off a somatic report: a
        Tier IV variant is a not-reported record and `reported` already says so, so the section
        would be a heading over "no reportable variants detected" on every case that has one """
    tiers = list(TIER_GROUP_ORDER)
    for variant in variants:
        tier = variant.tier or UNTIERED
        if tier not in tiers:
            tiers.append(tier)
    if UNTIERED not in tiers and any(not v.tier for v in variants):
        tiers.append(UNTIERED)

    tier_groups = []
    for tier in tiers:
        in_tier = [v for v in variants if (v.tier or UNTIERED) == tier]
        reported = [v for v in in_tier if v.reported]
        if not reported and tier not in ALWAYS_PRINTED_TIERS:
            continue
        label = TIER_GROUP_LABELS.get(tier, UNTIERED_LABEL)
        tier_groups.append(TierGroup(tier=tier, label=label, genes=build_gene_groups(reported),
                                     unreported_count=len(in_tier) - len(reported)))
    return tier_groups


def case_specimen(source_level: str, source) -> Optional[Specimen]:
    """ The specimen a case is about - the level it is reported from, or the one its arm came off.
        A patient or sample case has none, so it is reported without measures """
    if source_level == SampleSourceLevel.SPECIMEN:
        return source
    if source_level == SampleSourceLevel.EXTRACTION:
        return source.specimen
    return None


def case_combined_variant_output(samples: list[Sample],
                                 specimen: Optional[Specimen]) -> Optional[DragenTSO500CombinedVariantOutput]:
    """ The one analysis a case's TMB, MSI and GIS come off - the one whose arm is one of the case's samples,
        directly or through the sequencing sample it came off. Where none is linked yet, the newest
        analysis claiming the specimen, as specimen_library_qc picks """
    by_output = [F("output_datetime").desc(nulls_last=True), "-pk"]
    if samples:
        arm_q = Q()
        for arm in ("dna", "rna"):
            arm_q |= Q(**{f"{arm}_sample__in": samples})
            arm_q |= Q(**{f"{arm}_sequencing_sample__samplefromsequencingsample__sample__in": samples})
        if cvo := DragenTSO500CombinedVariantOutput.objects.filter(arm_q).distinct().order_by(*by_output).first():
            return cvo
    if specimen is None:
        return None
    return DragenTSO500CombinedVariantOutput.objects.filter(specimen=specimen).order_by(*by_output).first()


def case_measures(cvo: Optional[DragenTSO500CombinedVariantOutput],
                  specimen: Optional[Specimen]) -> dict[str, Measure]:
    """ The case's measures by context key - what the report prints, and what the build form starts its
        assay flags from: the analysis' numbers, and the pathologist's tumour content off the specimen """
    measures = {}
    if cvo:
        for key, (column, unit) in CVO_MEASURE_CONTEXT_KEYS.items():
            value = getattr(cvo, column)
            if value is None:
                continue
            measure_call = cvo.measure_call(key)
            measures[key] = Measure(value=value, unit=unit,
                                    call=measure_call.call if measure_call else None,
                                    threshold=measure_call.threshold if measure_call else None,
                                    threshold_source=measure_call.threshold_source if measure_call else None,
                                    method=cvo.method)
    if specimen:
        for specimen_measure in SpecimenMeasure.objects.filter(specimen=specimen):
            if key := MEASURE_CONTEXT_KEYS.get(specimen_measure.measure_type):
                measures[key] = Measure(value=specimen_measure.value, unit=specimen_measure.unit,
                                        call=specimen_measure.call, threshold=None, threshold_source=None,
                                        method=specimen_measure.method)
    return measures


def specimen_library_qc(specimen: Optional[Specimen]) -> dict[str, LibraryQC]:
    """ What the caller's QC said about the specimen's libraries, by context key - the newest run
        per category wins, since a repeat sequencing supersedes the QC of the one it replaced """
    if specimen is None:
        return {}
    library_qc = {}
    qs = LibraryQC.objects.filter(specimen=specimen) \
        .order_by(F("measured_date").asc(nulls_first=True), "pk")
    for qc in qs:
        if key := LIBRARY_QC_CONTEXT_KEYS.get(qc.category):
            library_qc[key] = qc
    return library_qc


def _sequencing_runs(samples: list[Sample]) -> list[str]:
    """ The runs the case's samples came off, where seqauto linked them """
    names = []
    for sample in samples:
        try:
            name = sample.samplefromsequencingsample.sequencing_sample.sequencing_run.name
        except ObjectDoesNotExist:
            continue  # seqauto isn't run everywhere, and a sample can be loaded straight from a VCF
        if name not in names:
            names.append(name)
    return names


DRAGEN_VERSION_PATTERN = re.compile(r"^##(?:DRAGEN(?:Version|CommandLine)|source)=.*$", re.MULTILINE)


def _caller_versions(samples: list[Sample]) -> list[str]:
    """ What called the VCFs - the Method section names it. Straight out of the header, because a
        header line is what the lab can point at when asked which pipeline produced a report """
    versions = []
    for vcf in {sample.vcf for sample in samples}:
        for line in DRAGEN_VERSION_PATTERN.findall(vcf.header or ""):
            line = line.strip()
            if line not in versions:
                versions.append(line)
    return versions


def build_versions(samples: list[Sample]) -> dict:
    annotation = {}
    for genome_build in {sample.genome_build for sample in samples}:
        if vav := VariantAnnotationVersion.latest(genome_build):
            annotation[genome_build.name] = str(vav)
    return {
        "variantgrid": get_cached_project_git_hash(),
        "annotation": annotation,
        "callers": _caller_versions(samples),
    }


def build_report_context(user: User, source_level: str, source,
                         modifications: list[ClassificationModification],
                         lab: Optional[Lab] = None,
                         reported_by_pk: Optional[dict[int, bool]] = None,
                         summary: str = "", case_values: Optional[dict] = None,
                         case_report=None) -> ReportContext:
    """ The whole case: who it is, what was sequenced, what was measured and what was classified """
    variants = build_report_variants(modifications, user, reported_by_pk=reported_by_pk)
    samples = get_sample_group(user, source_level, source).samples
    patient = get_patient_for_source(source_level, source)
    specimen = case_specimen(source_level, source)
    combined_variant_output = case_combined_variant_output(samples, specimen)
    extractions = []
    if source_level == SampleSourceLevel.EXTRACTION:
        extractions = [source]
    elif specimen is not None:
        extractions = list(specimen.extraction_set.order_by("pk"))

    return ReportContext(
        source_level=source_level,
        variants=variants,
        kind_groups=build_kind_groups(variants),
        tier_groups=build_tier_groups(variants),
        gene_groups=build_gene_groups(variants),
        lab=lab,
        user=user,
        case_report=case_report,
        patient=patient,
        specimen=specimen,
        extractions=extractions,
        samples=samples,
        sequencing_runs=_sequencing_runs(samples),
        combined_variant_output=combined_variant_output,
        measures=case_measures(combined_variant_output, specimen),
        summary=summary,
        splice_note=build_splice_note(variants),
        case_values=case_values or {},
        draft=case_report is None or case_report.status == CaseReportStatus.DRAFT,
        generated=timezone.now(),
        versions=build_versions(samples),
    )


def _model_as_dict(obj, fields: list[str]) -> Optional[dict]:
    """ Model instances reach a template as a small dict, so the same structure serialises straight
        into context_snapshot and the JSON output """
    if obj is None:
        return None
    as_dict = {"pk": obj.pk, "str": str(obj)}
    for name in fields:
        value = getattr(obj, name, None)
        as_dict[name] = value if not hasattr(value, "pk") else str(value)
    return as_dict


# What a template and the JSON record read off the report itself - the LIS fields can be entered
# after a report is final, so a rebuild refreshes this block over the stored snapshot
CASE_REPORT_DICT_FIELDS = ["status", "external_report_id", "report_date"]


def case_report_as_dict(case_report) -> Optional[dict]:
    return _model_as_dict(case_report, CASE_REPORT_DICT_FIELDS)


def _variant_as_dict(variant: ReportVariant) -> dict:
    return {
        "classification_modification_id": variant.modification.pk,
        "classification_id": variant.modification.classification_id,
        "kind": variant.kind,
        "alteration": variant.alteration,
        "gene_symbol": variant.gene_symbol,
        "gene_symbols": variant.gene_symbols,
        "gene_label": variant.gene_label,
        "splice_label": variant.splice_label,
        "tier": variant.tier,
        "amp_tier": variant.amp_tier,
        "tier_rank": variant.tier_rank,
        "vaf": variant.vaf,
        "vaf_percent": variant.vaf_percent,
        "vaf_percent_whole": variant.vaf_percent_whole,
        "copy_number": variant.copy_number,
        "fold_change": variant.fold_change,
        "reported": variant.reported,
        "sample": _model_as_dict(variant.sample, ["name"]),
        "evidence": variant.evidence,
        "warnings": variant.warnings,
    }


def _gene_group_as_dict(gene_group: GeneGroup) -> dict:
    return {
        "gene_symbol": gene_group.gene_symbol,
        "variants": [_variant_as_dict(v) for v in gene_group.variants],
        # The sapath#246 shape - a template written against it keeps working
        "classifications": gene_group.classifications,
        "gene_summary": gene_group.gene_summary,
        "gene_summary_source": gene_group.gene_summary_source,
        "warnings": gene_group.warnings,
    }


def _combined_variant_output_as_dict(cvo: Optional[DragenTSO500CombinedVariantOutput]) -> Optional[dict]:
    """ Which analysis the report's TMB, MSI and GIS came off """
    return _model_as_dict(cvo, ["sequencing_run_name", "pair_id", "dna_sample_name", "rna_sample_name",
                                "output_datetime", "module_version", "pipeline_version"])


def context_as_dict(report_context: ReportContext) -> dict:
    """ What the templates render over, and what CaseReport.context_snapshot stores - so a report
        can be re-rendered later without re-deriving numbers whose source rows may have changed """
    return {
        "source_level": report_context.source_level,
        "case_report": case_report_as_dict(report_context.case_report),
        "patient": _model_as_dict(report_context.patient,
                                  ["patient_code", "date_of_birth", "sex"]),
        "specimen": _model_as_dict(report_context.specimen,
                                   ["reference_id", "collection_date", "received_date"]),
        "extractions": [_model_as_dict(e, ["reference_id", "nucleic_acid_source"])
                        for e in report_context.extractions],
        "samples": [_model_as_dict(s, ["name"]) for s in report_context.samples],
        "sequencing_runs": report_context.sequencing_runs,
        "analysis": _combined_variant_output_as_dict(report_context.combined_variant_output),
        "measures": {key: {"value": measure.value, "unit": measure.unit, "call": measure.call,
                           "threshold": measure.threshold, "method": measure.method,
                           "str": measure.value_description}
                     for key, measure in report_context.measures.items()},
        "variants": [_variant_as_dict(v) for v in report_context.variants],
        "kind_groups": [{"kind": kg.kind, "label": kg.label,
                         "variants": [_variant_as_dict(v) for v in kg.variants]}
                        for kg in report_context.kind_groups],
        "tier_groups": [{"tier": tg.tier, "label": tg.label,
                         "genes": [_gene_group_as_dict(g) for g in tg.genes],
                         "unreported_count": tg.unreported_count}
                        for tg in report_context.tier_groups],
        "gene_groups": [_gene_group_as_dict(g) for g in report_context.gene_groups],
        "summary": report_context.summary,
        "splice_note": report_context.splice_note,
        "case_values": report_context.case_values,
        "draft": report_context.draft,
        "lab": _model_as_dict(report_context.lab, ["name"]),
        "user": _model_as_dict(report_context.user, ["username"]),
        "generated": report_context.generated,
        "versions": report_context.versions,
    }
