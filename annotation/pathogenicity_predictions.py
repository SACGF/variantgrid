from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from typing import Optional


class RawScoreDirection(Enum):
    """ Which way a raw score runs, so the DamageNode slider filters on the damaging side of it """
    HIGHER = "higher"        # higher = more damaging (every dbNSFP score): keep score >= slider
    LOWER = "lower"          # lower = more damaging (popEVE log-likelihood ratios): keep score <= slider
    MAGNITUDE = "magnitude"  # signed, strength is |score| (PromoterAI): keep |score| >= slider

    @property
    def comparison(self) -> str:
        """ How the slider value reads against the score, for labels and the thresholds page """
        return {
            RawScoreDirection.HIGHER: "≥",
            RawScoreDirection.LOWER: "≤",
            RawScoreDirection.MAGNITUDE: "|score| ≥",
        }[self]


@dataclass(frozen=True)
class PathogenicityTool:
    name: str                                          # display name
    rankscore_field: Optional[str]                     # VariantAnnotation rankscore col, may be None
    raw_field: Optional[str]                           # VariantAnnotation raw-score col, may be None
    pred_field: Optional[str]                          # VariantAnnotation pred char col, may be None
    raw_min: Optional[float]                           # slider min
    raw_max: Optional[float]                           # slider max
    raw_step: Optional[float]                          # slider step
    raw_max_benign_threshold: Optional[float]          # BP4-supporting upper bound (colour band lower)
    raw_pathogenic_threshold: Optional[float]          # PP3-supporting lower bound (colour band upper, raw count cutoff)
    pred_pathogenic_values: tuple[str, ...] = ()       # categorical pred values counted as pathogenic
    raw_direction: RawScoreDirection = RawScoreDirection.HIGHER
    node_label: str = ""                               # DamageNode slider label, where the auto field label would mislead
    source: str = ""                                   # short citation tag (table cell / colour-band tooltip)
    source_detail: str = ""                            # full citation, shown on the Pathogenicity Thresholds page
    source_url: str = ""                               # DOI / URL for the citation
    note: str = ""                                     # caveat shown on the Pathogenicity Thresholds page (e.g. threshold not clinically recommended)

    @property
    def node_threshold_field(self) -> str:
        """ DamageNode field holding this tool's slider value - named for the side it keeps """
        suffix = "max" if self.raw_direction == RawScoreDirection.LOWER else "min"
        return f"{self.raw_field}_{suffix}"


# Pejaver V, Byrne AB, Feng B-J, et al. Calibration of computational tools for missense variant
# pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria.
# Am J Hum Genet. 2022;109(12):2163-2177. DOI: 10.1016/j.ajhg.2022.10.013 (PMC9748256). Table 2.
#
# Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools
# expands ClinGen recommendation options for variant classification with PP3/BP4 criteria.
# Genet Med. 2025;27(3):101353. DOI: 10.1016/j.gim.2025.101353 (PMC11429929). Table 2.
# Full citations + DOIs, surfaced on the Pathogenicity Thresholds page. (detail, url) pairs shared by tools.
_PEJAVER_2022 = (
    "Pejaver V, Byrne AB, Feng B-J, et al. Calibration of computational tools for missense variant "
    "pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. "
    "Am J Hum Genet. 2022;109(12):2163-2177.",
    "https://doi.org/10.1016/j.ajhg.2022.10.013",
)
_BERGQUIST_2024 = (
    "Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools expands "
    "ClinGen recommendation options for variant classification with PP3/BP4 criteria. "
    "Genet Med. 2025;27(3):101353.",
    "https://doi.org/10.1016/j.gim.2025.101353",
)
_ALIREZAIE_2018 = (
    "Alirezaie N, Kernohan KD, Hartley T, Majewski J, Hocking TD. ClinPred: Prediction Tool to Identify "
    "Disease-Relevant Nonsynonymous Single-Nucleotide Variants. Am J Hum Genet. 2018;103(4):474-483. "
    "Author cutoff 0.5 (no ClinGen calibration).",
    "https://doi.org/10.1016/j.ajhg.2018.08.005",
)
_LI_2022 = (
    "Li C, Zhi D, Wang K, Liu X. MetaRNN: differentiating rare pathogenic and rare benign missense SNVs "
    "and InDels using deep learning. Genome Med. 2022;14(1):115. Author cutoff 0.5 (no ClinGen calibration).",
    "https://doi.org/10.1186/s13073-022-01120-z",
)
_WU_2021 = (
    "Wu Y, Li R, Sun S, et al. Improved pathogenicity prediction for rare human missense variants. "
    "Am J Hum Genet. 2021;108(10):1891-1906. No ClinGen calibration.",
    "https://doi.org/10.1016/j.ajhg.2021.08.012",
)


TOOLS: tuple[PathogenicityTool, ...] = (
    # AlphaMissense - Cheng 2023, Science (DOI 10.1126/science.adg7492). dbNSFP 5.3a readme field 135.
    # Bergquist 2024 PP3-supporting band [0.170, 0.791]; BP4-supporting upper 0.169. Note the
    # supporting band is essentially adjacent (0.169 vs 0.170) - colour gap will be tiny.
    PathogenicityTool(
        name="AlphaMissense",
        rankscore_field="alphamissense_rankscore",
        raw_field="alphamissense_score",
        pred_field="alphamissense_pred",
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=0.169,
        raw_pathogenic_threshold=0.170,
        pred_pathogenic_values=("p", "P"),
        source="Bergquist 2024", source_detail=_BERGQUIST_2024[0], source_url=_BERGQUIST_2024[1],
    ),
    # BayesDel_noAF - Feng 2017, Hum Mutat (DOI 10.1002/humu.23158). dbNSFP 5.3a readme fields 113-115.
    # Range -1.31914 to 0.840878. Pejaver 2022 PP3-supporting band [0.13, 0.27); BP4-supporting upper -0.18.
    PathogenicityTool(
        name="BayesDel (no AF)",
        rankscore_field="bayesdel_noaf_rankscore",
        raw_field="bayesdel_noaf_score",
        pred_field=None,
        raw_min=-1.32, raw_max=0.85, raw_step=0.05,
        raw_max_benign_threshold=-0.18,
        raw_pathogenic_threshold=0.13,
        source="Pejaver 2022", source_detail=_PEJAVER_2022[0], source_url=_PEJAVER_2022[1],
    ),
    # CADD - Kircher 2014, Nat Genet (DOI 10.1038/ng.2892). claude/pdfs/nihms555958.pdf.
    # Pejaver 2022 calibrated CADD_phred PP3-supporting band [25.3, 28.1); BP4-supporting upper 22.7.
    # We slider on cadd_phred (the user-facing scale). cadd_raw is not exposed as a slider.
    PathogenicityTool(
        name="CADD (phred)",
        rankscore_field="cadd_raw_rankscore",
        raw_field="cadd_phred",
        pred_field=None,
        raw_min=0.0, raw_max=99.0, raw_step=1.0,
        raw_max_benign_threshold=22.7,
        raw_pathogenic_threshold=25.3,
        source="Pejaver 2022", source_detail=_PEJAVER_2022[0], source_url=_PEJAVER_2022[1],
    ),
    # ClinPred - Alirezaie 2018, AJHG (DOI 10.1016/j.ajhg.2018.08.005). dbNSFP 5.3a readme field 118.
    # No ClinGen calibration; use author cutoff 0.5. Single threshold (no BP4 band).
    PathogenicityTool(
        name="ClinPred",
        rankscore_field="clinpred_rankscore",
        raw_field="clinpred_score",
        pred_field="clinpred_pred",
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=0.5,
        pred_pathogenic_values=("D",),
        source="dbNSFP readme 118", source_detail=_ALIREZAIE_2018[0], source_url=_ALIREZAIE_2018[1],
    ),
    # EVE - Frazer 2021, Nature (DOI 10.1038/s41586-021-04043-8). GRCh38, VEP >= 116; higher = more pathogenic.
    # Bergquist 2024 measured EVE supporting thresholds (pathogenic >= 0.684, benign <= 0.137) but explicitly
    # declined to recommend them clinically: EVE scored only ~half of their benign/likely-benign calibration
    # variants, so potential sampling bias left them unable to endorse the cutoffs. We therefore list EVE with
    # no active PP3/BP4 cutoff (no colour band, not counted in predictions_num_pathogenic) and say so via `note`.
    PathogenicityTool(
        name="EVE",
        rankscore_field=None,
        raw_field="eve_score",
        pred_field="eve_class",
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=None,
        node_label="EVE score",
        source="Bergquist 2024", source_detail=_BERGQUIST_2024[0], source_url=_BERGQUIST_2024[1],
        note="Bergquist 2024 measured supporting thresholds (pathogenic ≥ 0.684, benign ≤ 0.137) but "
             "declined to recommend them for clinical use due to potential sampling bias, so no "
             "ClinGen-calibrated PP3/BP4 cutoff is applied. EVE_CLASS uses the plugin's own "
             "Benign / Uncertain / Pathogenic call at a 75% uncertain threshold.",
    ),
    # popEVE - Orenbuch 2023 (population-adjusted EVE). GRCh38, VEP >= 116. Scores are
    # log-likelihood ratios, so they run the other way to every dbNSFP score: more negative =
    # more damaging, and the slider keeps everything at or below it. No ClinGen calibration.
    PathogenicityTool(
        name="popEVE",
        rankscore_field=None,
        raw_field="popeve_score",
        pred_field=None,
        raw_min=-15.0, raw_max=2.0, raw_step=0.5,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=None,
        raw_direction=RawScoreDirection.LOWER,
        node_label="popEVE score",
        source="popEVE (no calibration)",
        source_detail="Orenbuch R, Kollasch AW, Spinner HD, et al. Deep generative modeling of the human "
                      "proteome reveals over a hundred novel genes involved in rare genetic disorders. "
                      "medRxiv 2023.11.27.23299062. Population-adjusted EVE; no ClinGen calibration.",
        source_url="https://doi.org/10.1101/2023.11.27.23299062",
        note="Log-likelihood ratio scale - more negative is more damaging, so the slider keeps scores "
             "at or below the chosen value. No ClinGen-calibrated PP3/BP4 cutoff.",
    ),
    # MetaRNN - Li 2022, Genome Med (DOI 10.1186/s13073-022-01120-z). dbNSFP 5.3a readme field 82.
    # No ClinGen calibration; use author cutoff 0.5. Single threshold.
    PathogenicityTool(
        name="MetaRNN",
        rankscore_field=None,
        raw_field="metarnn_score",
        pred_field="metarnn_pred",
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=0.5,
        pred_pathogenic_values=("D",),
        source="dbNSFP readme 82", source_detail=_LI_2022[0], source_url=_LI_2022[1],
    ),
    # MPC - Samocha 2017, bioRxiv 148353 (DOI 10.1101/148353). claude/pdfs/148353v1.full.pdf.
    # Range 0-5 (page 14). Pejaver 2022 PP3-supporting band [1.360, 1.828). No BP4 band defined.
    PathogenicityTool(
        name="MPC",
        rankscore_field=None,
        raw_field="mpc_score",
        pred_field=None,
        raw_min=0.0, raw_max=5.0, raw_step=0.1,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=1.360,
        source="Pejaver 2022", source_detail=_PEJAVER_2022[0], source_url=_PEJAVER_2022[1],
    ),
    # MutPred2 - Pejaver 2020, Nat Commun (model); Pejaver 2022, AJHG (calibration). dbNSFP 5.3a readme field 90.
    # PP3-supporting band [0.737, 0.829); BP4-supporting upper 0.391.
    PathogenicityTool(
        name="MutPred2",
        rankscore_field=None,
        raw_field="mutpred2_score",
        pred_field=None,
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=0.391,
        raw_pathogenic_threshold=0.737,
        source="Pejaver 2022", source_detail=_PEJAVER_2022[0], source_url=_PEJAVER_2022[1],
    ),
    # PrimateAI - Sundaram 2018, Nat Genet (DOI 10.1038/s41588-018-0167-z). dbNSFP 5.3a readme fields 104, 106.
    # Pejaver 2022 PP3-supporting band [0.790, 0.867); BP4-supporting upper 0.483.
    PathogenicityTool(
        name="PrimateAI",
        rankscore_field=None,
        raw_field="primateai_score",
        pred_field="primateai_pred",
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=0.483,
        raw_pathogenic_threshold=0.790,
        pred_pathogenic_values=("D",),
        source="Pejaver 2022", source_detail=_PEJAVER_2022[0], source_url=_PEJAVER_2022[1],
    ),
    # REVEL - Ioannidis 2016, AJHG (DOI 10.1016/j.ajhg.2016.08.016). claude/pdfs/main.pdf.
    # Pejaver 2022 PP3-supporting band [0.644, 0.773); BP4-supporting upper 0.290.
    PathogenicityTool(
        name="REVEL",
        rankscore_field="revel_rankscore",
        raw_field="revel_score",
        pred_field=None,
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=0.290,
        raw_pathogenic_threshold=0.644,
        source="Pejaver 2022", source_detail=_PEJAVER_2022[0], source_url=_PEJAVER_2022[1],
    ),
    # VARITY_R - Wu 2021, AJHG (DOI 10.1016/j.ajhg.2021.08.012). claude/pdfs/1-s2.0-S0002929721003207-mainext.pdf.
    # Bergquist 2024 PP3-supporting band [0.252, 0.674]; BP4-supporting upper 0.251. Adjacent-band note
    # as for AlphaMissense.
    PathogenicityTool(
        name="VARITY_R",
        rankscore_field=None,
        raw_field="varity_r_score",
        pred_field=None,
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=0.251,
        raw_pathogenic_threshold=0.252,
        source="Bergquist 2024", source_detail=_BERGQUIST_2024[0], source_url=_BERGQUIST_2024[1],
    ),
    # VARITY_ER - Wu 2021, AJHG (same paper). No ClinGen calibration. Slider only; not in raw count.
    PathogenicityTool(
        name="VARITY_ER",
        rankscore_field=None,
        raw_field="varity_er_score",
        pred_field=None,
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=None,
        source="Wu 2021 (no calibration)", source_detail=_WU_2021[0], source_url=_WU_2021[1],
    ),
    # VEST4 - Carter 2013, BMC Genomics (DOI 10.1186/1471-2164-14-S3-S3). claude/pdfs/1471-2164-14-S3-S3.pdf.
    # Pejaver 2022 PP3-supporting band [0.764, 0.861); BP4-supporting upper 0.449.
    PathogenicityTool(
        name="VEST4",
        rankscore_field="vest4_rankscore",
        raw_field="vest4_score",
        pred_field=None,
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=0.449,
        raw_pathogenic_threshold=0.764,
        source="Pejaver 2022", source_detail=_PEJAVER_2022[0], source_url=_PEJAVER_2022[1],
    ),

    # ---- Not missense pathogenicity predictors: no PP3/BP4 cutoff, no colour band, and they sit
    # outside predictions_num_pathogenic. Offered as Effect node filters only (#1808). ----
    # ProtVar ddG - protein stability change in kcal/mol from the ProtVar plugin. Not a pathogenicity
    # predictor: it says whether a substitution destabilises the fold. The plugin header calls
    # ddG > 2 "likely to be unstabilising"; higher = more destabilising.
    PathogenicityTool(
        name="ProtVar ddG stability",
        rankscore_field=None,
        raw_field="protvar_stability",
        pred_field=None,
        raw_min=-10.0, raw_max=10.0, raw_step=0.5,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=None,
        node_label="ProtVar ddG stability",
        source="ProtVar (no calibration)",
        source_detail="Stephenson JD, Totoo P, Burke DF, et al. ProtVar: mapping and contextualizing human "
                      "missense variation. Nucleic Acids Res. 2024;52(W1):W140-W147. Stability is ddG in "
                      "kcal/mol, not a pathogenicity score; no ClinGen calibration.",
        source_url="https://doi.org/10.1093/nar/gkae413",
        note="Protein stability change (ddG, kcal/mol), not a pathogenicity prediction. The plugin calls "
             "ddG > 2 likely destabilising. No ClinGen-calibrated PP3/BP4 cutoff.",
    ),
    # PromoterAI - Illumina. Promoter variants' effect on expression, so unlike every other tool here it
    # is not missense. The score is signed (predicted expression up/down) and it is the magnitude that
    # says how strong the effect is - matching how the importer picks the strongest TSS entry.
    PathogenicityTool(
        name="PromoterAI",
        rankscore_field=None,
        raw_field="promoter_ai_score",
        pred_field=None,
        raw_min=0.0, raw_max=1.0, raw_step=0.05,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=None,
        raw_direction=RawScoreDirection.MAGNITUDE,
        node_label="PromoterAI",
        source="PromoterAI (no calibration)",
        source_detail="Jaganathan K, Panagiotopoulou SK, McRae JF, et al. PromoterAI: deep learning "
                      "predicts promoter variant effects on gene expression (Illumina). Signed score; "
                      "no ClinGen calibration.",
        source_url="https://github.com/Illumina/PromoterAI",
        note="Promoter (non-coding) expression effect, not a missense pathogenicity prediction. The score "
             "is signed for direction of expression change, so the slider filters on |score|. "
             "No ClinGen-calibrated PP3/BP4 cutoff.",
    ),
)


TOOLS_BY_RAW_FIELD = {t.raw_field: t for t in TOOLS if t.raw_field}
TOOLS_BY_PRED_FIELD = {t.pred_field: t for t in TOOLS if t.pred_field}


def raw_score_pathogenic_funcs() -> dict[str, Callable]:
    """Drives DamageNode raw-score filtering and predictions_num_pathogenic at v4.
    Includes only tools with a raw_pathogenic_threshold (BayesDel_noAF, CADD_phred, ClinPred,
    MetaRNN, MPC, MutPred2, PrimateAI, REVEL, VARITY_R, VEST4, AlphaMissense)."""
    funcs = {}
    for t in TOOLS:
        if t.raw_field and t.raw_pathogenic_threshold is not None:
            threshold = t.raw_pathogenic_threshold
            funcs[t.raw_field] = lambda v, th=threshold: float(v) >= th
    return funcs


def pred_pathogenic_funcs() -> dict[str, Callable]:
    """Pred-field categorical contributions to predictions_num_pathogenic at v4.
    AlphaMissense uses pred {p, P} (LP/P) per Bergquist 2024 / Cheng 2023."""
    funcs = {}
    for t in TOOLS:
        if t.pred_field and t.pred_pathogenic_values:
            values = t.pred_pathogenic_values
            funcs[t.pred_field] = lambda v, vs=values: v in vs
    return funcs
