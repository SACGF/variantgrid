from dataclasses import dataclass
from typing import Callable, Optional


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
    source: str = ""                                   # citation tag, expanded in the row's # comment


# Pejaver V, Byrne AB, Feng B-J, et al. Calibration of computational tools for missense variant
# pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria.
# Am J Hum Genet. 2022;109(12):2163-2177. DOI: 10.1016/j.ajhg.2022.10.013 (PMC9748256). Table 2.
#
# Bergquist T, Stenton SL, Nadeau EAW, et al. Calibration of additional computational tools
# expands ClinGen recommendation options for variant classification with PP3/BP4 criteria.
# Genet Med. 2025;27(3):101353. DOI: 10.1016/j.gim.2025.101353 (PMC11429929). Table 2.
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
        source="Bergquist 2024",
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
        source="Pejaver 2022",
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
        source="Pejaver 2022",
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
        source="dbNSFP readme 118",
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
        source="dbNSFP readme 82",
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
        source="Pejaver 2022",
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
        source="Pejaver 2022",
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
        source="Pejaver 2022",
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
        source="Pejaver 2022",
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
        source="Bergquist 2024",
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
        source="Wu 2021 (no calibration)",
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
        source="Pejaver 2022",
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
