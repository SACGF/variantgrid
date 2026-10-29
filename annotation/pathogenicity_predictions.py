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
    # CADD - Kircher 2014, Nat Genet (DOI 10.1038/ng.2892).
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
    # MPC - Samocha 2017, bioRxiv 148353 (DOI 10.1101/148353).
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
    # REVEL - Ioannidis 2016, AJHG (DOI 10.1016/j.ajhg.2016.08.016).
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
    # VARITY_R - Wu 2021, AJHG (DOI 10.1016/j.ajhg.2021.08.012).
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
    # VEST4 - Carter 2013, BMC Genomics (DOI 10.1186/1471-2164-14-S3-S3).
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

    # MetaLR - Dong 2015, Hum Mol Genet (DOI 10.1093/hmg/ddu733). dbNSFP ensemble score; we only
    # store the rankscore, so there is no raw slider - it bands off the site-wide rankscore settings.
    PathogenicityTool(
        name="MetaLR",
        rankscore_field="metalr_rankscore",
        raw_field=None,
        pred_field=None,
        raw_min=None, raw_max=None, raw_step=None,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=None,
        source="Dong 2015",
        source_detail="Dong C, Wei P, Jian X, et al. Comparison and integration of deleteriousness "
                      "prediction methods for nonsynonymous SNVs in whole exome sequencing studies. "
                      "Hum Mol Genet. 2015;24(8):2125-2137.",
        source_url="https://doi.org/10.1093/hmg/ddu733",
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
    # AVI (AlphaGenome Variant Impact) - Cheng 2026, medRxiv. GRCh38 SNVs, VEP >= 116. The PHRED scale is a
    # genome-wide rank (20 = top 1% of all SNVs, what their rare-disease triage used), not a pathogenicity call, and
    # the authors set no cutoff. Like CADD it isn't calibrated for non-coding variants, so it has no colour band and
    # stays out of predictions_num_pathogenic (#2561).
    PathogenicityTool(
        name="AVI (AlphaGenome Variant Impact)",
        rankscore_field=None,
        raw_field="avi_phred",
        pred_field=None,
        raw_min=0.0, raw_max=80.0, raw_step=1.0,
        raw_max_benign_threshold=None,
        raw_pathogenic_threshold=None,
        node_label="AVI (phred)",
        source="Cheng 2026 (no calibration)",
        source_detail="Cheng J, et al. AlphaGenome Atlas: in silico mutagenesis of the entire human genome "
                      "improves prioritization and interpretation of non-coding variants. medRxiv 2026 "
                      "(preprint, not peer reviewed). PHRED is a genome-wide rank: 10 = top 10%, 20 = top 1%.",
        source_url="https://doi.org/10.64898/2026.09.16.26363192",
        note="Coding and non-coding SNVs (GRCh38). PHRED 20 is the genome-wide top 1%, which the authors used "
             "to triage de novo variants in rare-disease cases - a ranking choice, not a calibrated "
             "pathogenicity threshold. The authors recommend region- or application-aware thresholds. "
             "No ClinGen-calibrated PP3/BP4 cutoff.",
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



@dataclass(frozen=True)
class ColourBand:
    """ Colour band for a variant details field that isn't a missense pathogenicity predictor, so it
        has no PP3/BP4 cutoff and takes no part in the damage count. Same green / amber / red bands as
        TOOLS, listed on the Pathogenicity Thresholds page so every coloured number has a stated source.
        max_benign_threshold None colours only the damaging end - the rest of the scale stays plain. """
    name: str
    fields: tuple[str, ...]
    max_benign_threshold: Optional[float]
    pathogenic_threshold: Optional[float]
    direction: RawScoreDirection = RawScoreDirection.HIGHER
    source: str = ""
    source_detail: str = ""
    source_url: str = ""
    note: str = ""

    @property
    def magnitude(self) -> bool:
        """ Band on |score| - the field is signed but only the size of the change matters """
        return self.direction == RawScoreDirection.MAGNITUDE


_VG_DEFAULT = "VariantGrid default"

COLOUR_BANDS: tuple[ColourBand, ...] = (
    # SpliceAI - Jaganathan 2019 gives 0.2 (high recall), 0.5 (recommended) and 0.8 (high precision)
    # delta score cutoffs; we band on the first two.
    ColourBand(
        name="SpliceAI delta scores",
        fields=("spliceai_pred_ds_ag", "spliceai_pred_ds_al", "spliceai_pred_ds_dg", "spliceai_pred_ds_dl"),
        max_benign_threshold=0.2,
        pathogenic_threshold=0.5,
        source="Jaganathan 2019",
        source_detail="Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting Splicing "
                      "from Primary Sequence with Deep Learning. Cell. 2019;176(3):535-548. Delta score "
                      "cutoffs 0.2 (high recall), 0.5 (recommended), 0.8 (high precision).",
        source_url="https://doi.org/10.1016/j.cell.2018.12.015",
    ),
    # dbscSNV - Jian 2014 recommends 0.6 for both ensemble scores. No benign cutoff published.
    ColourBand(
        name="dbscSNV (ada / rf)",
        fields=("dbscsnv_ada_score", "dbscsnv_rf_score"),
        max_benign_threshold=None,
        pathogenic_threshold=0.6,
        source="Jian 2014",
        source_detail="Jian X, Boerwinkle E, Liu X. In silico prediction of splice-altering single "
                      "nucleotide variants in the human genome. Nucleic Acids Res. 2014;42(22):13534-13544. "
                      "Both ensemble scores are called splice-altering at 0.6.",
        source_url="https://doi.org/10.1093/nar/gku1206",
    ),
    # MaxEntScan - the score is a signed % change against the reference site, and it is the size of the
    # change that matters, in either direction (a weakened native site or a strengthened cryptic one).
    ColourBand(
        name="MaxEntScan % diff vs ref",
        fields=("maxentscan_percent_diff_ref",),
        max_benign_threshold=70.0,
        pathogenic_threshold=150.0,
        direction=RawScoreDirection.MAGNITUDE,
        source=_VG_DEFAULT,
        source_detail="Yeo G, Burge CB. Maximum entropy modeling of short sequence motifs with applications "
                      "to RNA splicing signals. J Comput Biol. 2004;11(2-3):377-394. The band sizes are a "
                      "VariantGrid display choice - the paper sets no % change cutoff.",
        source_url="https://doi.org/10.1089/1066527041410418",
        note="Signed % change against the reference splice site, so the band is on |% diff|.",
    ),
    # Conservation - a conserved base is worth flagging, but a poorly conserved one is not benign
    # evidence on its own, so these colour the conserved end only.
    ColourBand(
        name="phastCons (100 way vertebrate, 46 way mammalian, 30 way mammalian)",
        fields=("phastcons_100_way_vertebrate", "phastcons_46_way_mammalian", "phastcons_30_way_mammalian"),
        max_benign_threshold=None,
        pathogenic_threshold=0.85,
        source=_VG_DEFAULT,
        source_detail="Siepel A, Bejerano G, Pedersen JS, et al. Evolutionarily conserved elements in "
                      "vertebrate, insect, worm, and yeast genomes. Genome Res. 2005;15(8):1034-1050. "
                      "phastCons is a probability of lying in a conserved element; the 0.85 cutoff is a "
                      "VariantGrid display choice.",
        source_url="https://doi.org/10.1101/gr.3715005",
    ),
    ColourBand(
        name="phyloP 100 way vertebrate",
        fields=("phylop_100_way_vertebrate",),
        max_benign_threshold=None,
        pathogenic_threshold=1.4,
        source=_VG_DEFAULT,
        source_detail="Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. Detection of nonneutral substitution "
                      "rates on mammalian phylogenies. Genome Res. 2010;20(1):110-121. Cutoff is a "
                      "VariantGrid display choice.",
        source_url="https://doi.org/10.1101/gr.097857.109",
    ),
    ColourBand(
        name="phyloP 30 way mammalian",
        fields=("phylop_30_way_mammalian",),
        max_benign_threshold=None,
        pathogenic_threshold=1.0,
        source=_VG_DEFAULT,
        source_detail="Pollard 2010 (as above). Cutoff is a VariantGrid display choice.",
        source_url="https://doi.org/10.1101/gr.097857.109",
    ),
    ColourBand(
        name="phyloP 46 way mammalian",
        fields=("phylop_46_way_mammalian",),
        max_benign_threshold=None,
        pathogenic_threshold=None,
        source=_VG_DEFAULT,
        source_detail="Pollard 2010 (as above).",
        note="Left uncoloured - the 46 way alignment scores run on a different scale to the 100 way and "
             "30 way sets, and we have no cutoff for it.",
        source_url="https://doi.org/10.1101/gr.097857.109",
    ),
    # Grantham - chemical distance between the reference and alternate amino acid, not a prediction.
    ColourBand(
        name="Grantham distance",
        fields=("grantham",),
        max_benign_threshold=57.0,
        pathogenic_threshold=105.0,
        source=_VG_DEFAULT,
        source_detail="Grantham R. Amino acid difference formula to help explain protein evolution. "
                      "Science. 1974;185(4154):862-864. Distance 0-215 between the reference and alternate "
                      "residue; the band edges are a VariantGrid display choice.",
        source_url="https://doi.org/10.1126/science.185.4154.862",
        note="Chemical difference between amino acids, not a pathogenicity prediction.",
    ),
    # GERP++ RS - conserved end only, as for phastCons / phyloP.
    ColourBand(
        name="GERP++ RS",
        fields=("gerp_pp_rs",),
        max_benign_threshold=None,
        pathogenic_threshold=3.0,
        source=_VG_DEFAULT,
        source_detail="Davydov EV, Goode DL, Sirota M, et al. Identifying a high fraction of the human "
                      "genome to be under selective constraint using GERP++. PLoS Comput Biol. "
                      "2010;6(12):e1001025. Rejected-substitutions score, range -12 to 6.8; the cutoff is a "
                      "VariantGrid display choice.",
        source_url="https://doi.org/10.1371/journal.pcbi.1001025",
    ),
)


@dataclass(frozen=True)
class EffectNodeFilter:
    """ An Effect (Damage) node filter over an annotation that has no per-tool PP3/BP4 score: a categorical
        call (SpliceVault, UTRAnnotator) or a score with no published cutoff (AlphaGenome splicing). Listed on
        the Pathogenicity Thresholds page with the rule the node applies and where it comes from (#2561). """
    name: str
    fields: tuple[str, ...]
    rule: str
    source: str
    source_detail: str
    source_url: str
    note: str = ""
    # further (detail, url) citations behind the rule, for the page's references
    references: tuple[tuple[str, str], ...] = ()


# AlphaGenome merged splicing score: splice site + splice site usage (each 0-1) + splice junction / 5, so 0 to ~3.
# The paper proposes no pathogenicity cutoff, so the slider covers the whole scale with no colour band
ALPHAGENOME_SPLICING_MIN = 0.0
ALPHAGENOME_SPLICING_MAX = 3.0
ALPHAGENOME_SPLICING_STEP = 0.1

_DAWES_2023 = (
    "Dawes R, Bournazos AM, Bryen SJ, et al. SpliceVault predicts the precise nature of variant-associated "
    "mis-splicing. Nat Genet. 2023;55(2):324-332.",
    "https://doi.org/10.1038/s41588-022-01293-8",
)
_MARTIN_GEARY_2025 = (
    "Martin-Geary AC, Blakes AJM, Dawes R, et al. Systematic identification of disease-causing promoter and "
    "untranslated region variants in 8040 undiagnosed individuals with rare disease. Genome Med. 2025;17:40.",
    "https://doi.org/10.1186/s13073-025-01464-2",
)
_WHIFFIN_2020 = (
    "Whiffin N, Karczewski KJ, Zhang X, et al. Characterising the loss-of-function impact of 5' untranslated "
    "region variants in 15,708 individuals. Nat Commun. 2020;11(1):2523.",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC7253449/",
)
_ZHANG_2021 = (
    "Zhang X, Wakeling M, Ware J, Whiffin N. Annotating high-impact 5'untranslated region variants with the "
    "UTRannotator. Bioinformatics. 2021;37(8):1171-1173.",
    "https://doi.org/10.1093/bioinformatics/btaa783",
)
_ELLINGFORD_2022 = (
    "Ellingford JM, Ahn JW, Bagnall RD, et al. Recommendations for clinical interpretation of variants found in "
    "non-coding regions of the genome. Genome Med. 2022;14(1):73.",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC9295495/",
)

EFFECT_NODE_FILTERS: tuple[EffectNodeFilter, ...] = (
    EffectNodeFilter(
        name="AlphaGenome splicing",
        fields=("alphagenome_splicing",),
        rule=f"Score ≥ slider ({ALPHAGENOME_SPLICING_MIN:g} to {ALPHAGENOME_SPLICING_MAX:g}), an OR term in the "
             "Splicing filter beside SpliceAI and dbscSNV, sharing its Required / Allow Null.",
        source="Avsec 2025 (no calibration)",
        source_detail="Avsec Ž, Latysheva N, Cheng J, et al. AlphaGenome: advancing regulatory variant effect "
                      "prediction with a unified DNA sequence model. bioRxiv 2025.06.25.661532, "
                      "published in Nature (doi 10.1038/s41586-025-10014-0). Merged score = splice site + "
                      "splice site usage + splice junction / 5.",
        source_url="https://doi.org/10.1038/s41586-025-10014-0",
        note="No published pathogenicity cutoff - the paper's thresholds are quantiles against common SNPs. "
             "GRCh38 SNVs only, and licensed for non-commercial use only, so it is offered only "
             "where the deployment allows non-commercial data.",
    ),
    EffectNodeFilter(
        name="SpliceVault out of frame",
        fields=("splicevault_out_of_frame_events",),
        rule="All top events out of frame: Frameshift:N/N (N > 0) - the authors' \"consistent with null "
             "outcomes\" case for PVS1 at Very Strong. Any top event out of frame: Frameshift > 0.",
        source="Dawes 2023",
        source_detail=_DAWES_2023[0],
        source_url=_DAWES_2023[1],
        note="Predicts what kind of mis-splicing follows from a variant already predicted to disrupt a splice "
             "site (the plugin only reports sites where SpliceAI predicts a loss), not whether splicing is "
             "disrupted. No sample count minimum: the paper counts an event from a single read. Low PPV, and "
             "less reliable for genes with tissue-specific or developmental alternative splicing.",
    ),
    EffectNodeFilter(
        name="High-impact 5' UTR (UTRAnnotator)",
        fields=("utr_annotator_consequence", "utr_annotator_annotation"),
        rule="Any one of: uAUG gain creating an oORF with strong or moderate Kozak; uSTOP loss with no "
             "alternative stop before the CDS (an oORF) with strong or moderate Kozak; uAUG loss with strong "
             "Kozak; uFrameshift creating an oORF with strong or moderate Kozak. Two options: \"any frame\" is "
             "the Martin-Geary 2025 rule as published; \"out-of-frame oORF only\" also requires the oORF of a "
             "uAUG gain, uSTOP loss or uFrameshift to be out of frame with the CDS (uAUG loss is unchanged).",
        source="Martin-Geary 2025",
        source_detail=_MARTIN_GEARY_2025[0],
        source_url=_MARTIN_GEARY_2025[1],
        note="The filter Martin-Geary 2025 used to find new diagnoses in undiagnosed rare-disease cases. uAUG "
             "gains forming an oORF are under selection similar to missense variants (Whiffin 2020). An "
             "oORF in frame with the CDS makes an N-terminal extension rather than disrupting translation, so is "
             "not necessarily loss of function (Ellingford 2022) - the reason for the out-of-frame option.",
        references=(_WHIFFIN_2020, _ZHANG_2021, _ELLINGFORD_2022),
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
