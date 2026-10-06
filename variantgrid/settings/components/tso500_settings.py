""" Illumina TSO 500 (DRAGEN) import and calls - the lab policy the upload/tso500 parsers and
    seqauto.models.models_seqauto.DragenTSO500CombinedVariantOutput read. Imported by default_settings,
    so every deployment has them; an env file overrides them by name """

# A TSO500 CombinedVariantOutput's 'Pair ID' is either the patient's code on its own or the pair's
# whole sample name with the code as a field inside it - SA Path's pipeline writes both forms, and
# which one a run carries is not something to rely on. This regex's 'patient_code' group reads either
# ('5_FMC080_FCUP_2619115319' -> 'FMC080', 'FMC080' -> 'FMC080'), the leading sequencing sample ID
# changing on re-sequencing while the code stays. A lab naming pairs some other way sets its own;
# None takes the whole Pair ID as the code, and a Pair ID the regex doesn't match accessions no patient
TSO500_PAIR_ID_PATIENT_CODE_REGEX = r"^(?:\d+_)?(?P<patient_code>[^_]+)(?:_|$)"
# Turning a pair's MSI and TMB numbers into a call is lab policy rather than vendor output - DRAGEN
# writes the numbers and no call. None leaves the call blank, so a deployment that has not set its
# policy reports the measure as not able to be determined, as it does today. A band list is
# [(lower bound, call), ...], the call being the first whose lower bound the number reaches - the
# words are the lab's too (SA Path: MSI-High / MSI-Low / MSS), and go into the Omico JSON as written
TSO500_MSI_MIN_USABLE_SITES = None   # fewer 'Usable MSI Sites' than this and MSI cannot be called
TSO500_MSI_CALL_BANDS = None         # over 'Percent Unstable MSI Sites', eg [(30, "MSI-High"), (10, "MSI-Low"), (0, "MSS")]
TSO500_TMB_CALL_BANDS = None         # over 'Total TMB' in mut/Mb, eg [(10, "High"), (0, "Low")]
TSO500_GIS_CALL_BANDS = None                # over 'Genomic Instability Score', eg [(42, "POSITIVE"), (0, "NEGATIVE")]
TSO500_GIS_MIN_TUMOR_FRACTION = None        # below this caller tumour fraction a GIS under the top band gets no call, eg 0.23
TSO500_MSI_HIGH_MIN_TUMOR_FRACTION = None   # below this caller tumour fraction the top MSI band gets no call, eg 0.20
# AllFusions rows DRAGEN did not keep are rescued when Score is strictly above the min and Filter is exactly
# the given string. None for either means no rescue
TSO500_FUSION_RESCUE_MIN_SCORE = None       # eg 0.5
TSO500_FUSION_RESCUE_FILTER = None          # eg "FAIL;LOW_MAPQ"
# DRAGEN's MetricsOutput carries its own LSL/USL guideline per metric, and that is the policy a
# library QC category is judged by. Where a lab's own methods paragraph quotes a different number
# this overrides it: {(section name, metric): (lsl, usl)}, eg
# {("RNA Library QC Metrics", "TOTAL_ON_TARGET_READS"): (9000000, None)} - keyed on the section too,
# since MEDIAN_INSERT_SIZE appears in two of them. A row records which of the two it was judged by
TSO500_LIBRARY_QC_GUIDELINES = None
