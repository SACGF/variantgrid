import io
import json
from enum import Enum
from typing import Optional

import pandas as pd
from django.core.management import BaseCommand
from pandas import Series

from classification.enums import EvidenceCategory, EvidenceKeyValueType
from classification.models import EvidenceKey

data = \
'''key,label,description,examples,options,see,evidence_category,value_type,mandatory,order,allow_custom_values,sub_label,hide,variantgrid_column_id,immutable,created,modified,copy_consensus,default_crit_evaluation,max_share_level,crit_allows_override_strengths,crit_uses_points,namespace_overrides
amp:level_a,Level A,"Biomarkers that predict response or resistance to US FDA-approved therapies for a specific type of tumor or have been included in professional guidelines as therapeutic, diagnostic, and/or prognostic biomarkers for specific types of tumors.

See <a href=""https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/"">AMP/ASCO/CAP 2017</a>",,"[{""key"": ""therapeutic"", ""index"": 1}, {""key"": ""diagnostic"", ""index"": 2}, {""key"": ""prognostic"", ""index"": 3}, {""key"": ""unspecified"", ""index"": 4}]",,SC,M,FALSE,2,FALSE,Evidence of strong clinical significance,FALSE,,FALSE,2024-02-28 15:10:42.327447+11,2024-02-28 15:10:42.327473+11,TRUE,,logged_in_users,FALSE,FALSE,
amp:level_b,Level B,"Biomarkers that predict response or resistance to a therapy, or have diagnostic and/or prognostic significance, for a specific type of tumor based on well-powered studies with consensus from experts in the field.

See <a href=""https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/"">AMP/ASCO/CAP 2017</a>",,"[{""key"": ""therapeutic"", ""index"": 1}, {""key"": ""diagnostic"", ""index"": 2}, {""key"": ""prognostic"", ""index"": 3}, {""key"": ""unspecified"", ""index"": 4}]",,SC,M,FALSE,2,FALSE,Evidence of strong clinical significance,FALSE,,FALSE,2024-02-28 15:10:42.331352+11,2024-02-28 15:10:42.331363+11,TRUE,,logged_in_users,FALSE,FALSE,
amp:level_c,Level C,"Biomarkers that predict response or resistance to therapies approved by FDA or professional societies for a different tumor type (i.e. off-label use of a drug), serve as inclusion criteria for clinical trials, or have diagnostic and/or prognostic significance based on the results of multiple small studies.

See <a href=""https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/"">AMP/ASCO/CAP 2017</a>",,"[{""key"": ""therapeutic"", ""index"": 1}, {""key"": ""diagnostic"", ""index"": 2}, {""key"": ""prognostic"", ""index"": 3}, {""key"": ""unspecified"", ""index"": 4}]",,SC,M,FALSE,2,FALSE,Evidence of potential clinical significance,FALSE,,FALSE,2024-02-28 15:10:42.332089+11,2024-02-28 15:10:42.332096+11,TRUE,,logged_in_users,FALSE,FALSE,
amp:level_d,Level D,"Biomarkers that show plausible therapeutic significance based on preclinical studies, or may assist disease diagnosis and/or prognosis themselves or along with other biomarkers based on small studies or multiple case reports with no consensus.

See <a href=""https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/"">AMP/ASCO/CAP 2017</a>",,"[{""key"": ""therapeutic"", ""index"": 1}, {""key"": ""diagnostic"", ""index"": 2}, {""key"": ""prognostic"", ""index"": 3}, {""key"": ""unspecified"", ""index"": 4}]",,SC,M,FALSE,2,FALSE,Evidence of potential clinical significance,FALSE,,FALSE,2024-02-28 15:10:42.332585+11,2024-02-28 15:10:42.332591+11,TRUE,,logged_in_users,FALSE,FALSE,
assertion_method,Assertion method,The curation method used to interpret the variant. Novel or in-house curation methods can be added as custom values.,[],"[{""key"": ""acmg"", ""index"": 1, ""label"": ""ACMG (PMID:25741868)"", ""aliases"": [""ACMG Guidelines""], ""namespaces"": [""acmg""]}, {""key"": ""bayesian_acmg"", ""index"": 2, ""label"": ""ACMG Bayesian (PMID:29300386)"", ""namespaces"": [""acmg""]}, {""key"": ""modified_ACMG"", ""index"": 3, ""label"": ""Modified ACMG"", ""namespaces"": [""acmg""]}, {""key"": ""horak"", ""index"": 5, ""label"": ""ClinGen/CGC/VICC 2022"", ""namespaces"": [""horak""]}, {""key"": ""amp"", ""label"": ""AMP/ASCO/CAP 2017"", ""namespaces"": [""amp""]}]",,HI,M,FALSE,5,TRUE,,FALSE,,FALSE,2019-10-23 09:46:59.64957+11,2024-02-22 12:26:17.473026+11,FALSE,,logged_in_users,FALSE,FALSE,
horak:om1,OM1,"Located in a critical and well-established part of a functional domain (e.g., active site of an enzyme).",,,,CP,C,FALSE,10000,FALSE,Located in a critical and well-established part of a functional domain,FALSE,,FALSE,2023-11-14 12:48:45.23003+11,2024-02-07 10:36:10.746366+11,TRUE,PM,logged_in_users,TRUE,TRUE,
horak:om2,OM2,"Protein length changes as a result of in-frame deletions\/insertions in a known oncogene, or tumor suppressor gene or stop-loss variants in a known tumor suppressor gene.",,,,CP,C,FALSE,10000,FALSE,Protein length changes as a result of in-frame deletions\/insertions ,FALSE,,FALSE,2023-11-14 12:48:45.230503+11,2024-02-07 10:36:10.748312+11,TRUE,PM,logged_in_users,TRUE,TRUE,
horak:om3,OM3,Missense variant at an amino acid residue where a different missense variant determined to be oncogenic (using this standard) has been documented. Amino acid difference from reference amino acid should be greater or at least approximately the same as for missense change determined to be oncogenic.,,,,CP,C,FALSE,10000,FALSE,Missense change at an amino acid residue where a different missense change determined to be oncogenic has been documented,FALSE,,FALSE,2023-11-14 12:48:45.230962+11,2024-02-07 10:36:10.749057+11,TRUE,PM,logged_in_users,TRUE,TRUE,
horak:om4,OM4,"Located in one of the hotspots in cancerhotspots.org with <50 samples with a somatic variant at the same amino acid position, and the same amino acid change count in cancerhotspots.org is at least 10.",,,,F,C,FALSE,10000,FALSE,Cancer hotspot with moderate frequency of recurrence,FALSE,,FALSE,2023-11-14 12:48:45.231444+11,2024-02-07 10:36:10.749738+11,TRUE,PM,logged_in_users,TRUE,TRUE,
horak:op1,OP1,"All used lines of computational evidence support an oncogenic effect of a variant (conservation\/evolutionary, splicing impact, etc.).",,,,CP,C,FALSE,10000,FALSE,All utilized lines of computational evidence support oncogenicity,FALSE,,FALSE,2023-11-14 12:48:45.232863+11,2024-02-07 10:36:10.750663+11,TRUE,PP,logged_in_users,TRUE,TRUE,
horak:op2,OP2,Somatic variant in a gene in a malignancy with a single genetic etiology. Example: retinoblastoma is caused by bi-allelic RB1 inactivation.,,,,CP,C,FALSE,10000,FALSE,Somatic variant in a gene in a malignancy with a single genetic etiology,FALSE,,FALSE,2023-11-14 12:48:45.233391+11,2024-02-07 10:36:10.739329+11,TRUE,PP,logged_in_users,TRUE,TRUE,
horak:op3,OP3,Located in one of the hotspots in cancerhotspots.org and the particular amino acid change count in cancerhotspots.org is below 10.,,,,F,C,FALSE,10000,FALSE,Cancer hotspots with low frequency of recurrence,FALSE,,FALSE,2023-11-14 12:48:45.233991+11,2024-02-07 10:36:10.740346+11,TRUE,PP,logged_in_users,TRUE,TRUE,
horak:op4,OP4,Absent from controls (or at an extremely low frequency) in Genome Aggregation Database (gnomAD).,,,,P,C,FALSE,10000,FALSE,Absent in population databases,FALSE,,FALSE,2023-11-14 12:48:45.234487+11,2024-02-07 10:36:10.741444+11,TRUE,PP,logged_in_users,TRUE,TRUE,
horak:os1,OS1,Same amino acid change as a previously established oncogenic variant (using this standard) regardless of nucleotide change. Example: Val→Leu caused by either G>C or G>T in the same codon.,,,,CP,C,FALSE,10000,FALSE,Same amino acid change as a previously established oncogenic variant,FALSE,,FALSE,2023-11-14 12:48:45.228323+11,2024-02-07 10:36:10.736163+11,TRUE,PS,logged_in_users,TRUE,TRUE,
horak:os2,OS2,"Well-established in vitro or in vivo functional studies, supportive of an oncogenic effect of the variant.
",,,,F,C,FALSE,10000,FALSE,Well-established functional studies supportive of an oncogenic effect,FALSE,,FALSE,2023-11-14 12:48:45.22903+11,2024-02-07 10:36:10.737464+11,TRUE,PS,logged_in_users,TRUE,TRUE,
horak:os3,OS3,"Located in one of the hotspots in cancerhotspots.org with at least 50 samples with a somatic variant at the same amino acid position, and the same amino acid change count in cancerhotspots.org in at least 10 samples.",,,,F,C,FALSE,10000,FALSE,Cancer hotspot with high frequency of recurrence,FALSE,,FALSE,2023-11-14 12:48:45.22952+11,2024-02-07 10:36:10.745501+11,TRUE,PS,logged_in_users,TRUE,TRUE,
horak:ovs1,OVS1,"Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multi-exon deletion) in a bona fide tumor suppressor gene.",,,,CP,C,FALSE,10000,FALSE,Null variant in tumor suppressor,FALSE,,FALSE,2023-11-14 12:48:45.226365+11,2024-02-07 10:36:10.734724+11,TRUE,PVS,logged_in_users,TRUE,TRUE,
horak:sbp1,SBP1,"All used lines of computational evidence suggest no effect of a variant (conservation\/evolutionary, splicing impact, etc.).",,,,CP,C,FALSE,10000,FALSE,All utilized lines of computational evidence suggest no impact of a variant,FALSE,,FALSE,2023-11-14 12:48:45.236328+11,2024-02-07 10:36:10.733039+11,TRUE,BS,logged_in_users,TRUE,TRUE,
horak:sbp2,SBP2,A synonymous (silent) variant for which splicing prediction algorithms predict no effect on the splice consensus sequence nor the creation of a new splice site and the nucleotide is not highly conserved.,,,,CP,C,FALSE,10000,FALSE,Silent mutation (no predicted impact on splicing),FALSE,,FALSE,2023-11-14 12:48:45.236782+11,2024-02-07 10:36:10.733919+11,TRUE,BS,logged_in_users,TRUE,TRUE,
horak:sbs1,SBS1,"Minor allele frequency is >1% in Genome Aggregation Database (gnomAD) in any of 5 general continental populations: African, East Asian, European (Non-Finnish), Latino, and South Asian.",,,,P,C,FALSE,10000,FALSE,MAF >1%,FALSE,,FALSE,2023-11-14 12:48:45.235415+11,2024-02-07 10:36:10.743367+11,TRUE,BS,logged_in_users,TRUE,TRUE,
horak:sbs2,SBS2,Well-established in vitro or in vivo functional studies show no oncogenic effects.,,,,F,C,FALSE,10000,FALSE,Well-established functional studies show no oncogenic effects,FALSE,,FALSE,2023-11-14 12:48:45.235869+11,2024-02-07 10:36:10.744255+11,TRUE,BS,logged_in_users,TRUE,TRUE,
horak:sbvs1,SBVS1,"Minor allele frequency is >5% in Genome Aggregation Database (gnomAD) in any of 5 general continental populations: African, East Asian, European (Non-Finnish), Latino, and South Asian.",,,,P,C,FALSE,10000,FALSE,MAF >5%,FALSE,,FALSE,2023-11-14 12:48:45.234959+11,2024-02-07 10:36:10.742394+11,TRUE,BA,logged_in_users,TRUE,TRUE,
sample_type,Sample Type,Sample type used for sequencing.,[],"[{""key"": ""amniocentesis"", ""index"": 3, ""label"": """"}, {""key"": ""blood"", ""index"": 4, ""label"": """"}, {""key"": ""blood_spot"", ""index"": 8, ""label"": """"}, {""key"": ""buccal"", ""index"": 5, ""label"": """"}, {""key"": ""cf_dna"", ""index"": 10, ""label"": ""cfDNA""}, {""key"": ""cvs"", ""index"": 2, ""label"": ""CVS""}, {""key"": ""hair"", ""index"": 7, ""label"": """"}, {""key"": ""saliva"", ""index"": 6, ""label"": """"}, {""key"": ""skin_biopsy"", ""index"": 9}, {""key"": ""tumor"", ""index"": 1, ""label"": """"}]",,HT,S,FALSE,0,TRUE,,FALSE,,FALSE,2019-10-23 09:46:59.662298+11,2024-02-22 12:26:17.476754+11,FALSE,,logged_in_users,FALSE,FALSE,
somatic:clinical_significance,Somatic Clinical Significance,"
Evidence-based variant categorisation. Somatic variants are classified into tiers based on their level of clinical significance in cancer diagnosis, prognosis, and/or therapeutics.<br/>
Tier I: Variants of Strong Clinical Significance<br/>
Tier II: Variants of Potential Clinical Significance<br/>
Tier III: Variants of Unknown Clinical Significance<br/>
Tier IV: Benign or Likely Benign Variants<br/>
Tier I/II: If specific tier of clinical significance is not denoted further than I/II as a group
",,"[{""key"": ""tier_1"", ""tier"": ""1"", ""index"": 1, ""label"": ""Tier I"", ""aliases"": [""1""]}, {""key"": ""tier_2"", ""tier"": ""2"", ""index"": 2, ""label"": ""Tier II"", ""aliases"": [""2""]}, {""key"": ""tier_3"", ""tier"": ""3"", ""index"": 3, ""label"": ""Tier III"", ""aliases"": [""1""]}, {""key"": ""tier_4"", ""tier"": ""4"", ""index"": 4, ""label"": ""Tier IV"", ""aliases"": [""4""]}, {""key"": ""tier_1_or_2"", ""tier"": ""1"", ""index"": 5, ""label"": ""Tier I / II"", ""aliases"": [""1 or 2""]}]",,HI,S,FALSE,1,FALSE,,FALSE,,FALSE,2023-11-16 11:20:11.986764+11,2024-02-22 12:26:17.475296+11,TRUE,,logged_in_users,FALSE,FALSE,
somatic:tumour_cellularity,Tumour Cellularity,,,,,HT,N,FALSE,0,FALSE,,FALSE,,FALSE,2024-02-22 11:59:49.852131+11,2024-02-22 11:59:49.852146+11,FALSE,,logged_in_users,FALSE,FALSE,'''


lab_config = {
    "namespaces": [
        "somatic",
        "amp",
        "horak"
    ],
    "clinical_significance": {
        "label": "Classification",
        "sub_label": "formerly Clinical Significance",
        "options": [
            {
                "vg": "1",
                "key": "B",
                "index": 1,
                "label": "Benign",
                "bucket": 1,
                "clinvar": "Benign"
            },
            {
                "vg": "2",
                "key": "LB",
                "index": 2,
                "label": "Likely Benign",
                "bucket": 1,
                "clinvar": "Likely benign"
            },
            {
                "vg": "3",
                "key": "VUS",
                "index": 3,
                "label": "VUS",
                "bucket": 2,
                "aliases": [
                    "US",
                    "VOUS"
                ],
                "clinvar": "Uncertain significance"
            },
            {
                "vg": "3",
                "key": "VUS_A",
                "index": 4,
                "label": "VUS A",
                "bucket": 2,
                "aliases": [
                    "US_A"
                ],
                "clinvar": "Uncertain significance"
            },
            {
                "vg": "3",
                "key": "VUS_B",
                "index": 5,
                "label": "VUS B",
                "bucket": 2,
                "aliases": [
                    "US_B"
                ],
                "clinvar": "Uncertain significance"
            },
            {
                "vg": "3",
                "key": "VUS_C",
                "index": 6,
                "label": "VUS C",
                "bucket": 2,
                "aliases": [
                    "US_C"
                ],
                "clinvar": "Uncertain significance"
            },
            {
                "vg": "4",
                "key": "LP",
                "index": 7,
                "label": "Likely Pathogenic",
                "bucket": 3,
                "clinvar": "Likely pathogenic",
                "namespace": "acmg"
            },
            {
                "vg": "5",
                "key": "P",
                "index": 8,
                "label": "Pathogenic",
                "bucket": 3,
                "clinvar": "Pathogenic",
                "namespace": "acmg"
            },
            {
                "vg": "4",
                "key": "LO",
                "index": 7,
                "label": "Likely Oncogenic",
                "bucket": 3,
                "clinvar": "Likely oncogenic",
                "namespace": "horak"
            },
            {
                "vg": "5",
                "key": "O",
                "index": 8,
                "label": "Oncogenic",
                "bucket": 3,
                "clinvar": "Oncogenic",
                "namespace": "horak"
            },
            {
                "vg": "0",
                "key": "D",
                "index": 9,
                "label": "Drug Response",
                "bucket": None,
                "clinvar": "drug response",
                "namespace": "germline"
            },
            {
                "vg": "0",
                "key": "R",
                "index": 10,
                "label": "Risk Factor",
                "bucket": 4,
                "clinvar": "risk factor",
                "namespace": "germline"
            }
        ]
    }
}

class CE(str, Enum):
    BENIGN_STANDALONE = 'BA'
    BENIGN_STRONG = 'BS'
    BENIGN_MODERATE = 'BM'  # Not a standard ACMG Strength
    BENIGN_SUPPORTING = 'BP'
    PATHOGENIC_SUPPORTING = 'PP'
    PATHOGENIC_MODERATE = 'PM'
    PATHOGENIC_STRONG = 'PS'
    PATHOGENIC_VERY_STRONG = 'PVS'


BENIGN = [
    CE.BENIGN_STANDALONE,
    CE.BENIGN_STRONG,
    CE.BENIGN_MODERATE,
    CE.BENIGN_SUPPORTING
]

PATHOGENIC = [
    CE.PATHOGENIC_SUPPORTING,
    CE.PATHOGENIC_MODERATE,
    CE.PATHOGENIC_STRONG,
    CE.PATHOGENIC_VERY_STRONG
]

POINTS = {
    CE.BENIGN_STANDALONE: -8,
    CE.BENIGN_STRONG: -4,
    CE.BENIGN_MODERATE: -2,
    CE.BENIGN_SUPPORTING: -1,
    CE.PATHOGENIC_SUPPORTING: 1,
    CE.PATHOGENIC_MODERATE: 2,
    CE.PATHOGENIC_STRONG: 4,
    CE.PATHOGENIC_VERY_STRONG: 8
}


def ce_as_option(ce: CE, index: int, is_default: bool = False):
    key = ce.value
    points = POINTS[key]
    plural = "" if points in {1,-1} else "s"
    label: str
    if is_default:
        label = f"⭐ Met: {points} Point{plural}"
    else:
        label = f"❗Override: {points} Point{plural}"
    return {
        "key": key,
        "label": label,
        "index": index
    }


class Command(BaseCommand):
    """
    Only use this to backport somatic curation functionality to VG3 (or to undo the backport)
    Not to be used with VG4 onwards
    """

    def add_arguments(self, parser):
        parser.add_argument('-direction', type=str, required=True, help='backport or forwardport (backport to add somatic functionality to VG3, forwardport once on VG4)')


    def process_row(self, row: Series, backporting: bool):
        key = row["key"]
        e_key, exists = EvidenceKey.objects.get_or_create(pk=key, defaults={"evidence_category": EvidenceCategory.INTERPRETATION})
        if exists:
            print(f"Insertin {key}")
        else:
            print(f"Updating {key}")
        e_key.label = row["label"]
        e_key.description = row["description"]
        e_key.examples = row["examples"]
        e_key.copy_consensus = row["copy_consensus"]
        e_key.see = row["see"]
        e_key.order = row["order"]
        e_key.allow_custom_values = row["allow_custom_values"]

        value_type = row["value_type"]
        default_crit_evaluation = row["default_crit_evaluation"]

        options: Optional[list] = None
        if row_options := row["options"]:
            options = json.loads(row_options)

        # BACKPORT
        if backporting:
            if key.startswith("horak:"):
                value_type = EvidenceKeyValueType.SELECT
                default_ce = CE(default_crit_evaluation)
                reference_values = list(BENIGN if default_ce in BENIGN else PATHOGENIC)

                options = [{
                    "key": "NM",
                    "label": "⚪ Not Met",
                    "index": 1
                }, {
                    "key": "NA",
                    "label": "⚪ Not Applicable",
                    "index": 2
                }, ce_as_option(default_ce, index=3, is_default=True)]
                use_index = 4
                for reference_value in reference_values:
                    if reference_value != default_ce:
                        options.append(ce_as_option(reference_value, index=use_index, is_default=False))
                        use_index += 1

            if key == "assertion_method":
                value_type = EvidenceKeyValueType.SELECT
        # END BACKPORT
        else:
            # Forward port, provide evidence keys with newer properties
            e_key.default_crit_evaluation = row["default_crit_evaluation"]
            e_key.crit_allows_override_strengths = row["crit_allows_override_strengths"]
            e_key.crit_uses_points = row["crit_uses_points"]
            e_key.namespace_overrides = row["crit_uses_points"]

        e_key.value_type = value_type
        e_key.options = options
        e_key.save()

    def handle(self, *args, **options):
        buffer = io.StringIO(data)
        df = pd.read_csv(buffer, sep=",", header=0).fillna("").convert_dtypes()

        backport: bool
        if options["direction"] == "backport":
            backport = True
        elif options["direction"] == "forwardport":
            backport = False
        else:
            raise ValueError("-direction must be backport of forwardport")

        for idx, row in df.iterrows():
            self.process_row(row, backporting=backport)
        print("done")

        if backport:
            print("Lab Config should be updated for the specific labs to be the following")
        else:
            print("Should be safe to remove Lab Config from a lab specifically setup")
        print(json.dumps(lab_config))