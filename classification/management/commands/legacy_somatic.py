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
'''key,mandatory,order,label,description,examples,options,see,evidence_category,value_type,allow_custom_values,sub_label,hide,variantgrid_column_id,immutable,created,modified,copy_consensus,default_crit_evaluation,max_share_level,crit_allows_override_strengths,crit_uses_points,namespace_overrides
amp:level_a,FALSE,2,Level A,"Biomarkers that predict response or resistance to US FDA-approved therapies for a specific type of tumour or have been included in professional guidelines as therapeutic, diagnostic, and/or prognostic biomarkers for specific types of tumours.

See <a href=""https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/"">AMP/ASCO/CAP 2017</a>",,"[{""key"": ""therapeutic"", ""index"": 1}, {""key"": ""diagnostic"", ""index"": 2}, {""key"": ""prognostic"", ""index"": 3}, {""key"": ""unspecified"", ""index"": 4}]",,SC,M,FALSE,Evidence of strong clinical significance,FALSE,,FALSE,2024-02-28 04:36:54.090957+00,2024-02-28 04:36:54.090982+00,TRUE,,logged_in_users,FALSE,FALSE,
amp:level_b,FALSE,2,Level B,"Biomarkers that predict response or resistance to a therapy, or have diagnostic and/or prognostic significance, for a specific type of tumour based on well-powered studies with consensus from experts in the field.

See <a href=""https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/"">AMP/ASCO/CAP 2017</a>",,"[{""key"": ""therapeutic"", ""index"": 1}, {""key"": ""diagnostic"", ""index"": 2}, {""key"": ""prognostic"", ""index"": 3}, {""key"": ""unspecified"", ""index"": 4}]",,SC,M,FALSE,Evidence of strong clinical significance,FALSE,,FALSE,2024-02-28 04:36:54.094618+00,2024-02-28 04:36:54.094635+00,TRUE,,logged_in_users,FALSE,FALSE,
amp:level_c,FALSE,2,Level C,"Biomarkers that predict response or resistance to therapies approved by FDA or professional societies for a different tumour type (i.e. off-label use of a drug), serve as inclusion criteria for clinical trials, or have diagnostic and/or prognostic significance based on the results of multiple small studies.

See <a href=""https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/"">AMP/ASCO/CAP 2017</a>",,"[{""key"": ""therapeutic"", ""index"": 1}, {""key"": ""diagnostic"", ""index"": 2}, {""key"": ""prognostic"", ""index"": 3}, {""key"": ""unspecified"", ""index"": 4}]",,SC,M,FALSE,Evidence of potential clinical significance,FALSE,,FALSE,2024-02-28 04:36:54.095854+00,2024-02-28 04:36:54.095868+00,TRUE,,logged_in_users,FALSE,FALSE,
amp:level_d,FALSE,2,Level D,"Biomarkers that show plausible therapeutic significance based on preclinical studies, or may assist disease diagnosis and/or prognosis themselves or along with other biomarkers based on small studies or multiple case reports with no consensus.

See <a href=""https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/"">AMP/ASCO/CAP 2017</a>",,"[{""key"": ""therapeutic"", ""index"": 1}, {""key"": ""diagnostic"", ""index"": 2}, {""key"": ""prognostic"", ""index"": 3}, {""key"": ""unspecified"", ""index"": 4}]",,SC,M,FALSE,Evidence of potential clinical significance,FALSE,,FALSE,2024-02-28 04:36:54.097161+00,2024-02-28 04:36:54.097178+00,TRUE,,logged_in_users,FALSE,FALSE,
assertion_method,FALSE,5,Assertion method,The curation method used to interpret the variant. Novel or in-house curation methods can be added as custom values.,[],"[{""key"": ""acmg"", ""index"": 1, ""label"": ""ACMG Guidelines, 2015"", ""aliases"": [""ACMG Guidelines""], ""namespaces"": [""acmg""]}, {""key"": ""bayesian_acmg"", ""index"": 2, ""label"": ""ACMG Bayesian (PMID:29300386)"", ""namespaces"": [""acmg""]}, {""key"": ""modified_ACMG"", ""index"": 3, ""label"": ""Modified ACMG"", ""namespaces"": [""acmg""]}, {""key"": ""sherloc"", ""index"": 4, ""label"": ""Sherloc"", ""namespaces"": [""sherloc""]}, {""key"": ""horak"", ""index"": 5, ""label"": ""ClinGen/CGC/VICC 2022"", ""namespaces"": [""horak""]}, {""key"": ""amp"", ""label"": ""AMP/ASCO/CAP 2017"", ""namespaces"": [""amp""]}]",,HI,M,TRUE,,FALSE,,FALSE,2019-10-21 03:25:22.375997+00,2024-02-22 01:55:26.26714+00,FALSE,,logged_in_users,FALSE,FALSE,
clinical_significance,FALSE,1,Classification,"The classification of the specified variant in the context of the condition analysed. Values are as specified by ACMG/AMP and ClinGen/CGC/VICC 2022 with an optional extension to VUS.<br/>
Options are
<ul class=""compact"">
<li>Benign</li>
<li>Likely Benign</li>
<li>VUS *</li>
<li>Likely Pathogenic/Oncogenic</li>
<li>Pathogenic/Oncogenic</li>
</ul>
* VUS can optionally be subdivided into the following
<ul class=""compact"">
<li>VUS-C : Variant of unknown significance with trending towards benign</li>
<li>VUS-B : Variant of unknown significance</li>
<li>VUS-A : Variant of unknown significance with trending towards pathogenic/oncogenic</li>
</ul>
Additional options may include
<ul class=""compact"">
<li>Artefact</li>
<li>Drug Response</li>
<li>Risk Factor</li>
</ul>

See <a href=""https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/"">ACMG/AMP</a> for widely used Germline guidelines or
See <a href=""https://pubmed.ncbi.nlm.nih.gov/35101336/"">ClinGen/CGC/VICC 2022</a> for widely used Somatic guidelines",[],"[{""vg"": ""1"", ""key"": ""B"", ""index"": 1, ""label"": ""Benign"", ""bucket"": 1, ""clinvar"": ""Benign""}, {""vg"": ""2"", ""key"": ""LB"", ""index"": 2, ""label"": ""Likely Benign"", ""bucket"": 1, ""clinvar"": ""Likely benign""}, {""vg"": ""3"", ""key"": ""VUS"", ""index"": 3, ""label"": ""VUS"", ""bucket"": 2, ""aliases"": [""US"", ""VOUS""], ""clinvar"": ""Uncertain significance""}, {""vg"": ""3"", ""key"": ""VUS_A"", ""index"": 4, ""label"": ""VUS A"", ""bucket"": 2, ""aliases"": [""US_A""], ""clinvar"": ""Uncertain significance""}, {""vg"": ""3"", ""key"": ""VUS_B"", ""index"": 5, ""label"": ""VUS B"", ""bucket"": 2, ""aliases"": [""US_B""], ""clinvar"": ""Uncertain significance""}, {""vg"": ""3"", ""key"": ""VUS_C"", ""index"": 6, ""label"": ""VUS C"", ""bucket"": 2, ""aliases"": [""US_C""], ""clinvar"": ""Uncertain significance""}, {""vg"": ""4"", ""key"": ""LP"", ""index"": 7, ""label"": ""Likely Pathogenic"", ""bucket"": 3, ""clinvar"": ""Likely pathogenic"", ""namespace"": ""acmg""}, {""vg"": ""5"", ""key"": ""P"", ""index"": 8, ""label"": ""Pathogenic"", ""bucket"": 3, ""clinvar"": ""Pathogenic"", ""namespace"": ""acmg""}, {""vg"": ""4"", ""key"": ""LO"", ""index"": 7, ""label"": ""Likely Oncogenic"", ""bucket"": 3, ""clinvar"": ""Likely oncogenic"", ""namespace"": ""horak""}, {""vg"": ""5"", ""key"": ""O"", ""index"": 8, ""label"": ""Oncogenic"", ""bucket"": 3, ""clinvar"": ""Oncogenic"", ""namespace"": ""horak""}, {""vg"": ""0"", ""key"": ""D"", ""index"": 9, ""label"": ""Drug Response"", ""bucket"": null, ""clinvar"": ""drug response"", ""namespace"": ""germline""}, {""vg"": ""0"", ""key"": ""R"", ""index"": 10, ""label"": ""Risk Factor"", ""bucket"": 4, ""clinvar"": ""risk factor"", ""namespace"": ""germline""}]",,HI,S,FALSE,formerly Clinical Significance,FALSE,,FALSE,2019-10-21 03:25:22.377726+00,2024-02-06 23:42:50.167267+00,FALSE,,public,FALSE,FALSE,
horak:om1,FALSE,10000,OM1,"Located in a critical and well-established part of a functional domain (e.g., active site of an enzyme).",,,,CP,C,FALSE,Located in a critical and well-established part of a functional domain,FALSE,,FALSE,2023-11-14 01:57:08.942507+00,2024-02-06 23:42:50.132968+00,TRUE,PM,logged_in_users,TRUE,TRUE,
horak:om2,FALSE,10000,OM2,"Protein length changes as a result of in-frame deletions/insertions in a known oncogene, or tumour suppressor gene or stop-loss variants in a known tumour suppressor gene.",,,,CP,C,FALSE,Protein length changes as a result of in-frame deletions/insertions ,FALSE,,FALSE,2023-11-14 01:57:08.943963+00,2024-02-06 23:42:50.134841+00,TRUE,PM,logged_in_users,TRUE,TRUE,
horak:om3,FALSE,10000,OM3,Missense variant at an amino acid residue where a different missense variant determined to be oncogenic (using this standard) has been documented. Amino acid difference from reference amino acid should be greater or at least approximately the same as for missense change determined to be oncogenic.,,,,CP,C,FALSE,Missense change at an amino acid residue where a different missense change determined to be oncogenic has been documented,FALSE,,FALSE,2023-11-14 01:57:08.945444+00,2024-02-06 23:42:50.136536+00,TRUE,PM,logged_in_users,TRUE,TRUE,
horak:om4,FALSE,10000,OM4,"Located in one of the hotspots in cancerhotspots.org with <50 samples with a somatic variant at the same amino acid position, and the same amino acid change count in cancerhotspots.org is at least 10.",,,,F,C,FALSE,Cancer hotspot with moderate frequency of recurrence,FALSE,,FALSE,2023-11-14 01:57:08.946874+00,2024-02-06 23:42:50.138282+00,TRUE,PM,logged_in_users,TRUE,TRUE,
horak:op1,FALSE,10000,OP1,"All used lines of computational evidence support an oncogenic effect of a variant (conservation/evolutionary, splicing impact, etc.).",,,,CP,C,FALSE,All utilised lines of computational evidence support oncogenicity,FALSE,,FALSE,2023-11-14 01:57:08.948358+00,2024-02-06 23:42:50.139898+00,TRUE,PP,logged_in_users,TRUE,TRUE,
horak:op2,FALSE,10000,OP2,Somatic variant in a gene in a malignancy with a single genetic etiology. Example: retinoblastoma is caused by bi-allelic RB1 inactivation.,,,,CP,C,FALSE,Somatic variant in a gene in a malignancy with a single genetic etiology,FALSE,,FALSE,2023-11-14 01:57:08.949784+00,2024-02-06 23:42:50.141481+00,TRUE,PP,logged_in_users,TRUE,TRUE,
horak:op3,FALSE,10000,OP3,Located in one of the hotspots in cancerhotspots.org and the particular amino acid change count in cancerhotspots.org is below 10.,,,,F,C,FALSE,Cancer hotspots with low frequency of recurrence,FALSE,,FALSE,2023-11-14 01:57:08.951251+00,2024-02-06 23:42:50.14303+00,TRUE,PP,logged_in_users,TRUE,TRUE,
horak:op4,FALSE,10000,OP4,Absent from controls (or at an extremely low frequency) in Genome Aggregation Database (gnomAD).,,,,P,C,FALSE,Absent in population databases,FALSE,,FALSE,2023-11-14 01:57:08.952645+00,2024-02-06 23:42:50.144638+00,TRUE,PP,logged_in_users,TRUE,TRUE,
horak:os1,FALSE,10000,OS1,Same amino acid change as a previously established oncogenic variant (using this standard) regardless of nucleotide change. Example: Val→Leu caused by either G>C or G>T in the same codon.,,,,CP,C,FALSE,Same amino acid change as a previously established oncogenic variant,FALSE,,FALSE,2023-11-14 01:57:08.93797+00,2024-02-06 23:42:50.151106+00,TRUE,PS,logged_in_users,TRUE,TRUE,
horak:os2,FALSE,10000,OS2,"Well-established in vitro or in vivo functional studies, supportive of an oncogenic effect of the variant.",,,,F,C,FALSE,Well-established functional studies supportive of an oncogenic effect,FALSE,,FALSE,2023-11-14 01:57:08.939477+00,2024-02-06 23:42:50.152715+00,TRUE,PS,logged_in_users,TRUE,TRUE,
horak:os3,FALSE,10000,OS3,"Located in one of the hotspots in cancerhotspots.org with at least 50 samples with a somatic variant at the same amino acid position, and the same amino acid change count in cancerhotspots.org in at least 10 samples.",,,,F,C,FALSE,Cancer hotspot with high frequency of recurrence,FALSE,,FALSE,2023-11-14 01:57:08.941037+00,2024-02-06 23:42:50.154287+00,TRUE,PS,logged_in_users,TRUE,TRUE,
horak:ovs1,FALSE,10000,OVS1,"Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multi-exon deletion) in a bona fide tumour suppressor gene.",,,,CP,C,FALSE,Null variant in tumour suppressor,FALSE,,FALSE,2023-11-14 01:57:08.934762+00,2024-02-06 23:42:50.155885+00,TRUE,PVS,logged_in_users,TRUE,TRUE,
horak:sbp1,FALSE,10000,SBP1,"All used lines of computational evidence suggest no effect of a variant (conservation/evolutionary, splicing impact, etc.).",,,,CP,C,FALSE,All utilised lines of computational evidence suggest no impact of a variant,FALSE,,FALSE,2023-11-14 01:57:08.959269+00,2024-02-06 23:42:50.157435+00,TRUE,BP,logged_in_users,TRUE,TRUE,
horak:sbp2,FALSE,10000,SBP2,A synonymous (silent) variant for which splicing prediction algorithms predict no effect on the splice consensus sequence nor the creation of a new splice site and the nucleotide is not highly conserved.,,,,CP,C,FALSE,Silent mutation (no predicted impact on splicing),FALSE,,FALSE,2023-11-14 01:57:08.960767+00,2024-02-06 23:42:50.15896+00,TRUE,BP,logged_in_users,TRUE,TRUE,
horak:sbs1,FALSE,10000,SBS1,"Minor allele frequency is >1% in Genome Aggregation Database (gnomAD) in any of 5 general continental populations: African, East Asian, European (Non-Finnish), Latino, and South Asian.",,,,P,C,FALSE,MAF >1%,FALSE,,FALSE,2023-11-14 01:57:08.955589+00,2024-02-06 23:42:50.160557+00,TRUE,BS,logged_in_users,TRUE,TRUE,
horak:sbs2,FALSE,10000,SBS2,Well-established in vitro or in vivo functional studies show no oncogenic effects.,,,,F,C,FALSE,Well-established functional studies show no oncogenic effects,FALSE,,FALSE,2023-11-14 01:57:08.957108+00,2024-02-06 23:42:50.16213+00,TRUE,BS,logged_in_users,TRUE,TRUE,
horak:sbvs1,FALSE,10000,SBVS1,"Minor allele frequency is >5% in Genome Aggregation Database (gnomAD) in any of 5 general continental populations: African, East Asian, European (Non-Finnish), Latino, and South Asian.",,,,P,C,FALSE,MAF >5%,FALSE,,FALSE,2023-11-14 01:57:08.954058+00,2024-02-06 23:42:50.163765+00,TRUE,BA,logged_in_users,TRUE,TRUE,
sample_type,FALSE,0,Sample Type,Sample type used for sequencing.,[],"[{""key"": ""amniocentesis"", ""index"": 3, ""label"": """"}, {""key"": ""blood"", ""index"": 4, ""label"": """"}, {""key"": ""blood_spot"", ""index"": 8, ""label"": """"}, {""key"": ""bone marrow"", ""index"": 11}, {""key"": ""buccal"", ""index"": 5, ""label"": """"}, {""key"": ""cf_dna"", ""index"": 10, ""label"": ""cfDNA""}, {""key"": ""cvs"", ""index"": 2, ""label"": ""CVS""}, {""key"": ""hair"", ""index"": 7, ""label"": """"}, {""key"": ""saliva"", ""index"": 6, ""label"": """"}, {""key"": ""skin_biopsy"", ""index"": 9}, {""key"": ""tumor"", ""index"": 1, ""label"": ""Tumour""}]",,HT,S,TRUE,,FALSE,,FALSE,2019-10-21 03:25:22.390169+00,2024-02-23 05:09:28.869313+00,FALSE,,logged_in_users,FALSE,FALSE,
somatic:clinical_significance,FALSE,1,Somatic Clinical Significance,"
Evidence-based variant categorisation. Somatic variants are classified into tiers based on their level of clinical significance in cancer diagnosis, prognosis, and/or therapeutics.<br/>
Tier I: Variants of Strong Clinical Significance<br/>
Tier II: Variants of Potential Clinical Significance<br/>
Tier III: Variants of Unknown Clinical Significance<br/>
Tier IV: Benign or Likely Benign Variants<br/>
Tier I/II: If specific tier of clinical significance is not denoted further than I/II as a group
",,"[{""key"": ""tier_1"", ""tier"": ""1"", ""index"": 1, ""label"": ""Tier I"", ""aliases"": [""1""]}, {""key"": ""tier_2"", ""tier"": ""2"", ""index"": 2, ""label"": ""Tier II"", ""aliases"": [""2""]}, {""key"": ""tier_3"", ""tier"": ""3"", ""index"": 3, ""label"": ""Tier III"", ""aliases"": [""1""]}, {""key"": ""tier_4"", ""tier"": ""4"", ""index"": 4, ""label"": ""Tier IV"", ""aliases"": [""4""]}, {""key"": ""tier_1_or_2"", ""tier"": ""1"", ""index"": 5, ""label"": ""Tier I or II"", ""aliases"": [""1 or 2""]}]",,HI,S,FALSE,,FALSE,,FALSE,2023-11-17 01:51:50.969851+00,2024-02-22 01:55:26.270375+00,FALSE,,logged_in_users,FALSE,FALSE,
somatic:tumor_cellularity,FALSE,0,Tumour Cellularity,,,,,HT,N,FALSE,,FALSE,,FALSE,2024-02-22 01:44:48.501533+00,2024-02-22 01:44:48.50155+00,FALSE,,logged_in_users,FALSE,FALSE,'''

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
        evidence_category = row["evidence_category"]

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
        e_key.sub_label = row["sub_label"]

        value_type = row["value_type"]
        default_crit_evaluation = row["default_crit_evaluation"]

        options: Optional[list] = None
        if row_options := row["options"]:
            options = json.loads(row_options)

        # BACKPORT
        if backporting:
            if key.startswith("horak:"):
                # horak codes need to be changed to regular selects instead of point based criteria
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

            if key == "assertion_method":  # don't want to risk backporting multiple assertion methods
                value_type = EvidenceKeyValueType.SELECT

            if evidence_category == "SC":  # somatic clinical evidence doesn't exist in old version
                evidence_category = EvidenceCategory.INTERPRETATION
        # END BACKPORT
        else:
            # Forward port, provide evidence keys with newer properties
            e_key.default_crit_evaluation = row["default_crit_evaluation"]
            e_key.crit_allows_override_strengths = row["crit_allows_override_strengths"]
            e_key.crit_uses_points = row["crit_uses_points"]
            e_key.namespace_overrides = row["crit_uses_points"]

        e_key.evidence_category = evidence_category
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