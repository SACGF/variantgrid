"""
Fake classifications for 'manage.py create_fake_data' (@see library.fake_data): the "classifications" step here -
current records from both fake labs, some agreeing, some discordant, some withdrawn, for the classification listing,
discordance and allele pages - and the "reclassifications" step (curation histories) in reclassifications.py.
"""
import random

from django.contrib.auth.models import User

from classification.enums import CriteriaEvaluation, ShareLevel, SpecialEKeys, SubmissionSource
from classification.fake_data import reclassifications  # registers the "reclassifications" step
from classification.fake_data.shared import create_fake_allele_infos, delete_fake_classifications
from classification.models.classification import Classification
from library.fake_data import FakeData, FakeDataContext, register

LAB_RECORD_PREFIX = "fake-class-"

# Clinical significance for the one lab, or for both - the first pair agree, the second disagree across the
# pathogenic / VUS / benign buckets, which is what raises a discordance
SINGLE_LAB_CALLS = ["P", "LP", "VUS", "VUS", "LB", "B"]
AGREEING_CALLS = [("P", "LP"), ("VUS", "VUS"), ("B", "LB"), ("LP", "LP")]
DISCORDANT_CALLS = [("P", "VUS"), ("LP", "LB"), ("VUS", "B")]
CRITERIA_FOR_CALL = {
    "P": {"acmg:pvs1": CriteriaEvaluation.PATHOGENIC_VERY_STRONG, "acmg:pm2": CriteriaEvaluation.PATHOGENIC_MODERATE},
    "LP": {"acmg:ps3": CriteriaEvaluation.PATHOGENIC_STRONG, "acmg:pp3": CriteriaEvaluation.PATHOGENIC_SUPPORTING},
    "VUS": {"acmg:pm2": CriteriaEvaluation.PATHOGENIC_MODERATE},
    "LB": {"acmg:bp4": CriteriaEvaluation.BENIGN_SUPPORTING, "acmg:bs1": CriteriaEvaluation.BENIGN_STRONG},
    "B": {"acmg:ba1": CriteriaEvaluation.BENIGN_STANDALONE},
}
WITHDRAWN_CHANCE = 0.08

# Conditions are given as ontology IDs, which condition matching resolves locally - free text falls through to
# the Monarch search, an external call per distinct text that can never match "Fake condition for X"
CONDITION_FOR_GENE = {
    "APC": "MONDO:0016613", "ASXL1": "MONDO:0011510", "ATM": "MONDO:0700270", "BMPR1A": "MONDO:0012405",
    "BRAF": "MONDO:0007265", "BRCA1": "MONDO:0700268", "BRCA2": "MONDO:0700269", "CDH1": "MONDO:0100488",
    "CFTR": "MONDO:0009061", "CHEK2": "MONDO:0700271", "COL4A5": "MONDO:0010520", "DMD": "MONDO:0700285",
    "DNMT3A": "MONDO:0014382", "EGFR": "MONDO:0014481", "FBN1": "MONDO:0007514", "GATA2": "MONDO:0013607",
    "IDH1": "MONDO:0013808", "IDH2": "MONDO:0013345", "JAK2": "MONDO:0009891", "KIT": "MONDO:0008244",
    "KRAS": "MONDO:0012371", "LDLR": "MONDO:0007750", "MLH1": "MONDO:0012249", "MSH2": "MONDO:0007356",
    "MSH6": "MONDO:0013710", "MUTYH": "MONDO:0012041", "MYBPC3": "MONDO:0007268", "NF1": "MONDO:0018975",
    "NRAS": "MONDO:0013186", "PALB2": "MONDO:0012565", "PIK3CA": "MONDO:0013038", "PMS2": "MONDO:0013699",
    "PTEN": "MONDO:0008021", "RB1": "MONDO:0018160", "RET": "MONDO:0008082", "RUNX1": "MONDO:0100083",
    "RYR1": "MONDO:0007783", "SCN1A": "MONDO:0011461", "SDHB": "MONDO:0007273", "SMAD4": "MONDO:0007688",
    "STK11": "MONDO:0008280", "TET2": "MONDO:0030858", "TP53": "MONDO:0011775", "TTN": "MONDO:0013412",
    "VHL": "MONDO:0009892",
}
UNDECIDED_CONDITION_CHANCE = 0.15
""" Two terms with no uncertain/co-occurring choice, which leaves the text on the condition matching page """


@register
class FakeClassifications(FakeData):
    name = "classifications"
    help = ("Current classifications from both fake labs - some agreeing, some discordant, some withdrawn - "
            "for the classification listing, discordance and allele pages")
    requires = ("people", "variants")

    @classmethod
    def add_arguments(cls, parser):
        parser.add_argument("--alleles", type=int, default=40, help="Variants to classify")

    def create(self, context: FakeDataContext, **options):
        if Classification.objects.filter(lab_record_id__startswith=LAB_RECORD_PREFIX).exists():
            context.stdout.write("Fake classifications already exist")
            return

        rng = random.Random(context.seed)
        genes_and_variants = sorted((variant_id, gene_symbol)
                                    for gene_symbol, variant_ids in context.variant_ids_by_gene.items()
                                    for variant_id in variant_ids)
        chosen = rng.sample(genes_and_variants, min(options["alleles"], len(genes_and_variants)))
        allele_infos = create_fake_allele_infos(context.genome_build, chosen)

        germline_lab, somatic_lab = context.labs
        users_by_lab = {lab: list(User.objects.filter(groups__name=lab.group_name,
                                                      pk__in=[user.pk for user in context.users]))
                        for lab in context.labs}
        classifications = []
        for (_variant_id, gene_symbol), allele_info in zip(chosen, allele_infos):
            roll = rng.random()
            if roll < 0.4:
                calls = [(rng.choice(context.labs), rng.choice(SINGLE_LAB_CALLS))]
            else:
                pair = rng.choice(AGREEING_CALLS if roll < 0.75 else DISCORDANT_CALLS)
                calls = list(zip((germline_lab, somatic_lab), pair))

            for lab, call in calls:
                user = rng.choice(users_by_lab[lab])
                lab_record_id = f"{LAB_RECORD_PREFIX}{len(classifications) + 1}"
                condition = _condition(rng, gene_symbol)
                evidence = _evidence(context, gene_symbol, allele_info.imported_c_hgvs, call, condition)
                classification = Classification.create(user=user, lab=lab, lab_record_id=lab_record_id,
                                                       data=evidence, source=SubmissionSource.API,
                                                       allele_info=allele_info)
                classification.apply_allele_info_to_classification()
                classification.publish_latest(user, share_level=ShareLevel.ALL_USERS)
                if rng.random() < WITHDRAWN_CHANCE:
                    classification.set_withdrawn(user, withdraw=True)
                classifications.append(classification)

        num_withdrawn = sum(c.withdrawn for c in classifications)
        context.stdout.write(f"Created {len(classifications)} classifications ({num_withdrawn} withdrawn) "
                             f"over {len(chosen)} alleles")

    def delete(self, context: FakeDataContext, **options):
        classifications_qs = Classification.objects.filter(lab_record_id__startswith=LAB_RECORD_PREFIX)
        context.stdout.write(delete_fake_classifications(classifications_qs))


def _condition(rng: random.Random, gene_symbol: str) -> str:
    """ The gene's condition, or (sometimes, and always for a gene without one) two terms for a curator to decide """
    term_id = CONDITION_FOR_GENE.get(gene_symbol)
    if term_id and rng.random() >= UNDECIDED_CONDITION_CHANCE:
        return term_id
    other_term_ids = sorted(set(CONDITION_FOR_GENE.values()) - {term_id})
    return "; ".join(filter(None, [term_id, *rng.sample(other_term_ids, 2 if term_id is None else 1)]))


def _evidence(context: FakeDataContext, gene_symbol: str, c_hgvs: str, call: str, condition: str) -> dict:
    evidence = {
        SpecialEKeys.GENOME_BUILD: context.genome_build.name,
        SpecialEKeys.C_HGVS: c_hgvs,
        SpecialEKeys.GENE_SYMBOL: gene_symbol,
        SpecialEKeys.ALLELE_ORIGIN: "germline",
        SpecialEKeys.CLINICAL_SIGNIFICANCE: call,
        SpecialEKeys.CONDITION: condition,
        "interpretation_summary": f"Fake {call} call on a {gene_symbol} variant",
    }
    evidence.update(CRITERIA_FOR_CALL[call])
    return {key: {"value": value} for key, value in evidence.items()}
