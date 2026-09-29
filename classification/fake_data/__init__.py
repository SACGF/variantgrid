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
                evidence = _evidence(context, gene_symbol, allele_info.imported_c_hgvs, call)
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


def _evidence(context: FakeDataContext, gene_symbol: str, c_hgvs: str, call: str) -> dict:
    evidence = {
        SpecialEKeys.GENOME_BUILD: context.genome_build.name,
        SpecialEKeys.C_HGVS: c_hgvs,
        SpecialEKeys.GENE_SYMBOL: gene_symbol,
        SpecialEKeys.ALLELE_ORIGIN: "germline",
        SpecialEKeys.CLINICAL_SIGNIFICANCE: call,
        SpecialEKeys.CONDITION: f"Fake condition for {gene_symbol}",
        "interpretation_summary": f"Fake {call} call on a {gene_symbol} variant",
    }
    evidence.update(CRITERIA_FOR_CALL[call])
    return {key: {"value": value} for key, value in evidence.items()}
