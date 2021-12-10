from collections import defaultdict
from typing import Dict

from django.core.management.base import BaseCommand

from classification.models import ConditionText, sync_all_condition_resolutions_to_classifications
from classification.models.condition_text_matching import ConditionTextMatch


class Command(BaseCommand):
    """
    Updates Condition Text Matches to match the state of the classifications
    """

    def add_arguments(self, parser):
        parser.add_argument('--reset', action='store_true', default=False)
        parser.add_argument('--clear', action='store_true', default=False)
        parser.add_argument('--classifications', action='store_true', default=False)
        parser.add_argument('--orphans', actions='store_true', default=False)

    def handle(self, *args, **options):
        if options["classifications"]:
            print("Updating classifications")
            sync_all_condition_resolutions_to_classifications()
            print("Complete")
            return
        if options["orphans"]:
            print("Removing orphans")
            for ct in ConditionText.objects.all():
                classifications_for_gene_symbol: Dict[str, int] = defaultdict(int)
                for ctm in list(ct.conditiontextmatch_set.all()):
                    if gene_symbol := ctm.gene_symbol:
                        if ctm.classification_id:
                            classifications_for_gene_symbol[gene_symbol] += 1
                        else:
                            classifications_for_gene_symbol[gene_symbol] += 0  # register that we have the gene symbol
                for gene_symbol, count in classifications_for_gene_symbol.items():
                    if count == 0:
                        ct.conditiontextmatch_set.filter(gene_symbol=gene_symbol).delete()
            return

        if options["reset"]:
            print("Deleting old records")
            ConditionText.objects.all().delete()
        if options["clear"]:
            print("Clearing all existing values")
            ct: ConditionText
            for ct in ConditionText.objects.all():
                ct.clear()

        print("Syncing")
        ConditionTextMatch.sync_all()
        print("Complete")
