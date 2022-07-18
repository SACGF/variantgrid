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
        parser.add_argument('--obsolete', action='store_true', default=False)
        parser.add_argument('--reset', action='store_true', default=False)
        parser.add_argument('--clear', action='store_true', default=False)
        parser.add_argument('--classifications', action='store_true', default=False)
        parser.add_argument('--orphans', action='store_true', default=False)

    def check_obsoletes(self):
        for ctm in ConditionTextMatch.objects.filter(condition_xrefs__isnull=False):
            for term in ctm.condition_xref_terms:
                if term.is_obsolete:
                    print(
                        f"ConditionTextMatch ({ctm.condition_text.lab.name} '{ctm.condition_text.normalized_text}' {ctm.condition_text.pk}) refers to obsolete term {term}")

    def handle(self, *args, **options):
        if options["obsolete"]:
            print("Checking for obsoletes")
            self.check_obsoletes()
            print("Complete")
            return
        if options["classifications"]:
            print("Updating classifications")
            sync_all_condition_resolutions_to_classifications()
            print("Complete")
            return
        if options["orphans"]:
            print("Removing orphans")
            total_deleted = 0
            ct_total_deleted = 0
            for ct in ConditionText.objects.all():
                classifications_for_gene_symbol: Dict[str, int] = defaultdict(int)
                for ctm in list(ct.conditiontextmatch_set.all()):
                    if gene_symbol := ctm.gene_symbol:
                        if ctm.classification_id:
                            classifications_for_gene_symbol[gene_symbol] += 1
                        else:
                            classifications_for_gene_symbol[gene_symbol] += 0  # register that we have the gene symbol
                found_records = False
                for gene_symbol, count in classifications_for_gene_symbol.items():
                    if count == 0:
                        delete_qs = ct.conditiontextmatch_set.filter(gene_symbol=gene_symbol)
                        delete_count = delete_qs.count()
                        total_deleted += delete_count
                        delete_qs.delete()
                        found_records = True

                if not ct.conditiontextmatch_set.filter(classification__isnull=False).exists():
                    found_records += 1
                    ct.delete()
                    ct_total_deleted += 1

                if found_records:
                    print(f"Total parent / child records deleted - {ct_total_deleted} / {total_deleted}")
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
