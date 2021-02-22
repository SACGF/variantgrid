from collections import defaultdict

from django.core.management import BaseCommand
from django.db.models import Q

from genes.gene_matching import GeneSymbolMatcher
from genes.models import GeneListGeneSymbol


class Command(BaseCommand):
    def handle(self, *args, **options):
        # Some legacy gene list may not have been stripped properly.
        stripped = []
        q_strip = Q(original_name__startswith=' ') | Q(original_name__endswith=' ')
        gene_lists = set()
        for glgs in GeneListGeneSymbol.objects.filter(q_strip):
            glgs.original_name = glgs.original_name.strip()
            stripped.append(glgs)
            gene_lists.add(glgs.gene_list)

        if stripped:
            print(f"{len(stripped)} records stripped")
            # Need to make sure that stripped record doesn't already exist
            original_names_per_list = defaultdict(set)
            glgs_qs = GeneListGeneSymbol.objects.filter(gene_list__in=gene_lists)
            for gene_list_id, original_name in glgs_qs.values_list("gene_list_id", "original_name"):
                original_names_per_list[gene_list_id].add(original_name)

            stripped_ok = []
            dupes_to_delete = []
            for glg in stripped:
                if glg.original_name in original_names_per_list[glg.gene_list_id]:
                    dupes_to_delete.append(glg.pk)
                else:
                    stripped_ok.append(glg)

            if dupes_to_delete:
                print(f"{len(dupes_to_delete)} records stripped to dupes, deleting")
                GeneListGeneSymbol.objects.filter(pk__in=dupes_to_delete).delete()
            GeneListGeneSymbol.objects.bulk_update(stripped_ok,
                                                   fields=["original_name"],
                                                   batch_size=2000)

        gsm = GeneSymbolMatcher()
        modified_records = []
        num_unmatched = 0
        for glgs in GeneListGeneSymbol.objects.filter(gene_symbol__isnull=True):
            gene_symbol_id, alias_id = gsm.get_gene_symbol_id_and_alias_id(glgs.original_name)
            if gene_symbol_id:
                glgs.gene_symbol_id = gene_symbol_id
                glgs.gene_symbol_alias_id = alias_id
                modified_records.append(glgs)
            else:
                num_unmatched += 1

        print(f"Matched {len(modified_records)}, Unmatched: {num_unmatched}")
        if modified_records:
            GeneListGeneSymbol.objects.bulk_update(modified_records,
                                                   fields=["gene_symbol_id", "gene_symbol_alias"],
                                                   batch_size=2000)
            print("Matching to Gene Annotation Release")
            gsm._match_symbols_to_genes_in_releases()
