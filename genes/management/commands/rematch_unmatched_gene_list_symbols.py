from django.core.management import BaseCommand
from django.db.models import Q

from genes.gene_matching import GeneSymbolMatcher
from genes.models import GeneListGeneSymbol


class Command(BaseCommand):
    def handle(self, *args, **options):
        # Some legacy gene list may not have been stripped properly.
        stripped = []
        q_strip = Q(gene_symbol__startswith=' ') | Q(gene_symbol__endswith=' ')
        for glgs in GeneListGeneSymbol.objects.filter(q_strip):
            glgs.original_name = glgs.original_name.strip()
            stripped.append(glgs)
        if stripped:
            print(f"{len(stripped)} records stripped")
            GeneListGeneSymbol.objects.bulk_update(stripped,
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
