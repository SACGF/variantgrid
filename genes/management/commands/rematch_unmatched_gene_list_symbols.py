from django.core.management import BaseCommand

from genes.gene_matching import GeneSymbolMatcher
from genes.models import GeneListGeneSymbol


class Command(BaseCommand):
    def handle(self, *args, **options):
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
