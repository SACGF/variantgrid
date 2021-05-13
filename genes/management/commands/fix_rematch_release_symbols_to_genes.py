from django.core.management import BaseCommand
from genes.models import GeneAnnotationRelease, ReleaseGeneSymbolGene, ReleaseGeneSymbol
from genes.gene_matching import GeneMatcher
from library.log_utils import log_traceback


class Command(BaseCommand):
    def handle(self, *args, **options):
        for gar in GeneAnnotationRelease.objects.all():
            no_match_qs = ReleaseGeneSymbol.objects.filter(release=gar, releasegenesymbolgene__isnull=True)

            qs = ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=gar)
            num_genes_original = qs.count()
            num_no_match_original = no_match_qs.count()
            print(f"{gar} - symbols w/o gene: {num_no_match_original}")
            gm = GeneMatcher(gar)
            release_gene_symbols = gar.releasegenesymbol_set.all()
            gm.match_symbols_to_genes(release_gene_symbols)
            num_genes = qs.count()

            print(f"{gar} - matched {num_genes - num_genes_original} genes")
            print(f"{gar} - {no_match_qs.count() - num_no_match_original} less symbols w/o genes")


