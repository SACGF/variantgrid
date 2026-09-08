"""
Resyncs ReleaseGeneSymbolGene (the symbol -> gene matches for a GeneAnnotationRelease) with what
ReleaseGeneMatcher produces today: inserts missing matches, updates changed match_info and deletes
matches the current rules no longer make.

The delete is what makes this more than a rematch: matching only ever inserted, so rows written by
the old multi-hop alias traversal (#1669) and rows derived from alias rows since deleted by a
re-import survive a rematch. Run with --dry-run first to see the diff.
"""
from django.core.management import BaseCommand

from genes.gene_matching import ReleaseGeneMatcher
from genes.models import GeneAnnotationRelease, ReleaseGeneSymbol, ReleaseGeneSymbolGene

BATCH_SIZE = 2000


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report the diff without writing anything")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        verbosity = options["verbosity"]

        for gar in GeneAnnotationRelease.objects.all():
            release_gene_symbols = list(gar.releasegenesymbol_set.all())
            gm = ReleaseGeneMatcher(gar)
            expected = gm._get_gene_id_and_match_info_for_symbol(rgs.gene_symbol_id for rgs in release_gene_symbols)
            expected_matches = {}  # (gene_symbol_id, gene_id) -> match_info
            for gene_symbol_id, gene_id_and_match_info in expected.items():
                for gene_id, match_info in gene_id_and_match_info:
                    expected_matches[(gene_symbol_id, gene_id)] = match_info

            existing_qs = ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=gar)
            existing = existing_qs.values_list("pk", "release_gene_symbol__gene_symbol_id", "gene_id", "match_info")

            delete_pks = []
            update_records = []
            for pk, gene_symbol_id, gene_id, match_info in existing:
                key = (gene_symbol_id, gene_id)
                if key in expected_matches:
                    expected_match_info = expected_matches.pop(key)  # What's left over needs inserting
                    if expected_match_info != match_info:
                        update_records.append(ReleaseGeneSymbolGene(pk=pk, match_info=expected_match_info))
                else:
                    delete_pks.append(pk)
                    if verbosity >= 2:
                        print(f"{gar} - delete {gene_symbol_id} -> {gene_id} ({match_info})")

            num_insert = len(expected_matches)
            print(f"{gar} - delete: {len(delete_pks)}, update: {len(update_records)}, insert: {num_insert}")

            if not dry_run:
                for i in range(0, len(delete_pks), BATCH_SIZE):
                    ReleaseGeneSymbolGene.objects.filter(pk__in=delete_pks[i:i + BATCH_SIZE]).delete()
                if update_records:
                    ReleaseGeneSymbolGene.objects.bulk_update(update_records, ["match_info"], batch_size=BATCH_SIZE)
                if num_insert:
                    gm.match_symbols_to_genes(release_gene_symbols)

            no_match_qs = ReleaseGeneSymbol.objects.filter(release=gar, releasegenesymbolgene__isnull=True)
            print(f"{gar} - symbols w/o gene: {no_match_qs.count()}")
