from django.core.management import BaseCommand

from genes.enums import AnnotationConsortium
from genes.gene_matching import ReleaseGeneMatcher
from genes.models import GeneAnnotationRelease, GeneSymbol, GeneVersion

"""
One off fix - can be deleted once all legacy environments migrate past genes.migrations.0015

Ensures:
    * HGNCGeneNames - symbols are matched to genes in each annotation release
    * Gene.hgnc is populated
"""


class Command(BaseCommand):
    category = "one-off"

    def handle(self, *args, **options):

        # Make sure gene symbols are matched to genes in each release
        for release in GeneAnnotationRelease.objects.all():
            gm = ReleaseGeneMatcher(release)
            gm.match_unmatched_in_hgnc_and_gene_lists()

        for annotation_consortium in [AnnotationConsortium.REFSEQ, AnnotationConsortium.ENSEMBL]:
            modified_records = []
            num_unmatched_hgnc = 0
            gene_versions = list(GeneVersion.objects.filter(gene__annotation_consortium=annotation_consortium))
            resolutions = GeneSymbol.resolve_names(gv.gene_symbol_id for gv in gene_versions)
            for gv in gene_versions:
                resolution, _ambiguous = resolutions[gv.gene_symbol_id]
                if hgnc := resolution.hgnc:
                    gv.hgnc = hgnc
                    modified_records.append(gv)
                else:
                    num_unmatched_hgnc += 1

            num_matched = len(modified_records)
            print(f"{AnnotationConsortium(annotation_consortium).label} - HGNC matched: {num_matched}, unmatched: {num_unmatched_hgnc}")
            if num_matched:
                GeneVersion.objects.bulk_update(modified_records, ["hgnc"], 2000)
