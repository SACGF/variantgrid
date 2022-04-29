from django.core.management import BaseCommand

from genes.gene_matching import ReleaseGeneMatcher, HGNCMatcher
from genes.models import GeneAnnotationRelease, GeneVersion
from genes.models_enums import AnnotationConsortium

"""
One off fix - can be deleted once all legacy environments migrate past genes.migrations.0015

Ensures:
    * HGNCGeneNames - symbols are matched to genes in each annotation release
    * Gene.hgnc is populated
"""


class Command(BaseCommand):
    def handle(self, *args, **options):

        # Make sure gene symbols are matched to genes in each release
        for release in GeneAnnotationRelease.objects.all():
            gm = ReleaseGeneMatcher(release)
            gm.match_unmatched_in_hgnc_and_gene_lists()

        hgnc_matcher = HGNCMatcher()
        for annotation_consortium in [AnnotationConsortium.REFSEQ, AnnotationConsortium.ENSEMBL]:
            modified_records = []
            num_unmatched_hgnc = 0
            for gv in GeneVersion.objects.filter(gene__annotation_consortium=annotation_consortium):
                if hgnc := hgnc_matcher.match_hgnc(gv.gene_symbol_id):
                    gv.hgnc = hgnc
                    modified_records.append(gv)
                else:
                    num_unmatched_hgnc += 1

            num_matched = len(modified_records)
            print(f"{AnnotationConsortium(annotation_consortium).label} - HGNC matched: {num_matched}, unmatched: {num_unmatched_hgnc}")
            if num_matched:
                GeneVersion.objects.bulk_update(modified_records, ["hgnc"], 2000)
