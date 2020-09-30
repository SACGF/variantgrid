from django.core.management.base import BaseCommand
import itertools
from random import randrange
import sys

from annotation.vcf_files.vcf_export_utils import get_vcf_header_from_contigs
from snpdb.models import GenomeBuild

MAX_SIZE = 1000000


def random_contig():
    human_contigs = itertools.chain(range(1, 22), ['X', 'Y'])
    contigs = [str(i) for i in human_contigs]
    i = randrange(len(contigs))
    return contigs[i]


def random_base(not_base=None):
    bases = 'GATC'

    while True:
        i = randrange(len(bases))
        base = bases[i]
        if base != not_base:
            return base


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument('--num-variants', type=int, default=10)

    def handle(self, *args, **options):
        f = sys.stdout

        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        header_lines = get_vcf_header_from_contigs(genome_build, {})
        for line in header_lines:
            f.write(line + '\n')

        num_variants = options["num_variants"]
        rows = []
        for _ in range(num_variants):
            ref = random_base()
            alt = random_base(not_base=ref)
            chrom = random_contig()
            pos = randrange(MAX_SIZE)
            rows.append([chrom, pos, '.', ref, alt, '.', '.', '.'])

        rows = sorted(rows, key=lambda x: (x[0], int(x[1])))
        for data in rows:
            line = '\t'.join(map(str, data))
            f.write(line + '\n')
