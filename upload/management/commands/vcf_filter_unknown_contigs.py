"""

Uploaded VCFs are first passed through this command to make sure contigs are ok

(VT will die on bad contigs)

"""
from collections import Counter
from django.core.management.base import BaseCommand
import sys

from snpdb.models import GenomeBuild, GenomeFasta


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--vcf', help='VCF file, default: - (stdin)', default="-")
        parser.add_argument('--genome-build', help='GenomeBuild name', required=True)
        parser.add_argument('--unknown-contigs-stats-file', help='File name to write stats')

    def handle(self, *args, **options):
        vcf_filename = options["vcf"]
        build_name = options["genome_build"]
        skipped_contigs_stats_file = options.get("unknown_contigs_stats_file")

        genome_build = GenomeBuild.get_name_or_alias(build_name)
        genome_fasta = GenomeFasta.get_for_genome_build(genome_build)
        chrom_to_contig_id = genome_build.get_chrom_contig_id_mappings()
        contig_to_fasta_names = genome_fasta.get_contig_id_to_name_mappings()

        if vcf_filename == '-':
            f = sys.stdin
        else:
            f = open(vcf_filename)

        skipped_contigs = Counter()
        for line in f:
            if line and line[0] != '#':
                chrom, rest_of_line = line.split("\t", 1)
                contig_id = chrom_to_contig_id.get(chrom)
                fasta_chrom = None
                if contig_id:
                    fasta_chrom = contig_to_fasta_names.get(contig_id)

                if not fasta_chrom:
                    if skipped_contigs_stats_file:
                        skipped_contigs[chrom] += 1
                    continue

                sys.stdout.write("\t".join([fasta_chrom, rest_of_line]))
            else:
                sys.stdout.write(line)

        if skipped_contigs_stats_file:
            with open(skipped_contigs_stats_file, "w") as f:
                for contig, count in skipped_contigs.items():
                    f.write('%s\t%d\n' % (contig, count))
