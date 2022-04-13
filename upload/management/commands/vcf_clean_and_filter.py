"""
Uploaded VCFs are first passed through this command to fix things that will cause VT to die (eg bad contigs)
"""
import re
import sys
from collections import Counter

from django.conf import settings
from django.core.management.base import BaseCommand

from snpdb.models import GenomeBuild, GenomeFasta


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--vcf', help='VCF file, default: - (stdin)', default="-")
        parser.add_argument('--genome-build', help='GenomeBuild name', required=True)
        parser.add_argument('--skipped-contigs-stats-file', help='File name')
        parser.add_argument('--skipped-records-stats-file', help='File name')

    def handle(self, *args, **options):
        vcf_filename = options["vcf"]
        build_name = options["genome_build"]
        skipped_contigs_stats_file = options.get("skipped_contigs_stats_file")
        skipped_records_stats_file = options.get("skipped_records_stats_file")

        genome_build = GenomeBuild.get_name_or_alias(build_name)
        genome_fasta = GenomeFasta.get_for_genome_build(genome_build)
        chrom_to_contig_id = genome_build.get_chrom_contig_id_mappings()
        contig_to_fasta_names = genome_fasta.get_contig_id_to_name_mappings()

        if vcf_filename == '-':
            f = sys.stdin
        else:
            f = open(vcf_filename)

        skipped_contigs = Counter()
        skipped_records = Counter()

        ref_standard_bases_pattern = re.compile("[GATC]")
        alt_standard_bases_pattern = re.compile("[GATC,\.]")  # Can be multi-alts, or "." for reference

        skip_patterns = {}
        if skip_regex := getattr(settings, "VCF_IMPORT_SKIP_RECORD_REGEX", {}):
            for name, regex in skip_regex.items():
                skip_patterns[name] = re.compile(regex)

        for line in f:
            if line and line[0] != '#':
                columns = line.split("\t")
                chrom = columns[0]
                contig_id = chrom_to_contig_id.get(chrom)
                fasta_chrom = None
                if contig_id:
                    fasta_chrom = contig_to_fasta_names.get(contig_id)

                if not fasta_chrom:
                    if skipped_contigs_stats_file:
                        skipped_contigs[chrom] += 1
                    continue

                if skip_patterns:
                    if skip_reason := self._check_skip_line(skip_patterns, line):
                        skipped_records[skip_reason] += 1
                        continue

                # Check ref / alt bases are ok
                if settings.VARIANT_STANDARD_BASES_ONLY:
                    ref = columns[3]
                    alt = columns[4]
                    if ref_standard_bases_pattern.sub("", ref):
                        if ref.startswith("<") and ref.endswith(">"):
                            skip_reason = f"REF = {ref}"
                        else:
                            skip_reason = "non-standard bases in REF sequence"
                        skipped_records[skip_reason] += 1
                        continue

                    if alt_standard_bases_pattern.sub("", alt):
                        if alt.startswith("<") and alt.endswith(">"):
                            skip_reason = f"ALT = {alt}"
                        else:
                            skip_reason = "non-standard bases in ALT sequence"
                        skipped_records[skip_reason] += 1
                        continue

                columns[0] = fasta_chrom
                # Zero out INFO (makes file size smaller and causes bcftools issues)
                columns[7] = "."
                # If (7) INFO was the last column, we just stripped the newline - might need to add it back
                if len(columns) == 8:
                    columns[7] += "\n"
                sys.stdout.write("\t".join(columns))
            else:
                sys.stdout.write(line)

        self._write_skip_counts(skipped_contigs, skipped_contigs_stats_file)
        self._write_skip_counts(skipped_records, skipped_records_stats_file)

    @staticmethod
    def _write_skip_counts(counts, filename):
        if counts and filename:
            with open(filename, "w") as f:
                for name, count in counts.items():
                    f.write('%s\t%d\n' % (name, count))

    @staticmethod
    def _check_skip_line(skip_patterns: dict, line: str):
        for name, pattern in skip_patterns.items():
            if pattern.findall(line):
                return name
        return None
