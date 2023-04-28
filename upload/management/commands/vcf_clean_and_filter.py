"""
Uploaded VCFs are first passed through this command to fix things that will cause VT to die (eg bad contigs)

"""
import re
import sys
from collections import Counter
from io import StringIO

from django.conf import settings
from django.core.management.base import BaseCommand
from vcf import Reader

from snpdb.models import GenomeBuild


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--vcf', help='VCF file, default: - (stdin)', default="-")
        parser.add_argument('--genome-build', help='GenomeBuild name', required=True)
        parser.add_argument('--remove-info', action='store_true', help='clear INFO field')
        parser.add_argument('--skipped-contigs-stats-file', help='File name')
        parser.add_argument('--skipped-records-stats-file', help='File name')
        parser.add_argument('--skipped-filters-stats-file', help='File name')

    def handle(self, *args, **options):
        QUICK_ACCEPT_FILTERS = {".", "PASS"}
        vcf_filename = options["vcf"]
        build_name = options["genome_build"]
        remove_info = options["remove_info"]
        skipped_contigs_stats_file = options.get("skipped_contigs_stats_file")
        skipped_records_stats_file = options.get("skipped_records_stats_file")
        skipped_filters_stats_file = options.get("skipped_filters_stats_file")

        genome_build = GenomeBuild.get_name_or_alias(build_name)
        chrom_to_contig_id = genome_build.get_chrom_contig_id_mappings()
        contig_lengths = dict(genome_build.contigs.values_list("pk", "length"))

        if vcf_filename == '-':
            f = sys.stdin
        else:
            f = open(vcf_filename)

        skipped_contigs = Counter()
        skipped_records = Counter()
        skipped_filters = Counter()

        skip_patterns = {}
        if skip_regex := getattr(settings, "VCF_IMPORT_SKIP_RECORD_REGEX", {}):
            for name, regex in skip_regex.items():
                skip_patterns[name] = re.compile(regex)

        vcf_header_lines = []
        defined_filters = None
        for line in f:
            if line[0] == '#':
                vcf_header_lines.append(line)
                sys.stdout.write(line)
            else:
                if defined_filters is None:
                    defined_filters = self._get_defined_vcf_filters(vcf_header_lines)
                columns = line.split("\t")
                chrom = columns[0]
                if contig_id := chrom_to_contig_id.get(chrom):
                    valid_position = False
                    try:
                        position = int(columns[1])
                        length = contig_lengths[contig_id]
                        valid_position = 0 < position < length
                    except ValueError:
                        pass
                    if not valid_position:
                        skipped_records["position out of range"] += 1
                        continue
                else:
                    if skipped_contigs_stats_file:
                        skipped_contigs[chrom] += 1
                    continue

                if skip_patterns:
                    if skip_reason := self._check_skip_line(skip_patterns, line):
                        skipped_records[skip_reason] += 1
                        continue

                # Remove filters not in header
                filter_column = columns[6]
                if filter_column not in QUICK_ACCEPT_FILTERS:
                    cleaned_filters = []
                    for fc in filter_column.split(";"):
                        if fc in defined_filters:
                            cleaned_filters.append(fc)
                        else:
                            skipped_filters[fc] += 1

                    if cleaned_filters:
                        filter_column = ";".join(cleaned_filters)
                    else:
                        filter_column = "."
                    columns[6] = filter_column

                if remove_info:
                    # Zero out INFO (makes file size smaller and causes bcftools issues)
                    columns[7] = "."
                    # If (7) INFO was the last column, we just stripped the newline - might need to add it back
                    if len(columns) == 8:
                        columns[7] += "\n"
                sys.stdout.write("\t".join(columns))

        self._write_skip_counts(skipped_contigs, skipped_contigs_stats_file)
        self._write_skip_counts(skipped_records, skipped_records_stats_file)
        self._write_skip_counts(skipped_filters, skipped_filters_stats_file)

    @staticmethod
    def _get_defined_vcf_filters(vcf_header_lines) -> set:
        defined_filters = {"PASS"}
        stream = StringIO("".join(vcf_header_lines))
        reader = Reader(stream)
        defined_filters.update(reader.filters.keys())
        return defined_filters

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
