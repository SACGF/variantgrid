"""
Uploaded VCFs are first passed through this command to fix things that will cause VT to die (eg bad contigs)

We also rename chromosomes to the RefSeq contig IDs used in our fasta file

We don't use "bcftools annotate --rename-chrs" to do this as it fails if contigs are not defined in the header
Which we unfortunately see sometimes.

"""
import re
import sys
from collections import Counter
from io import StringIO

from django.conf import settings
from django.core.management.base import BaseCommand
from vcf import Reader

from library.genomics.vcf_enums import VCFColumns
from snpdb.models import GenomeBuild, GenomeFasta


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--vcf', help='VCF file, default: - (stdin)', default="-")
        parser.add_argument('--replace-header', help='VCF header file', required=False)
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
        replace_header = options["replace_header"]
        skipped_contigs_stats_file = options.get("skipped_contigs_stats_file")
        skipped_records_stats_file = options.get("skipped_records_stats_file")
        skipped_filters_stats_file = options.get("skipped_filters_stats_file")

        genome_build = GenomeBuild.get_name_or_alias(build_name)
        genome_fasta = GenomeFasta.get_for_genome_build(genome_build)
        chrom_to_contig_id = genome_build.get_chrom_contig_id_mappings()
        contig_lengths = dict(genome_build.contigs.values_list("pk", "length"))
        contig_to_fasta_names = genome_fasta.get_contig_id_to_name_mappings()

        if vcf_filename == '-':
            f = sys.stdin
        else:
            f = open(vcf_filename)

        skipped_contigs = Counter()
        skipped_records = Counter()
        skipped_filters = Counter()

        ref_standard_bases_pattern = re.compile(r"[GATCN]")  # Reference can be N (and FreeBayes often writes these)
        alt_standard_bases_pattern = re.compile(r"[GATC,\.]")  # Can be multi-alts, or "." for reference

        skip_patterns = {}
        if skip_regex := getattr(settings, "VCF_IMPORT_SKIP_RECORD_REGEX", {}):
            for name, regex in skip_regex.items():
                skip_patterns[name] = re.compile(regex)

        vcf_header_lines = []
        if replace_header:
            with open(replace_header, "r") as rh_f:
                vcf_header_lines = rh_f.readlines()

        first_non_header_line = True
        defined_filters = None
        for line in f:
            if line[0] == '#':
                if not replace_header:
                    vcf_header_lines.append(line)
            else:
                if first_non_header_line:
                    # Dump out the header
                    for header_line in vcf_header_lines:
                        sys.stdout.write(header_line)
                    defined_filters = self._get_defined_vcf_filters(vcf_header_lines)
                    first_non_header_line = False

                columns = line.split("\t")
                chrom = columns[VCFColumns.CHROM]
                fasta_chrom = None
                if contig_id := chrom_to_contig_id.get(chrom):
                    valid_position = False
                    try:
                        position = int(columns[VCFColumns.POS])
                        length = contig_lengths[contig_id]
                        valid_position = 0 < position < length
                    except ValueError:
                        pass
                    if not valid_position:
                        skipped_records["position out of range"] += 1
                        continue
                    fasta_chrom = contig_to_fasta_names.get(contig_id)

                if not fasta_chrom:
                    if skipped_contigs_stats_file:
                        skipped_contigs[chrom] += 1
                    continue
                columns[VCFColumns.CHROM] = fasta_chrom

                if skip_patterns:
                    if skip_reason := self._check_skip_line(skip_patterns, line):
                        skipped_records[skip_reason] += 1
                        continue

                # Check ref / alt bases are ok
                ref = columns[VCFColumns.REF]
                alt = columns[VCFColumns.ALT]
                if ref_standard_bases_pattern.sub("", ref):
                    if ref.startswith("<") and ref.endswith(">"):
                        skip_reason = f"REF = {ref}"
                    else:
                        skip_reason = "non-standard bases in REF sequence"
                    skipped_records[skip_reason] += 1
                    continue

                if alt_standard_bases_pattern.sub("", alt):
                    skip_reason = None
                    if alt.startswith("<") and alt.endswith(">"):
                        if alt not in settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES:
                            skip_reason = f"ALT = {alt}"
                    else:
                        skip_reason = "non-standard bases in ALT sequence"
                    if skip_reason:
                        skipped_records[skip_reason] += 1
                        continue

                # Remove filters not in header
                filter_column = columns[VCFColumns.FILTER]
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
                    columns[VCFColumns.FILTER] = filter_column

                if remove_info:
                    # Zero out INFO (makes file size smaller and causes bcftools issues)
                    columns[VCFColumns.INFO] = "."
                    # If (7) INFO was the last column, we just stripped the newline - might need to add it back
                    if len(columns) == 8:
                        columns[VCFColumns.INFO] += "\n"
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
