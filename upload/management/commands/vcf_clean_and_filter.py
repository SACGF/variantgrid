"""
Uploaded VCFs are first passed through this command to fix things that will cause bcftools to die (eg bad contigs)

It also fails the import early on unsorted VCFs, which would otherwise only be noticed by 'bcftools index'
several pipeline stages later. When '--unsorted-marker-file' is given, the reason is written there before
failing, so the caller can tell an unsorted VCF apart from any other failure and retry through 'bcftools sort'
(pass '--allow-unsorted' on that second pass).

We also rename chromosomes to the RefSeq contig IDs used in our fasta file

We don't use "bcftools annotate --rename-chrs" to do this as it fails if contigs are not defined in the header
Which we unfortunately see sometimes.

"""
import re
import sys
from collections import Counter

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from library.genomics.vcf_enums import VCFColumns
from library.genomics.vcf_utils import UnsortedVCFError, VCFSortChecker, vcf_header_filter_ids
from snpdb.models import GenomeBuild, GenomeFasta


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--vcf', help='VCF file, default: - (stdin)', default="-")
        parser.add_argument('--replace-header', help='VCF header file', required=False)
        parser.add_argument('--genome-build', help='GenomeBuild name', required=True)
        parser.add_argument('--skipped-contigs-stats-file', help='File name')
        parser.add_argument('--skipped-records-stats-file', help='File name')
        parser.add_argument('--skipped-filters-stats-file', help='File name')
        parser.add_argument('--unsorted-marker-file', help='Write reason here if VCF is unsorted')
        parser.add_argument('--allow-unsorted', action='store_true', help='Skip the sorted check')

    @staticmethod
    def _write_header(vcf_header_lines: list[str]):
        for header_line in vcf_header_lines:
            sys.stdout.write(header_line)

    def handle(self, *args, **options):
        QUICK_ACCEPT_FILTERS = {".", "PASS"}
        vcf_filename = options["vcf"]
        build_name = options["genome_build"]
        replace_header = options["replace_header"]
        skipped_contigs_stats_file = options.get("skipped_contigs_stats_file")
        skipped_records_stats_file = options.get("skipped_records_stats_file")
        skipped_filters_stats_file = options.get("skipped_filters_stats_file")
        unsorted_marker_file = options.get("unsorted_marker_file")
        allow_unsorted = options.get("allow_unsorted")

        genome_build = GenomeBuild.get_name_or_alias(build_name)
        genome_fasta = GenomeFasta.get_for_genome_build(genome_build)
        chrom_to_contig_id = genome_build.get_chrom_contig_id_mappings(standard_contigs_only=True)
        contig_lengths = {contig_id: int(length) for contig_id, length in genome_build.contigs.values_list("pk", "length")}
        contig_to_fasta_names = genome_fasta.get_contig_id_to_name_mappings()

        if vcf_filename == '-':
            f = sys.stdin
        else:
            f = open(vcf_filename)

        skipped_contigs = Counter()
        skipped_records = Counter()
        skipped_filters = Counter()

        ref_standard_bases_pattern = re.compile(r"[GATCN]")  # Reference can be N (and FreeBayes often writes these)

        skip_patterns = {}
        if skip_regex := getattr(settings, "VCF_IMPORT_SKIP_RECORD_REGEX", {}):
            for name, regex in skip_regex.items():
                skip_patterns[name] = re.compile(regex)

        vcf_header_lines = []
        if replace_header:
            with open(replace_header) as rh_f:
                vcf_header_lines = rh_f.readlines()

        first_non_header_line = True
        defined_filters = None
        sort_checker = None if allow_unsorted else VCFSortChecker()
        for line_number, line in enumerate(f, start=1):
            if line[0] == '#':
                if not replace_header:
                    vcf_header_lines.append(line)
            else:
                if first_non_header_line:
                    self._write_header(vcf_header_lines)
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

                # Check ref bases are ok
                # Alts are checked in 'vcf_remove_non_standard_alts', which happens after split multi-allelic
                ref = columns[VCFColumns.REF]
                if ref_standard_bases_pattern.sub("", ref):
                    if ref.startswith("<") and ref.endswith(">"):
                        skip_reason = f"REF = {ref}"
                    else:
                        skip_reason = "non-standard bases in REF sequence"
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

                # Only the records we write need to be sorted - skipped ones can't break the downstream index
                if sort_checker:
                    try:
                        sort_checker.check(contig_id, position, chrom, line_number)
                    except UnsortedVCFError as e:
                        reason = f"VCF is not sorted. {e}"
                        self._write_unsorted_marker(unsorted_marker_file, reason)
                        raise CommandError(reason) from e

                sys.stdout.write("\t".join(columns))

        if first_non_header_line:
            # Handle case where there were no records at all - still need to write header
            self._write_header(vcf_header_lines)

        self._write_skip_counts(skipped_contigs, skipped_contigs_stats_file)
        self._write_skip_counts(skipped_records, skipped_records_stats_file)
        self._write_skip_counts(skipped_filters, skipped_filters_stats_file)

    @staticmethod
    def _get_defined_vcf_filters(vcf_header_lines) -> set:
        return {"PASS"} | vcf_header_filter_ids(vcf_header_lines)

    @staticmethod
    def _write_unsorted_marker(filename, reason: str):
        if filename:
            with open(filename, "w") as f:
                f.write(reason + "\n")

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
