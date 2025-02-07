"""
    This should be run after variant normalization + splitting multi-alts
    so we can remove only the bad alt not all on the line
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
    ALT_CONVERSION = {
        # This should all be in upper case
        "<DUP:TANDEM>": "<DUP>",
    }

    def add_arguments(self, parser):
        parser.add_argument('--vcf', help='VCF file, default: - (stdin)', default="-")
        parser.add_argument('--skipped-records-stats-file', help='File name')
        parser.add_argument('--converted-records-stats-file', help='File name')

    def handle(self, *args, **options):
        vcf_filename = options["vcf"]
        skipped_records_stats_file = options.get("skipped_records_stats_file")
        converted_records_stats_file = options.get("converted_records_stats_file")

        if vcf_filename == '-':
            f = sys.stdin
        else:
            f = open(vcf_filename)

        skipped_records = Counter()
        converted_records = Counter()

        alt_standard_bases_pattern = re.compile(r"[GATC,.]")  # Can "." for reference
        skipped_alt_first_examples = None

        for line in f:
            if line[0] == '#':
                sys.stdout.write(line)
            else:
                columns: list[str] = line.split("\t")
                # Check alt is ok
                alt = columns[VCFColumns.ALT].upper()
                if to_value := self.ALT_CONVERSION.get(alt):
                    converted_msg = f"'{alt=}' to '{to_value}'"
                    converted_records[converted_msg] += 1
                    columns[VCFColumns.ALT] = to_value
                elif alt_standard_bases_pattern.sub("", alt):
                    skip_reason = None
                    if alt.startswith("<") and alt.endswith(">"):
                        if settings.VARIANT_SYMBOLIC_ALT_ENABLED:
                            if alt not in settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES:
                                skip_reason = f"ALT = {alt}"
                        else:
                            skip_reason = "Symbolic variants disabled via settings."
                    else:
                        if not skipped_alt_first_examples:
                            skipped_alt_first_examples = f"non-standard bases in ALT sequence. First example: '{alt}'"
                        skip_reason = skipped_alt_first_examples

                    if skip_reason:
                        skipped_records[skip_reason] += 1
                        continue

                sys.stdout.write("\t".join(columns))

        self._write_vcf_stats(skipped_records, skipped_records_stats_file)
        self._write_vcf_stats(converted_records, converted_records_stats_file)

    @staticmethod
    def _write_vcf_stats(counts, filename):
        if counts and filename:
            with open(filename, "w") as f:
                for name, count in counts.items():
                    f.write('%s\t%d\n' % (name, count))
