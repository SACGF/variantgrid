"""
Reading what VEP told us it skipped (#1701).

Two sources, because a run we are annotating now and a run annotated months ago leave behind
different things:

    --skipped_variants_file   one tab-separated row per skipped record, written by VEP at the very end
                              of its process (Runner::finish -> BaseRunner::dump_skipped_variants).
                              Exact - counting is len() over the rows.
    vep_warnings              the `_warnings.txt` blob captured on the AnnotationRun. All a run
                              annotated before that flag was turned on has.

VEP >= 112 routes every skipped record through BaseVEP::skipped_variant_msg, which formats it as

    WARNING: line <N> skipped (<input line, whitespace-squeezed, cut to 50 chars>): <description>

where N is the Nth *data* record of the input VCF (the parser starts its line_number at 0 and advances
once per record in create_VariationFeatures; header lines are consumed by htslib and never counted).
The pre-112 wording is recognised too, because a scan of historical rows meets blobs written by
whichever VEP was current at the time.

Nothing here touches the ORM, so a migration can import it.
"""
import logging
import re
from dataclasses import dataclass
from typing import Optional

# Every skip goes through skipped_variant_msg, so matching its common prefix counts all of them -
# including reasons nobody has seen yet. Matching on descriptions counts only the ones we thought to
# enumerate, which is how the previous parser came to return zero for every VEP 116 run.
SKIPPED_LINE_PATTERN = re.compile(r"^WARNING: line (\d+) skipped\b")
# Wording shared by both eras - the contig name is the one thing handle_vep_skipped needs by name.
CONTIG_NOT_FOUND_PATTERN = re.compile(r"Chromosome (\S+) not found")
# Pre-112 wording. The contig warning carried its own line number rather than a "line N skipped" prefix.
LEGACY_ON_LINE_PATTERN = re.compile(r"\bon line (\d+)")
LEGACY_INCOMPLETE_PATTERN = re.compile(r"WARNING: VCF line on line (\d+) looks incomplete, skipping:")
LEGACY_VARIANT_ID_PATTERN = re.compile(r"variant_id=(\d+)")

SKIPPED_VARIANTS_FILE_HEADER_PREFIX = "[VEP skipped the following variants from "
SKIPPED_VARIANTS_FILE_LINE_PREFIX = "Line "


class TruncatedVEPOutputError(ValueError):
    """ #1701: VEP's annotated VCF holds fewer records than we handed it, and VEP did not report the
        difference as skipped - so records were lost rather than dropped on purpose. Importing it would
        write vep_skipped_reason=UNKNOWN rows that make the missing variants look annotated. """


@dataclass(frozen=True)
class VEPWarnings:
    skipped_line_numbers: set[int]  # dump-VCF data lines VEP dropped
    skipped_contigs: set[str]  # contigs VEP could not find
    incomplete_variant_ids: set[int]  # legacy pre-112 recovery

    @property
    def skipped_count(self) -> int:
        return len(self.skipped_line_numbers)


@dataclass(frozen=True)
class VEPSkippedVariant:
    line_number: Optional[int]
    line: str
    description: str


def parse_vep_warnings(vep_warnings: Optional[str]) -> VEPWarnings:
    """ Pull the skipped records out of a VEP `_warnings.txt` blob. """
    skipped_line_numbers: set[int] = set()
    skipped_contigs: set[str] = set()
    incomplete_variant_ids: set[int] = set()

    awaiting_incomplete_variant_id = False
    for line in (vep_warnings or "").splitlines():
        if m := SKIPPED_LINE_PATTERN.search(line):
            skipped_line_numbers.add(int(m.group(1)))
        elif m := LEGACY_INCOMPLETE_PATTERN.search(line):
            # The variant_id lives on the following line, which is the offending input line echoed back
            skipped_line_numbers.add(int(m.group(1)))
            awaiting_incomplete_variant_id = True
            continue
        elif awaiting_incomplete_variant_id:
            if m := LEGACY_VARIANT_ID_PATTERN.search(line):
                incomplete_variant_ids.add(int(m.group(1)))
        awaiting_incomplete_variant_id = False

        if m := CONTIG_NOT_FOUND_PATTERN.search(line):
            skipped_contigs.add(m.group(1))
            if m := LEGACY_ON_LINE_PATTERN.search(line):
                skipped_line_numbers.add(int(m.group(1)))

    return VEPWarnings(skipped_line_numbers=skipped_line_numbers,
                       skipped_contigs=skipped_contigs,
                       incomplete_variant_ids=incomplete_variant_ids)


def parse_skipped_variants_file(filename: str) -> list[VEPSkippedVariant]:
    r""" Read VEP's --skipped_variants_file, written as

            [VEP skipped the following variants from <input>]
            Line <n>\t<input line, cut to 40 chars>\t<description>

        (printf "Line %-5s\t%-40.40s\t%s\n", so both leading fields are blank-padded). """
    skipped_variants = []
    with open(filename) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith(SKIPPED_VARIANTS_FILE_HEADER_PREFIX):
                continue
            fields = line.split("\t")
            if len(fields) < 3 or not fields[0].startswith(SKIPPED_VARIANTS_FILE_LINE_PREFIX):
                logging.warning("Unexpected row in VEP skipped variants file '%s': %r", filename, line)
                continue
            line_number_str = fields[0][len(SKIPPED_VARIANTS_FILE_LINE_PREFIX):].strip()
            skipped_variants.append(VEPSkippedVariant(
                line_number=int(line_number_str) if line_number_str.isdigit() else None,
                line=fields[1].strip(),
                description="\t".join(fields[2:]).strip(),
            ))
    return skipped_variants


def annotation_run_shortfall(dump_count: int, annotated_count: int, vep_warnings: Optional[str]) -> int:
    """ #1701: records handed to VEP that neither came back annotated nor were reported skipped.

        A plain function over row values (no ORM) so the Part C command and its ManualOperation
        migration run identical arithmetic. """
    return dump_count - (annotated_count + parse_vep_warnings(vep_warnings).skipped_count)
