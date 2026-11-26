"""
A ploidy caller's whole-sample results, which it writes into its VCF header rather than any record: DRAGEN
ploidy's ##estimatedSexKaryotype, ##referenceSexKaryotype and ##autosomeDepthOfCoverage. Entry point:
read_ploidy_header, on the raw header text (VCF.header), so import and a backfill read it the same way.
"""
import re
from dataclasses import dataclass
from typing import Optional

ESTIMATED_SEX_KARYOTYPE = "estimatedSexKaryotype"
REFERENCE_SEX_KARYOTYPE = "referenceSexKaryotype"
AUTOSOME_DEPTH_OF_COVERAGE = "autosomeDepthOfCoverage"
_HEADER_LINE_PATTERN = re.compile(r"^##(\w+)=(.*)$", re.MULTILINE)


@dataclass(frozen=True)
class PloidyHeader:
    sex_karyotype: str
    reference_sex_karyotype: Optional[str]
    autosome_depth_of_coverage: Optional[float]


def _to_float(value: Optional[str]) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_ploidy_header(header: Optional[str]) -> Optional[PloidyHeader]:
    """ None unless the header has an estimated sex karyotype - the other two only qualify it """
    values = dict(_HEADER_LINE_PATTERN.findall(header or ""))
    if not (sex_karyotype := values.get(ESTIMATED_SEX_KARYOTYPE, "").strip()):
        return None
    return PloidyHeader(sex_karyotype=sex_karyotype,
                        reference_sex_karyotype=values.get(REFERENCE_SEX_KARYOTYPE, "").strip() or None,
                        autosome_depth_of_coverage=_to_float(values.get(AUTOSOME_DEPTH_OF_COVERAGE)))
