"""
The LibraryQC rows of a pair's arm as one line, for the pages that list QC beside something else - the
specimen page and the sequencing run's Data tab. The metrics themselves live on the pair's page
(seqauto.views.view_tso500_pair), which every summary links to.

Entry point is summarise_library_qc, rendered by seqauto/templates/seqauto/library_qc_table.html.
"""
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Optional

from patients.models import Specimen
from patients.models_enums import NucleicAcid
from seqauto.models import LibraryQC, SequencingRun


@dataclass
class LibraryQCArmSummary:
    """ Every QC category of one arm of one pair on one run - the five DNA categories or the RNA one """
    sequencing_run_name: str
    sequencing_run: Optional[SequencingRun]
    pair_id: str
    nucleic_acid: str  # NucleicAcid's label - 'DNA' / 'RNA'
    specimen: Optional[Specimen]
    rows: list[LibraryQC]

    @property
    def url(self) -> str:
        return self.rows[0].get_absolute_url()

    @property
    def passed(self) -> Optional[bool]:
        """ The arm as a whole - False where any category failed, None where none of them judged """
        judged = [row.passed for row in self.rows if row.passed is not None]
        if not judged:
            return None
        return all(judged)

    @property
    def completed(self) -> Optional[bool]:
        """ [Analysis Status] for the pair, which every row carries """
        return self.rows[0].completed

    @property
    def method(self) -> str:
        return self.rows[0].method

    @property
    def measured_date(self):
        return self.rows[0].measured_date


def summarise_library_qc(library_qc: Iterable[LibraryQC]) -> list[LibraryQCArmSummary]:
    """ One summary per (run, pair, arm), in the order the rows come in """
    summaries: dict[tuple[str, str, str], LibraryQCArmSummary] = {}
    for row in library_qc:
        key = (row.sequencing_run_name, row.pair_id, row.nucleic_acid)
        summary = summaries.get(key)
        if summary is None:
            summary = LibraryQCArmSummary(sequencing_run_name=row.sequencing_run_name,
                                          sequencing_run=row.sequencing_run,
                                          pair_id=row.pair_id,
                                          nucleic_acid=NucleicAcid(row.nucleic_acid).label,
                                          specimen=row.specimen,
                                          rows=[])
            summaries[key] = summary
        summary.rows.append(row)
    return list(summaries.values())
