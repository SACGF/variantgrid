""" Premature termination codon (PTC) maths for frameshift variants - @see #579

    VEP's NMD.pm (and LOFTEE's 50bp rule) evaluate the NMD escape rules against the
    variant's own position. For a frameshift the clinically relevant position is the new
    stop codon, which can be hundreds of bases further along the transcript and fall the
    other side of the 50nt boundary. These are the same four rules, anchored on the PTC.

    Pure functions over primitives - the transcript geometry they're combined with lives in
    annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:TranscriptGeometry
"""
import re
from typing import Optional

from annotation.models.models_enums import NMDEscapeStatus

# HGVS numbers the frameshift stop from the first changed residue, counting it as 1, eg
# "NP_000242.1:p.Cys843ValfsTer49". "fsTer?" means VEP couldn't locate the new stop.
_FRAMESHIFT_TER_PATTERN = re.compile(r"fsTer(\d+)$")
# A frameshift whose changed codon is itself a stop has no "fs", eg "p.Tyr535Ter" (#665)
_TER_PATTERN = re.compile(r"[A-Za-z]{3}\d+Ter$")

# NMD.pm: a PTC within 50nt of the last exon-exon junction escapes decay
NMD_LAST_JUNCTION_RULE_NT = 50
# NMD.pm: a PTC in the first 100 coding bases escapes decay (reinitiation downstream)
NMD_START_PROXIMAL_CDS = 100


def parse_ptc_distance_codons(hgvs_p: str) -> Optional[int]:
    """ Codons from the first changed residue to the new stop, counting that residue as 1.
        None when there's no locatable stop. """
    if not hgvs_p:
        return None
    protein_change = hgvs_p.rsplit(":", maxsplit=1)[-1]
    if m := _FRAMESHIFT_TER_PATTERN.search(protein_change):
        return int(m.group(1))
    if _TER_PATTERN.search(protein_change):
        return 1
    return None


def calculate_ptc_position(protein_position: int, ptc_distance_codons: int, indel_offset: int) -> int:
    """ First base of the PTC in reference CDS coordinates.

        @param indel_offset len(ref) - len(alt): everything downstream of the indel shifts by
        this much between the mutant transcript HGVS describes and the reference the exon
        junctions are defined on. """
    ptc_protein_position = protein_position + ptc_distance_codons - 1
    ptc_cds = 3 * (ptc_protein_position - 1) + 1
    return ptc_cds + indel_offset


def predict_nmd_escape(ptc_cds: int, ptc_last_junction_distance: int, single_exon: bool) -> NMDEscapeStatus:
    """ @param ptc_last_junction_distance nucleotides from the PTC to the last exon-exon
        junction - negative means the PTC lies in the final exon. """
    if single_exon:
        return NMDEscapeStatus.ESCAPING
    if ptc_last_junction_distance <= NMD_LAST_JUNCTION_RULE_NT:
        return NMDEscapeStatus.ESCAPING
    if ptc_cds <= NMD_START_PROXIMAL_CDS:
        return NMDEscapeStatus.ESCAPING
    return NMDEscapeStatus.PREDICTED_NMD
