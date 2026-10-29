"""
The one canonical p.HGVS string, whatever produced it (VEP, our own c_to_p, ClinGen or a lab's typing), and its
short display form. Canonical keeps the protein accession the source gave (none if it gave none) and writes the p-dot
3-letter, Ter for stop and always in predicted brackets - NP_000483.3:p.(Gly480Cys), p.(Met246IlefsTer14) - except
p.? and p.(=). A change to the initiation codon is p.(Met1?) whatever the source wrote (ClinGen gives p.Met1Val):
it stops translation starting rather than substituting a residue (#1383).

Entry points: canonical_p_hgvs (any string), canonical_from_sequence_variant (a biocommons p. variant),
short_p_hgvs (canonical -> p.(G12D), no accession) and is_frameshift_in_disguise. Parsing is biocommons', so
ranges, dups and delins read as well as single residue changes (#1383).
"""
import re
from typing import Optional
from urllib.parse import unquote

from hgvs.edit import AARefAlt
from hgvs.parsers import Parser
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant

P_HGVS_UNKNOWN = "p.?"
P_HGVS_NO_CHANGE = "p.(=)"

_CANONICAL_CONF = {"p_3_letter": True, "p_term_asterisk": False}
_SHORT_CONF = {"p_3_letter": False, "p_term_asterisk": True}

_INITIATION_CODON_UNKNOWN = AARefAlt(ref="", alt="", init_met=True)  # formats as Met1? / M1?
_INITIATION_CODON_LOSS_EDIT_TYPES = {"sub", "del", "delins", "fs"}  # ext (a new upstream start) and dup keep Met1

# The p-dot starts the string, or follows the protein accession's ':' (gene in brackets allowed) or whitespace
_P_DOT = re.compile(r"(?:^|(?:(?P<accession>[^\s:()]+)(?:\([^)]*\))?:)|[)\s])p\.(?P<p_dot>.*)$")
# Legacy 'X' for a stop - biocommons reads a bare X as Xaa (unknown residue), so p.Trp26X would be p.(Trp26Xaa)
_LEGACY_STOP_X = re.compile(r"X(?!aa)")

_parser = Parser()


def _is_initiation_codon_change(posedit: PosEdit | AARefAlt) -> bool:
    """ c_to_p gives a bare init_met edit with no position, a parsed string a change starting at Met1 """
    if getattr(posedit, "init_met", False):
        return True
    if not (pos := getattr(posedit, "pos", None)):
        return False
    edit = posedit.edit
    if getattr(edit, "type", None) not in _INITIATION_CODON_LOSS_EDIT_TYPES or edit.ref == getattr(edit, "alt", None):
        return False
    return pos.start.base == 1 and pos.start.aa == "M"


def _format_posedit(posedit: Optional[PosEdit | AARefAlt], conf: dict) -> str:
    if posedit is None:
        return P_HGVS_UNKNOWN
    if _is_initiation_codon_change(posedit):
        return f"p.({_INITIATION_CODON_UNKNOWN.format(conf)})"
    if posedit.pos is None and str(posedit.edit) == "=":
        return P_HGVS_NO_CHANGE
    posedit.uncertain = True  # always predicted: we can't confirm a lab's p.HGVS was experimentally determined
    return f"p.{posedit.format(conf)}"


def _with_accession(accession: Optional[str], p_dot: str) -> str:
    return f"{accession}:{p_dot}" if accession else p_dot


def _parse(raw: str) -> tuple[Optional[str], Optional[PosEdit | str]]:
    """ (accession, PosEdit) - the PosEdit is P_HGVS_UNKNOWN for p.?, or None if it doesn't parse """
    if not raw:
        return None, None
    if not (m := _P_DOT.search(unquote(raw.strip()))):  # VEP writes p.%3D for p.=
        return None, None
    accession = m.group("accession")
    p_dot = re.sub(r"\s+", "", m.group("p_dot"))
    if p_dot in ("?", "(?)"):
        return accession, P_HGVS_UNKNOWN
    p_dot = _LEGACY_STOP_X.sub("Ter", p_dot)
    try:
        return accession, _parser.parse_p_posedit(p_dot)
    except Exception:  # Not just HGVSParseError - a residue the grammar knows but can't map ('Glx') is a KeyError
        return None, None


def canonical_p_hgvs(raw: str) -> Optional[str]:
    """ The canonical form of any p.HGVS string, or None if it isn't one biocommons can parse """
    accession, posedit = _parse(raw)
    if posedit is None:
        return None
    if isinstance(posedit, str):
        return _with_accession(accession, posedit)
    return _with_accession(accession, _format_posedit(posedit, _CANONICAL_CONF))


def short_p_hgvs(canonical: Optional[str]) -> Optional[str]:
    """ 1-letter, '*' for stop, no accession: NP_1.1:p.(Met246IlefsTer14) -> p.(M246Ifs*14).
        Returns the input if it doesn't parse """
    _, posedit = _parse(canonical)
    if posedit is None:
        return canonical
    if isinstance(posedit, str):
        return posedit
    return _format_posedit(posedit, _SHORT_CONF)


def canonical_from_sequence_variant(variant: SequenceVariant) -> str:
    """ Canonical form of a biocommons p. variant (eg AssemblyMapper.c_to_p output) on its protein accession """
    return _with_accession(variant.ac, _format_posedit(variant.posedit, _CANONICAL_CONF))


def is_frameshift_in_disguise(variant: SequenceVariant, protein_length: int) -> bool:
    """ biocommons can't name some frameshifts (eg from a dup) and writes a delins running to the end of the protein,
        p.(Ser1321_Val2843delinsIleLeuTer) where VEP writes p.Ser1321IlefsTer3.
        protein_length is the residue count excluding the stop """
    posedit = variant.posedit
    if getattr(posedit, "pos", None) is None:  # p.? or an init_met edit
        return False
    edit_type = getattr(posedit.edit, "type", None)
    return edit_type == "delins" and posedit.pos.end.pos >= protein_length
