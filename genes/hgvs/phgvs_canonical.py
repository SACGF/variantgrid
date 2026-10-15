"""
The one canonical p.HGVS string, whatever produced it (VEP, our own c_to_p, ClinGen or a lab's typing), and its
short display form. Canonical is p-dot only (no protein accession), 3-letter, Ter for stop and always in predicted
brackets - p.(Gly12Asp), p.(Met246IlefsTer14) - except p.? and p.(=).

Entry points: canonical_p_hgvs (any string), canonical_from_sequence_variant (a biocommons p. variant),
short_p_hgvs (canonical -> p.(G12D)) and is_frameshift_in_disguise. Parsing is biocommons' - genes/hgvs/phgvs.py:PHGVS
is a regex for single residue changes and can't read ranges, dups or delins (#1383).
"""
import re
from typing import Optional
from urllib.parse import unquote

from hgvs.parsers import Parser
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant

P_HGVS_UNKNOWN = "p.?"
P_HGVS_NO_CHANGE = "p.(=)"

_CANONICAL_CONF = {"p_3_letter": True, "p_term_asterisk": False}
_SHORT_CONF = {"p_3_letter": False, "p_term_asterisk": True}

# The p-dot starts the string, or follows the protein accession's ':' (or a gene in brackets / whitespace)
_P_DOT = re.compile(r"(?:^|[:)\s])p\.(.*)$")
# Legacy 'X' for a stop - biocommons reads a bare X as Xaa (unknown residue), so p.Trp26X would be p.(Trp26Xaa)
_LEGACY_STOP_X = re.compile(r"X(?!aa)")

_parser = Parser()


def _format_posedit(posedit: Optional[PosEdit], conf: dict) -> str:
    if posedit is None:
        return P_HGVS_UNKNOWN
    if posedit.pos is None and str(posedit.edit) == "=":
        return P_HGVS_NO_CHANGE
    posedit.uncertain = True  # predicted, matching EvidenceMixin.p_parts' override_is_confirmed_to=False
    return f"p.{posedit.format(conf)}"


def _parse_p_dot(raw: str) -> Optional[PosEdit | str]:
    """ PosEdit, or P_HGVS_UNKNOWN for p.?, or None if it doesn't parse """
    if not raw:
        return None
    if not (m := _P_DOT.search(unquote(raw.strip()))):  # VEP writes p.%3D for p.=
        return None
    p_dot = re.sub(r"\s+", "", m.group(1))
    if p_dot in ("?", "(?)"):
        return P_HGVS_UNKNOWN
    p_dot = _LEGACY_STOP_X.sub("Ter", p_dot)
    try:
        return _parser.parse_p_posedit(p_dot)
    except Exception:  # Not just HGVSParseError - a residue the grammar knows but can't map ('Glx') is a KeyError
        return None


def canonical_p_hgvs(raw: str) -> Optional[str]:
    """ The canonical form of any p.HGVS string, or None if it isn't one biocommons can parse """
    posedit = _parse_p_dot(raw)
    if posedit is None or isinstance(posedit, str):
        return posedit
    return _format_posedit(posedit, _CANONICAL_CONF)


def short_p_hgvs(canonical: Optional[str]) -> Optional[str]:
    """ 1-letter, '*' for stop: p.(Met246IlefsTer14) -> p.(M246Ifs*14). Returns the input if it doesn't parse """
    posedit = _parse_p_dot(canonical)
    if posedit is None:
        return canonical
    if isinstance(posedit, str):
        return posedit
    return _format_posedit(posedit, _SHORT_CONF)


def canonical_from_sequence_variant(variant: SequenceVariant) -> str:
    """ Canonical form of a biocommons p. variant (eg AssemblyMapper.c_to_p output) - drops the protein accession """
    return _format_posedit(variant.posedit, _CANONICAL_CONF)


def is_frameshift_in_disguise(variant: SequenceVariant, protein_length: int) -> bool:
    """ biocommons can't name some frameshifts (eg from a dup) and writes a delins running to the end of the protein,
        p.(Ser1321_Val2843delinsIleLeuTer) where VEP writes p.Ser1321IlefsTer3.
        protein_length is the residue count excluding the stop """
    posedit = variant.posedit
    if posedit is None or posedit.pos is None:
        return False
    edit_type = getattr(posedit.edit, "type", None)
    return edit_type == "delins" and posedit.pos.end.pos >= protein_length
