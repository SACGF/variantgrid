"""
Normalize a ranged del/dup/inv of any length by reading only the bases around its breakpoints (#2103).
biocommons hgvs normalization builds the explicit sequence of the whole span, which is quadratic in its
length - a 100kb dup takes 0.5s and a whole chromosome arm never finishes.

TODO: Remove this module once biocommons hgvs Normalizer can shuffle/trim a sequence-less del/dup/inv
without building the span - then normalize the parsed g. with it instead (#2103).

Entry points: shuffle_interval (del/dup), inversion_palindromic_trim (inv). Both take contig_sequence,
anything sliced 0-based half-open as contig_sequence[start:end] that returns short past the contig end
(eg FastaFileContigWrapper).
"""
import os

from bioutils.sequences import reverse_complement

# Bases read per fetch when shuffling a long SV's breakpoints along a repeat
SHUFFLE_WINDOW_SIZE = 10_000


def shuffle_interval(contig_sequence, start: int, end: int, shift_right: bool) -> tuple[int, int]:
    """ Slide a deleted/duplicated 1-based inclusive interval along the reference while the base leaving one
        end equals the base entering the other - ie as far as normalization would move it """
    shift = 0
    while True:
        if shift_right:
            # bases from start leave as bases from end+1 enter - entering comes back short at the contig end
            n = SHUFFLE_WINDOW_SIZE
            leaving = contig_sequence[start - 1 + shift:start - 1 + shift + n]
            entering = contig_sequence[end + shift:end + shift + n]
        else:
            # bases from end leave as bases from start-1 enter, read outwards from the breakpoints
            n = min(SHUFFLE_WINDOW_SIZE, start - 1 - shift)
            leaving = contig_sequence[end - shift - n:end - shift][::-1]
            entering = contig_sequence[start - 1 - shift - n:start - 1 - shift][::-1]
        matched = len(os.path.commonprefix([leaving.upper(), entering.upper()]))
        shift += matched
        if n <= 0 or matched < n:
            break
    if not shift_right:
        shift = -shift
    return start + shift, end + shift


def inversion_palindromic_trim(contig_sequence, start: int, end: int) -> int:
    """ Bases normalization trims off each end of a 1-based inclusive inversion - while the first base
        is the complement of the last, inverting them changes nothing """
    trim = 0
    while True:
        n = min(SHUFFLE_WINDOW_SIZE, (end - start + 1) // 2 - trim)
        leading = contig_sequence[start - 1 + trim:start - 1 + trim + n]
        trailing = contig_sequence[end - trim - n:end - trim]
        matched = len(os.path.commonprefix([leading.upper(), reverse_complement(trailing).upper()]))
        trim += matched
        if n <= 0 or matched < n:
            return trim
