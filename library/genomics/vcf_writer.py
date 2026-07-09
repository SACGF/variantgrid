"""
Shared low-level VCF writing:

  * ``build_header_lines`` (+ the header dataclasses) renders VCF header lines.
  * ``VCFWriter`` writes VCF records as text to any handle with a ``.write(str)`` method
    (a text file, a ``StashFile`` streaming buffer, ...). Callers writing to a binary
    destination (a bgzip stream, a raw pipe) wrap it in ``io.TextIOWrapper`` first, so the
    text/bytes concern stays at the handle, not in the writer.

Parsing / header-rewriting of existing VCFs lives in ``library.genomics.vcf_utils``.
"""
from dataclasses import dataclass
from typing import Iterable, Optional, Union

DEFAULT_FILE_FORMAT = "VCFv4.1"

VCF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


@dataclass
class VCFInfoHeader:
    """ A single ``##INFO=<…>`` header line. """
    id: str
    type: str
    description: str = ""
    number: Union[int, str] = 1

    def __str__(self) -> str:
        # Description is rendered inside double quotes, so any double quotes in the text
        # must become single quotes (kept in one place - see module docstring).
        description = (self.description or "").replace('"', "'")
        return f'##INFO=<ID={self.id},Number={self.number},Type={self.type},Description="{description}">'


def build_header_lines(*, meta_lines: Iterable[str] = None,
                       info: Iterable[Union[VCFInfoHeader, str]] = None,
                       formats: Iterable[str] = None,
                       contig_lines: Iterable[str] = None,
                       samples: Iterable[str] = None,
                       file_format: str = DEFAULT_FILE_FORMAT) -> list[str]:
    """ The single VCF header builder.

        :param meta_lines: raw ``##…`` lines placed straight after ``##fileformat`` (e.g. ``##source=VariantGrid``)
        :param info: ``VCFInfoHeader`` objects (or already-rendered ``##INFO`` strings)
        :param formats: already-rendered ``##FORMAT`` lines - only emitted when ``samples`` are also given
        :param contig_lines: already-rendered ``##contig`` / ``##reference`` lines
        :param samples: sample names - when given together with ``formats`` a ``FORMAT`` column + sample columns are added
        :return: list of header lines (no trailing newlines)
    """
    header_lines = [f"##fileformat={file_format}"]
    if meta_lines:
        header_lines.extend(meta_lines)
    if info:
        header_lines.extend(str(info_header) for info_header in info)

    use_format = bool(samples and formats)
    if use_format:
        header_lines.extend(formats)

    if contig_lines:
        header_lines.extend(contig_lines)

    colnames = list(VCF_COLUMNS)
    if use_format:
        colnames += ["FORMAT"] + list(samples)

    header_lines.append("#" + "\t".join(colnames))
    return header_lines


def symbolic_alt_info(alt: str, *, svlen, end) -> dict:
    """ Single source of truth for the symbolic-alt INFO fields.

        :return: ordered ``{"END", "SVLEN", "SVTYPE"}`` dict for a symbolic ``alt`` (e.g. ``<DEL>``),
                 else an empty dict. Callers may re-order the keys to preserve historical output. """
    if not (alt.startswith("<") and alt.endswith(">")):
        return {}
    return {
        "END": end,
        "SVLEN": svlen,
        "SVTYPE": alt[1:-1],  # Strip off the < > brackets
    }


class VCFWriter:
    """ Writes VCF records (and optionally the header) as text to a handle.

        :param handle: any object with a ``.write(str)`` method (text file, ``StashFile``, ...).
            For a binary destination (bgzip / raw pipe) wrap it in ``io.TextIOWrapper`` first.
        :param header_lines: header lines to write immediately (optional)
        :param encode_info: optional ``value -> value`` hook applied to each INFO value before joining
    """

    def __init__(self, handle, header_lines: Iterable[str] = None, *, encode_info=None):
        self.handle = handle
        self.encode_info = encode_info
        if header_lines:
            for line in header_lines:
                self._write_line(line)

    def _write_line(self, line: str) -> None:
        self.handle.write(line + "\n")

    def _format_info(self, info: Optional[dict]) -> str:
        if not info:
            return "."
        parts = []
        for k, v in info.items():
            if self.encode_info is not None:
                v = self.encode_info(v)
            parts.append(f"{k}={v}")
        return ";".join(parts)

    def write_record(self, chrom, pos, ref, alt, *, vcf_id=".", qual=".", vcf_filter=".",
                     info: Optional[dict] = None, fmt: Optional[str] = None,
                     sample_calls: Optional[Iterable[str]] = None) -> None:
        columns = [str(chrom), str(pos), str(vcf_id), str(ref), str(alt),
                   str(qual), str(vcf_filter), self._format_info(info)]
        if fmt is not None:
            columns.append(fmt)
            if sample_calls:
                columns.extend(str(s) for s in sample_calls)
        self._write_line("\t".join(columns))
