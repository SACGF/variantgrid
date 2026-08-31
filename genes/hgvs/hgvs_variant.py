from typing import Optional

from bioutils.sequences import reverse_complement
from hgvs.edit import NARefAlt
from hgvs.sequencevariant import SequenceVariant


def _looks_like_transcript(accession: str) -> bool:
    return accession.startswith(('NM_', 'NR_', 'XM_', 'XR_', 'ENST'))


class HGVSVariant:
    """ Wraps a BioCommons SequenceVariant, providing the accessors our code uses """

    def __init__(self, sequence_variant: SequenceVariant):
        self._sequence_variant = sequence_variant

    @property
    def contig_accession(self) -> str:
        _genomic_kinds = ('g', 'm')
        if self.kind not in _genomic_kinds:
            raise ValueError(f"'{self}' can only request contig for genomic kinds '{','.join(_genomic_kinds)}'")
        return self._sequence_variant.ac

    @property
    def gene(self) -> str:
        return self._sequence_variant.gene

    @gene.setter
    def gene(self, value):
        self._sequence_variant.gene = value

    @property
    def length(self) -> int:
        return self._sequence_variant.posedit.length_change()

    @property
    def transcript(self) -> str:
        return self._sequence_variant.ac

    @transcript.setter
    def transcript(self, value):
        self._sequence_variant.ac = value

    @property
    def kind(self) -> str:
        return self._sequence_variant.type

    @kind.setter
    def kind(self, value):
        self._sequence_variant.type = value

    @property
    def ref_allele(self):
        return self._sequence_variant.posedit.edit.ref

    @ref_allele.setter
    def ref_allele(self, value):
        self._sequence_variant.posedit.edit.ref = value

    @property
    def mutation_type(self) -> str:
        biocommons_type = self._sequence_variant.posedit.edit.type
        if biocommons_type == "sub":
            mutation_type = ">"
        else:
            mutation_type = biocommons_type
        return mutation_type

    def get_ref_alt(self):
        edit = self._sequence_variant.posedit.edit
        ref = edit.ref or ''
        alt = edit.alt or ''
        return ref, alt

    def get_cdna_coords(self) -> str:
        return str(self._sequence_variant.posedit.pos.start)

    def format(self, use_delins_for_inv: bool = False, max_ref_length=None):
        sv: SequenceVariant = self._sequence_variant
        if use_delins_for_inv:
            if sv.posedit.edit.type == "inv":
                ref = sv.posedit.edit.ref
                sv.posedit.edit = NARefAlt(ref=ref, alt=reverse_complement(ref))

        conf = {}
        if max_ref_length:
            conf["max_ref_length"] = max_ref_length
        return sv.format(conf)

    def get_gene_symbol_if_no_transcript(self) -> Optional[str]:
        # Biocommons SequenceVariant works like:
        # NM_001145661.2:c.1113dup          transcript=NM_001145661.2, gene=None
        # NM_001145661.2(GATA2):c.1113dup   transcript=NM_001145661.2, gene=GATA2
        # GATA2:c.1113dup                   transcript=GATA2
        gene_symbol = None
        if not (self.gene and self.transcript):
            # Will always be transcript
            if not _looks_like_transcript(self.transcript):
                gene_symbol = self.transcript
        return gene_symbol

    def __repr__(self):
        return self.format()

    def __eq__(self, other):
        if isinstance(other, HGVSVariant):
            return self.format() == other.format()
        return False
