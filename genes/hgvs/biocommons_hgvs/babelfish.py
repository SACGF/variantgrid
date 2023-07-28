"""translate between HGVS and other formats"""
import os

import hgvs
import hgvs.normalizer
from bioutils.assemblies import make_ac_name_map, make_name_ac_map
from bioutils.sequences import reverse_complement
from hgvs.edit import NARefAlt
from hgvs.location import Interval, SimplePosition
from hgvs.normalizer import Normalizer
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant


def _as_interbase(posedit):
    if posedit.edit.type == "ins":
        # ins coordinates (only) exclude left position
        start_i = posedit.pos.start.base
        end_i = posedit.pos.end.base - 1
    else:
        start_i = posedit.pos.start.base - 1
        end_i = posedit.pos.end.base
    return (start_i, end_i)


class Babelfish:
    def __init__(self, hdp, assembly_name, symbolic_size_threshold=None):
        """ symbolic_size_threshold - size of alt where you use symbolic alts eg <DEL> or <DUP>
                                    - you need to call hgvs_to_vcf_start_end
        """
        self.assembly_name = assembly_name
        self.hdp = hdp
        self.hn = hgvs.normalizer.Normalizer(hdp, cross_boundaries=False, shuffle_direction=5, validate=False)
        self.ac_to_name_map = make_ac_name_map(assembly_name)
        self.name_to_ac_map = make_name_ac_map(assembly_name)
        self.symbolic_size_threshold = symbolic_size_threshold

    def hgvs_to_vcf(self, var_g):
        if self.symbolic_size_threshold is not None:
            raise ValueError("You must call call 'hgvs_to_vcf_with_end' or 'hgvs_to_vcf_and_info' "
                             "when 'symbolic_size_threshold' is set")
        chrom, start_i, _, ref, alt, typ = self.hgvs_to_vcf_with_end(var_g)
        return chrom, start_i, ref, alt, typ

    @staticmethod
    def allele_is_symbolic(allele):
        return allele.startswith("<") and allele.endswith(">")

    def hgvs_to_vcf_with_info(self, var_g):
        chrom, start_i, end_i, ref, alt, typ = self.hgvs_to_vcf_with_end(var_g)
        info = {}
        return chrom, start_i, ref, alt, info, typ

    def hgvs_to_vcf_with_end(self, var_g):
        """**EXPERIMENTAL**

        converts a single hgvs allele to (chr, start, end, ref, alt, type) using
        the given assembly_name. The chr name uses official chromosome
        name (i.e., without a "chr" prefix).
        """

        if var_g.type != "g":
            raise RuntimeError("Expected g. variant, got {var_g}".format(var_g=var_g))

        vleft = self.hn.normalize(var_g)

        start_i, end_i = _as_interbase(vleft.posedit)

        chrom = self.ac_to_name_map[vleft.ac]

        typ = vleft.posedit.edit.type

        # Symbolic
        if self.symbolic_size_threshold is not None and typ in ('dup', 'del'):
            length = end_i - start_i
            if length >= self.symbolic_size_threshold:
                # Just need 1 base REF
                ref = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i-1, start_i)
                alt = "<" + typ.upper() + ">"
                return chrom, start_i, end_i, ref, alt, typ

        if typ == "dup":
            start_i -= 1
            alt = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
            ref = alt[0]
        elif typ == 'inv':
            ref = vleft.posedit.edit.ref
            alt = reverse_complement(ref)
        else:
            alt = vleft.posedit.edit.alt or ""

            if typ in ('del', 'ins'):  # Left anchored
                start_i -= 1
                ref = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
                alt = ref[0] + alt
            else:
                ref = vleft.posedit.edit.ref
                if ref == alt:
                    alt = '.'
        return chrom, start_i + 1, end_i, ref, alt, typ

    def vcf_to_g_hgvs(self, chrom, position, ref, alt):
        ac = self.name_to_ac_map[chrom]

        keep_left_anchor = False
        if not keep_left_anchor:
            pfx = os.path.commonprefix([ref, alt])
            lp = len(pfx)
            if lp > 0:
                ref = ref[lp:]
                alt = alt[lp:]
                position += lp

        if ref == '':  # Insert
            # Insert uses coordinates around the insert point.
            start = position - 1
            end = position
        else:
            start = position
            end = position + len(ref) - 1

        if alt == '.':
            alt = ref

        var_g = SequenceVariant(ac=ac,
                                type='g',
                                posedit=PosEdit(Interval(start=SimplePosition(start),
                                                         end=SimplePosition(end),
                                                         uncertain=False),
                                                NARefAlt(ref=ref or None,
                                                         alt=alt or None,
                                                         uncertain=False)))
        n = Normalizer(self.hdp)
        return n.normalize(var_g)
