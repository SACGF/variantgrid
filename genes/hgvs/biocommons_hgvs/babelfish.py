# This is from BioCommons, with a few fixes applied:
# https://github.com/biocommons/hgvs/blob/main/src/hgvs/extras/babelfish.py
# Will remove this once this issue / pull request is closed
# https://github.com/biocommons/hgvs/issues/653

"""translate between HGVS and other formats"""
import os

import hgvs
import hgvs.normalizer
from bioutils.assemblies import make_ac_name_map, make_name_ac_map
from bioutils.sequences import reverse_complement
from hgvs.assemblymapper import AssemblyMapper
from hgvs.edit import NARefAlt
from hgvs.exceptions import HGVSInvalidVariantError
from hgvs.location import Interval, SimplePosition
from hgvs.normalizer import Normalizer
from hgvs.parser import Parser
from hgvs.posedit import PosEdit
from hgvs.sequencevariant import SequenceVariant
from hgvs.validator import ExtrinsicValidator


def _as_interbase(posedit):
    if posedit.edit.type == "ins":
        # ins coordinates (only) exclude left position
        start_i = posedit.pos.start.base
        end_i = posedit.pos.end.base - 1
    else:
        start_i = posedit.pos.start.base - 1
        end_i = posedit.pos.end.base
    return (start_i, end_i)


class ParserSingleton:
    __instance = None

    def __init__(self):
        self._parser = Parser()

    @classmethod
    def parser(cls):
        if not cls.__instance:
            cls.__instance = ParserSingleton()
        return cls.__instance._parser


class Babelfish:
    """ The plan here is to get this into biocommons HGVS via pull request
        Once this is done, we can delete this file and use the biocommons library """
    def __init__(self, hdp, assembly_name):
        self.assembly_name = assembly_name
        self.hdp = hdp
        self.hn = hgvs.normalizer.Normalizer(hdp, cross_boundaries=False, shuffle_direction=5, validate=False)
        self.am = AssemblyMapper(self.hdp,
                                 assembly_name=assembly_name,
                                 alt_aln_method='splign', replace_reference=True)
        self.ev = ExtrinsicValidator(self.hdp)
        self.ac_to_name_map = make_ac_name_map(assembly_name)
        self.name_to_ac_map = make_name_ac_map(assembly_name)

    def hgvs_to_vcf(self, var_g):
        """**EXPERIMENTAL**

        converts a single hgvs allele to (chr, pos, ref, alt) using
        the given assembly_name. The chr name uses official chromosome
        name (i.e., without a "chr" prefix).

        Returns None for non-variation (e.g., NC_000006.12:g.49949407=)

        """

        if var_g.type != "g":
            raise RuntimeError("Expected g. variant, got {var_g}".format(var_g=var_g))

        vleft = self.hn.normalize(var_g)

        (start_i, end_i) = _as_interbase(vleft.posedit)

        chr = self.ac_to_name_map[vleft.ac]

        typ = vleft.posedit.edit.type

        if typ == "dup":
            start_i -= 1
            alt = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
            ref = alt[0]
            end_i = start_i
            return (chr, start_i + 1, ref, alt, typ)
        elif typ == 'inv':
            ref = vleft.posedit.edit.ref
            alt = reverse_complement(ref)
        else:
            if vleft.posedit.edit.ref == vleft.posedit.edit.alt:
                return None

            alt = vleft.posedit.edit.alt or ""

            if typ in ('del', 'ins'):  # Left anchored
                # When making biocommons pull request, add comment:
                # The original test here was: end_i - start_i == 1 and vleft.posedit.length_change() == 0
                # The only difference is the old way caught delins and gave them a common prefix
                # example NC_000003.12:g.128483940_128483945delinsC
                # old way: 3:128202782 TGGCCGG>TC
                # new way: 3:128202783 GGCCGG>C (this is what VT normalizes the above to)

                start_i -= 1
                ref = self.hdp.seqfetcher.fetch_seq(vleft.ac, start_i, end_i)
                alt = ref[0] + alt
            else:
                ref = vleft.posedit.edit.ref

        return chr, start_i + 1, ref, alt, typ

    def vcf_to_hgvs(self, chrom, position, ref, alt, transcript_accession=None):
        if transcript_accession:
            return self.vcf_to_c_hgvs(chrom, position, ref, alt, transcript_accession)
        return self.vcf_to_g_hgvs(chrom, position, ref, alt)

    def vcf_to_c_hgvs(self, chrom, position, ref, alt, transcript_accession):
        var_g = self.vcf_to_g_hgvs(chrom, position, ref, alt)
        return self.am.g_to_c(var_g, transcript_accession)

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

    def hgvs_to_g_hgvs(self, hgvs_string: str):
        CONVERT_TO_G = {
            'c': self.am.c_to_g,
            'n': self.am.n_to_g,
        }

        parser = ParserSingleton.parser()
        var_x = parser.parse_hgvs_variant(hgvs_string)
        try:
            self.ev.validate(var_x, strict=True)  # Validate in transcript range
        except HGVSInvalidVariantError as hgvs_e:
            ACCEPTABLE_VALIDATION_MESSAGES = [
                'Cannot validate sequence of an intronic variant',
                'does not agree with reference sequence'
            ]
            ok = False
            exception_str = str(hgvs_e)
            for msg in ACCEPTABLE_VALIDATION_MESSAGES:
                if msg in exception_str:
                    ok = True
                    break
            if not ok:
                raise

        if converter := CONVERT_TO_G.get(var_x.type):
            var_x = converter(var_x)
        return var_x


if __name__ == "__main__":
    """
      49949___  400       410       420
                  |123456789|123456789|
    NC_000006.12  GACCAGAAAGAAAAATAAAAC

    """

    import hgvs.easy
    import hgvs.normalizer
    from hgvs.extras.babelfish import Babelfish

    babelfish38 = Babelfish(hgvs.easy.hdp, assembly_name="GRCh38")
    hnl = hgvs.normalizer.Normalizer(hgvs.easy.hdp, cross_boundaries=False, shuffle_direction=5, validate=False)

    def _h2v(h):
        return babelfish38.hgvs_to_vcf(hgvs.easy.parser.parse(h))

    def _v22(*v):
        return babelfish38.vcf_to_hgvs(*v)

    def _vp(h):
        v = hgvs.easy.parser.parse(h)
        vl = hnl.normalize(v)
        return (v, vl)

    for h in (
        # Non-variation
        "NC_000006.12:g.49949407=",
        # SNV
        "NC_000006.12:g.49949407A>T",
        # delins
        "NC_000006.12:g.49949413_49949414delinsCC",
        # del
        "NC_000006.12:g.49949415del",
        "NC_000006.12:g.49949413del",
        "NC_000006.12:g.49949414del",
        "NC_000006.12:g.49949413_49949414del",
        # ins
        "NC_000006.12:g.49949413_49949414insC",
        "NC_000006.12:g.49949414_49949415insC",
        "NC_000006.12:g.49949414_49949415insCC",
        # ins (dup)
        "NC_000006.12:g.49949413_49949414insA",
        "NC_000006.12:g.49949414_49949415insA",
        "NC_000006.12:g.49949414_49949415insAA",
    ):
        v = _h2v(h)
        _v22(*v)  # Try converting it back
        print('assert _h2v("{h}") == {res}'.format(res=str(v), h=h))
