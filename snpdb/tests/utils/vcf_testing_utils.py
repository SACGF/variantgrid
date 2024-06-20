import vcf

from library.utils import sha256sum_str
from snpdb.models import Locus, Variant, Sequence, GenomeBuild, Allele, VariantAllele, AlleleOrigin, \
    AlleleConversionTool, VariantCoordinate


def slowly_create_test_variant(chrom: str, position: int, ref: str, alt: str, genome_build: GenomeBuild) -> Variant:
    """ For test only - doesn't use VariantPKLookup """
    vc = VariantCoordinate.from_explicit_no_svlen(chrom, position, ref, alt)
    vc = vc.as_internal_symbolic(genome_build)
    contig = genome_build.contigs.get(name=vc.chrom)
    uref = vc.ref.upper()
    ualt = vc.alt.upper()
    ref_seq, _ = Sequence.objects.get_or_create(seq=uref, seq_sha256_hash=sha256sum_str(uref))
    alt_seq, _ = Sequence.objects.get_or_create(seq=ualt, seq_sha256_hash=sha256sum_str(ualt))
    locus, _ = Locus.objects.get_or_create(contig=contig, position=position, ref=ref_seq)
    defaults = {"end": vc.end}
    variant, _ = Variant.objects.get_or_create(locus=locus, alt=alt_seq, svlen=vc.svlen, defaults=defaults)
    return variant


def create_mock_allele(variant: Variant, genome_build: GenomeBuild):
    allele = Allele.objects.create()
    VariantAllele.objects.create(
        variant=variant,
        genome_build=genome_build,
        allele=allele,
        origin=AlleleOrigin.IMPORTED_TO_DATABASE,
        allele_linking_tool=AlleleConversionTool.DBSNP
    )
    return allele


def slowly_create_loci_and_variants_for_vcf(genome_build, vcf_filename, get_variant_id_from_info=False):
    """ For tests - doesn't use VariantPKLookup """

    pk_by_seq = Sequence.get_pk_by_seq()
    for v in vcf.Reader(filename=vcf_filename):
        ref = str(v.REF)
        alt = str(v.ALT[0])
        if svlen := v.INFO.get("SVLEN"):
            vc = VariantCoordinate(chrom=v.CHROM, position=int(v.POS), ref=ref, alt=alt, svlen=svlen)
        else:
            vc = VariantCoordinate.from_explicit_no_svlen(v.CHROM, int(v.POS), ref, alt)
        vc = vc.as_internal_symbolic(genome_build)
        ref_id = pk_by_seq.get(ref)
        if ref_id is None:
            sequence = Sequence.objects.create(seq=ref, seq_sha256_hash=sha256sum_str(ref))
            ref_id = sequence.pk
            pk_by_seq[ref] = ref_id

        contig = genome_build.chrom_contig_mappings[vc.chrom]
        locus, _ = Locus.objects.get_or_create(contig=contig,
                                               position=vc.position,
                                               ref_id=ref_id)

        alt_id = pk_by_seq.get(alt)
        if alt_id is None:
            sequence = Sequence.objects.create(seq=alt, seq_sha256_hash=sha256sum_str(alt))
            alt_id = sequence.pk
            pk_by_seq[alt] = alt_id

        kwargs = {"locus": locus,
                  "end": vc.end,
                  "alt_id": alt_id}
        if get_variant_id_from_info:
            kwargs["id"] = v.INFO.get("variant_id")

        Variant.objects.get_or_create(**kwargs)
