import vcf

from snpdb.models import Locus, Variant, Sequence, GenomeBuild, Allele, VariantAllele, AlleleOrigin, \
    AlleleConversionTool


def slowly_create_test_variant(chrom: str, position: int, ref: str, alt: str, genome_build: GenomeBuild) -> Variant:
    """ For test only - doesn't use VariantPKLookup """
    contig = genome_build.contigs.get(name=chrom)
    ref_seq, _ = Sequence.objects.get_or_create(seq=ref.upper(), length=1)
    alt_seq, _ = Sequence.objects.get_or_create(seq=alt.upper(), length=1)
    locus, _ = Locus.objects.get_or_create(contig=contig, position=position, ref=ref_seq)
    variant, _ = Variant.objects.get_or_create(locus=locus, alt=alt_seq)
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
        ref_id = pk_by_seq.get(ref)
        if ref_id is None:
            sequence = Sequence.objects.create(seq=ref, length=len(ref))
            ref_id = sequence.pk
            pk_by_seq[ref] = ref_id

        contig = genome_build.chrom_contig_mappings[v.CHROM]
        locus, _ = Locus.objects.get_or_create(contig=contig,
                                               position=int(v.POS),
                                               ref_id=ref_id)

        alt_id = pk_by_seq.get(alt)
        if alt_id is None:
            sequence = Sequence.objects.create(seq=alt, length=len(alt))
            alt_id = sequence.pk
            pk_by_seq[alt] = alt_id

        kwargs = {"locus": locus,
                  "alt_id": alt_id}
        if get_variant_id_from_info:
            kwargs["id"] = v.INFO.get("variant_id")

        Variant.objects.get_or_create(**kwargs)
