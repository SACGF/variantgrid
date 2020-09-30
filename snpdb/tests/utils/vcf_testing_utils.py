import vcf

from snpdb.models import Locus, Variant, Sequence


def slowly_create_loci_and_variants_for_vcf(genome_build, vcf_filename, get_variant_id_from_info=False):
    """ This is very slow - only use for tests """

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
