def format_chrom(chrom, want_chr):
    """ Pass in a chromosome (unknown format), return in your format (CONTIG UPPERCASED)
        @param chrom: chromosome ID (eg 1 or 'chr1')
        @param want_chr: Boolean - whether you want "CHR" at the beginning of chrom
        @return: "chr1" or "1" (for want_chr True/False)
    """

    chrom_no_chr = chrom.upper().replace("CHR", "")
    if want_chr:
        return f"chr{chrom_no_chr}"
    return chrom_no_chr


def get_genomic_size_description(genomic_size):
    if genomic_size == 1000000:
        genomic_size_description = "mb"
    elif genomic_size == 1000:
        genomic_size_description = "kb"
    else:
        genomic_size_description = f"{genomic_size} bases"

    return genomic_size_description
