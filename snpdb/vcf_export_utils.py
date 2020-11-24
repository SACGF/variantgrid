
def get_vcf_header_lines(top_lines=None, info_dict=None, formats=None, contigs=None, samples=None):
    """ info_dict - values of ('number', 'type', 'description')
        contigs - tuples of (contig, length, assembly)
    """

    if samples is None:
        samples = []

    header_lines = ['##fileformat=VCFv4.1']
    if top_lines:
        header_lines.extend(top_lines)
    if info_dict:
        for (info_id, data) in info_dict.items():
            data['id'] = info_id
            line_template = '##INFO=<ID=%(id)s,Number=%(number)s,Type=%(type)s,Description="%(description)s">'
            line = line_template % data
            header_lines.append(line)

    use_format = samples and formats
    if use_format:
        header_lines.extend(formats)

    if contigs:
        for (contig, length, assembly) in contigs:
            line = f"##contig=<ID={contig},length={length},assembly={assembly}>"
            header_lines.append(line)

    colnames = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    if use_format:
        colnames += ['FORMAT'] + samples

    header_lines.append('#' + '\t'.join(colnames))

    return header_lines


def get_vcf_header_from_contigs(genome_build, info_dict=None, samples=None):
    """ info_dict which contains ('number', 'type', 'description') """

    contigs = []
    for contig in genome_build.contigs:
        contigs.append((contig.name, contig.length, genome_build.name))

    formats = ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
               '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
               '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
               '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">',
               '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">']

    return get_vcf_header_lines(info_dict=info_dict, formats=formats, contigs=contigs, samples=samples)

