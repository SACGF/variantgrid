from library.genomics.vcf_utils import get_contigs_header_lines


def get_vcf_header_lines(top_lines=None, info_dict=None, formats=None, contig_lines=None, samples=None):
    """ info_dict - values = dict with keys of 'type', 'description' (optional 'number' default = 1) """

    if samples is None:
        samples = []

    header_lines = ['##fileformat=VCFv4.1']
    if top_lines:
        header_lines.extend(top_lines)
    if info_dict:
        for info_id, data in info_dict.items():
            data['id'] = info_id
            if data.get("number") is None:
                data["number"] = 1
            line_template = '##INFO=<ID=%(id)s,Number=%(number)s,Type=%(type)s,Description="%(description)s">'
            line = line_template % data
            header_lines.append(line)

    use_format = samples and formats
    if use_format:
        header_lines.extend(formats)

    if contig_lines:
        header_lines.extend(contig_lines)

    colnames = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    if use_format:
        colnames += ['FORMAT'] + samples

    header_lines.append('#' + '\t'.join(colnames))

    return header_lines


def get_vcf_header_from_contigs(genome_build, info_dict=None, samples=None, use_accession=True):
    """ info_dict which contains ('number', 'type', 'description') """

    formats = ['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
               '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
               '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
               '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">',
               '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">']

    contig_lines = get_contigs_header_lines(genome_build, use_accession=use_accession)
    return get_vcf_header_lines(info_dict=info_dict, formats=formats, contig_lines=contig_lines, samples=samples)
