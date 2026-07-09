from library.genomics.vcf_utils import get_contigs_header_lines
from library.genomics.vcf_writer import VCFInfoHeader, build_header_lines


def get_vcf_header_lines(top_lines=None, info_dict=None, formats=None, contig_lines=None, samples=None):
    """ info_dict - values = dict with keys of 'type', 'description' (optional 'number' default = 1)

        Thin wrapper over library.genomics.vcf_writer.build_header_lines (issue #1068) """

    info = None
    if info_dict:
        info = []
        for info_id, data in info_dict.items():
            data['id'] = info_id
            if data.get("number") is None:
                data["number"] = 1
            info.append(VCFInfoHeader(id=info_id, type=data["type"],
                                      description=data.get("description") or "", number=data["number"]))

    return build_header_lines(meta_lines=top_lines, info=info, formats=formats,
                              contig_lines=contig_lines, samples=samples)


def get_vcf_header_from_contigs(genome_build, info_dict=None, samples=None, use_accession=True):
    """ info_dict which contains ('number', 'type', 'description') """

    formats = [
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=A,Type=Integer,Description="Allelic Depths">',
        '##FORMAT=<ID=PL,Number=1,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">',
        '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1)">'
    ]
    contig_lines = get_contigs_header_lines(genome_build, use_accession=use_accession)
    return get_vcf_header_lines(info_dict=info_dict, formats=formats, contig_lines=contig_lines, samples=samples)
