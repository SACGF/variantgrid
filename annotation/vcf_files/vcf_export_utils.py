from django.http.response import StreamingHttpResponse
from io import StringIO
from vcf.model import _Call, _Record, make_calldata_tuple, _Substitution
from vcf.parser import Reader, Writer

from annotation.models import ColumnVCFInfo, VCFInfoTypes
from library.jqgrid_export import StashFile
from patients.models_enums import Zygosity


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


def get_column_vcf_info():
    column_vcf_info = {}
    for cvi in ColumnVCFInfo.objects.all().values('column__variant_column', 'info_id', 'number', 'type', 'description'):
        lookup = dict(VCFInfoTypes.CHOICES)
        type_display = lookup[cvi['type']]
        cvi['type'] = type_display

        index = cvi['column__variant_column']
        column_vcf_info[index] = cvi
    return column_vcf_info


def grid_item_to_vcf_record(info_dict, obj, sample_ids, sample_names):  # , get_genotype_from_expanded_zygosity):
    CHROM = obj.get("locus__contig__name", ".")
    POS = obj.get("locus__position", ".")
    ID = obj.get("variantannotation__dbsnp_rs_id")
    REF = obj.get("locus__ref__seq", ".")
    ALT = obj.get("alt__seq", ".")
    QUAL = '.'  # QUAL = obj.get("annotation__quality", ".")
    FILTER = None
    INFO = {}

    for info_id, data in info_dict.items():
        col = data['column__variant_column']
        val = obj.get(col)
        if val:
            INFO[info_id] = val

    FORMAT = None
    MY_FORMAT = ['GT', 'AD', 'AF', 'PL', 'DP', 'GQ']
    CallData = make_calldata_tuple(MY_FORMAT)
    sample_indexes = {}
    samples = []

    if sample_ids:
        FORMAT = ':'.join(MY_FORMAT)

    alts = [_Substitution(ALT)]
    ALT = alts
    record = _Record(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, sample_indexes)

    if sample_ids:
        for i, (sample_id, sample) in enumerate(zip(sample_ids, sample_names)):
            ad = obj[f"{sample_id}_samples_allele_depth"]
            zygosity = obj[f"{sample_id}_samples_zygosity"]
            gt = Zygosity.get_genotype_from_expanded_zygosity(zygosity)
            pl = obj[f"{sample_id}_samples_phred_likelihood"]
            dp = obj[f"{sample_id}_samples_read_depth"]
            gq = obj[f"{sample_id}_samples_genotype_quality"]
            af = obj[f"{sample_id}_samples_allele_frequency"]
            # TODO: Need to grab information for reference base to be able to properly fill in this data.
            data_args = {'AD': ['.', ad],
                         'GT': gt,
                         'PL': ['.', pl],
                         'DP': ['.', dp],
                         'GQ': ['.', gq],
                         'AF': ['.', af]}

            data = CallData(**data_args)
            call = _Call(record, sample, data)
            samples.append(call)
            sample_indexes[sample] = i

        record.samples = samples

    return record


def get_colmodel_info_dict(colmodels):
    column_vcf_info = get_column_vcf_info()

    info_dict = {}
    for c in colmodels:
        name = c['name']
        col_info = column_vcf_info.get(name)
        if col_info:
            col_info['number'] = col_info['number'] or '.'

            info_id = col_info['info_id']
            info_dict[info_id] = col_info
    return info_dict


def colmodels_to_vcf_header(genome_build, info_dict, samples):
    """ returns file which contains header """

    header_lines = get_vcf_header_from_contigs(genome_build, info_dict, samples)
    return StringIO('\n'.join(header_lines))


def grid_export_vcf(filename, genome_build, colmodels, items, sample_ids, sample_names_by_id):
    samples = [sample_names_by_id[s_id] for s_id in sample_ids]

    info_dict = get_colmodel_info_dict(colmodels)
    vcf_template_file = colmodels_to_vcf_header(genome_build, info_dict, samples)
    vcf_reader = Reader(vcf_template_file)

    pseudo_buffer = StashFile()

    vcf_writer = Writer(pseudo_buffer, vcf_reader)

    def iter_row_writer():

        for obj in items:
            record = grid_item_to_vcf_record(info_dict, obj, sample_ids, samples)
            vcf_writer.write_record(record)
            yield pseudo_buffer.value

    response = StreamingHttpResponse(iter_row_writer(), content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="{filename}.vcf"'
    return response
