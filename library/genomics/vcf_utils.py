import logging
import operator
import os
import re
from collections import defaultdict
from typing import Optional, Iterable, IO, Union

import cyvcf2
import vcf

from library.genomics.vcf_enums import VCFSymbolicAllele
from library.utils import open_handle_gzip, open_file_or_filename
from snpdb.models import Variant, Sequence, GenomeFasta, SequenceRole, VariantCoordinate


def cyvcf2_header_types(cyvcf2_reader) -> defaultdict:
    header_types = defaultdict(dict)
    for h in cyvcf2_reader.header_iter():
        info = h.info()
        h_id = info.get("ID")
        if h_id:  # Not much use w/o this
            header_types[h.type][h_id] = info
    return header_types


def cyvcf2_header_get(cyvcf2_reader, key, default=None):
    try:
        header_dict = cyvcf2_reader.get_header_type(key)
    except:
        header_dict = {}
    return header_dict.get(key, default)


def cyvcf2_get_contig_lengths_dict(cyvcf2_reader):
    try:
        # For some versions of HTSLib it seems that if the VCF has eg contigs as "chr1" and records as "1" it will
        # add in fake contigs into the header - so skip those
        return {k: v for k,v in zip(cyvcf2_reader.seqnames, cyvcf2_reader.seqlens) if v != -1}
    except AttributeError:
        return {}


def write_vcf_from_variant_coordinates(file_or_filename: Union[str, IO], variant_coordinates: Iterable[VariantCoordinate],
                                       vcf_ids: Optional[Iterable[str]] = None, header_lines: list[str] = None):
    """ If vcf_ids is supplied it must be in sync (record for record) with variant_coordinates """

    if header_lines is None:
        header_lines = []

    if vcf_ids is None:
        vcf_ids = ("." for _ in variant_coordinates)

    # Put them together and then sort
    vc_and_ids = ((vc, vcf_id) for vc, vcf_id in zip(variant_coordinates, vcf_ids))
    vc_and_ids = sorted(vc_and_ids, key=operator.itemgetter(0))

    info = [
        '##INFO=<ID=END,Number=.,Type=Integer,Description="Stop position of the interval">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    ]
    columns = "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    header = "\n".join(["##fileformat=VCFv4.1", "##source=VariantGrid"] + info + header_lines + ["#" + columns])
    with open_file_or_filename(file_or_filename, mode="wt", encoding="utf-8") as f:
        f.write(header + "\n")
        for variant_coordinate, vcf_id in vc_and_ids:
            (chrom, position, ref, alt, svlen) = variant_coordinate
            if Sequence.allele_is_symbolic(alt):
                svtype = alt[1:-1]  # Strip off brackets
                info_dict = {
                    "SVLEN": svlen,
                    "SVTYPE": svtype,
                    "END": variant_coordinate.end,
                }
                record_info = ";".join([f"{k}={v}" for k,v in info_dict.items()])
            else:
                record_info = "."
                if alt == Variant.REFERENCE_ALT:
                    alt = "."
            line = "\t".join((chrom, str(position), str(vcf_id), ref, alt, ".", ".", record_info))
            f.write(line + "\n")


def vcf_to_variant_coordinates_and_records(filename: str) -> Iterable[tuple[VariantCoordinate, cyvcf2.Variant]]:
    reader = cyvcf2.Reader(filename)  # Can only handle filenames
    for v in reader:
        for alt in v.ALT:
            variant_coordinate = VariantCoordinate(chrom=v.CHROM, position=v.POS, ref=v.REF, alt=alt)
            yield variant_coordinate, v


def vcf_to_variant_coordinates(filename: str) -> Iterable[VariantCoordinate]:
    for vc, _ in vcf_to_variant_coordinates_and_records(filename):
        yield vc


def get_variant_caller_and_version_from_vcf(filename) -> tuple[str, str]:
    variant_caller = None
    version = None

    if os.path.exists(filename):
        reader = vcf.Reader(filename=filename)

        if source_list := reader.metadata.get("source"):
            for source in source_list:
                # Match source = "freeBayes v1.3.5" or "VarDict_v1.8.2"
                if m := re.match(r"(.*?)[ _]v([\d\\.]+)", source):
                    variant_caller, version = m.groups()
                    break

        if gatk_commandline := reader.metadata.get("GATKCommandLine"):
            variant_caller = "GATK"
            for commandline in gatk_commandline:
                if caller_id := commandline.get("ID"):
                    if caller_id == 'HaplotypeCaller':
                        caller_id = "GATK"  # Just stay with GATK
                    variant_caller = caller_id

                if version := commandline.get("Version"):
                    version = version.replace('"', "")  # Strip quotes

    return variant_caller, version


def vcf_allele_is_symbolic(allele: str) -> bool:
    return allele.startswith("<") and allele.endswith(">")


def vcf_get_ref_alt_svlen_and_modification(variant: cyvcf2.Variant, old_variant_info: str) -> tuple[str, str, Optional[int], Optional[str]]:
    """ old_variant_info: name of INFO created via bcftools norm --old-rec-tag=INFO_NAME

        Looks like 'NC_000022.11|1544211|N|<DUP>,<DEL>|2' - we use the final index (1/2)


    """
    ref = variant.REF.strip().upper()
    if variant.ALT:
        alt = variant.ALT[0].strip().upper()
    else:
        alt = Variant.REFERENCE_ALT

    modification = None
    if Sequence.allele_is_symbolic(ref) or Sequence.allele_is_symbolic(alt):
        # Need to provide END or SVLEN
        if svlen_info := variant.INFO.get('SVLEN'):
            # issue #1268 - bcftools norm doesn't split multi-allelic SVLEN
            if isinstance(svlen_info, tuple):
                if not old_variant_info:
                    raise ValueError("SVLEN is a tuple, you need to pass in 'old_variant_info' to resolve these")

                if old_variant := variant.INFO.get(old_variant_info):
                    old_index_str = old_variant.split("|")[-1]
                    try:
                        old_index = int(old_index_str)
                    except ValueError as ve:
                        raise ValueError(f"SVLEN is a tuple (from bcftools norm multi-allelic split). Could not obtain"
                                         f" old alt index from INFO/{old_variant_info}='{old_variant}'") from ve
                    try:
                        svlen_info = svlen_info[old_index - 1]
                    except IndexError as ie:
                        raise ValueError(f"SVLEN is a tuple of length={len(svlen_info)}. Obtained 1 based index "
                                         f"{old_index=} from INFO/{old_variant_info}='{old_variant}' "
                                         f"which was out of range") from ie
                else:
                    raise ValueError(f"SVLEN is a tuple. This requires INFO='{old_variant_info}', "
                                     f"expected split multi-allelic by bcftools with that as --old-rec-tag")

            svlen = int(svlen_info)
            if alt == VCFSymbolicAllele.DEL and svlen > 0:
                # issue #1245 - Manta SV produces <DEL> with positive SVLEN
                svlen = -svlen
                modification = f"SVLEN - inverted positive value for {alt=}"
        elif end_info := variant.INFO.get('END'):
            svlen = int(end_info) - variant.POS
        else:
            raise ValueError(f"SVLEN or END info field MUST be provided for symbolic (ie '<x>') {ref=},{alt=}")
    else:
        svlen = None
    return ref, alt, svlen, modification


def get_vcf_header_contig_lines(contigs: list[tuple]) -> list[str]:
    header_lines = []
    for contig, length, assembly in contigs:
        line = f"##contig=<ID={contig},length={length},assembly={assembly}>"
        header_lines.append(line)
    return header_lines


def get_contigs_header_lines(genome_build, standard_only=True, use_accession=True, contig_allow_list: set = None) -> list[str]:
    """ use_accession: If True - write contigs like 'NC_000004.12' if False then '4' """
    if standard_only:
        contig_qs = genome_build.standard_contigs
    else:
        contig_qs = genome_build.contigs

    if contig_allow_list:
        logging.info("Writing contigs header for contigs: %s", contig_allow_list)

    reference_lines = [
        f"##reference={genome_build.name}>"
    ]

    contigs = []
    for contig in contig_qs:
        if use_accession:
            contig_name = contig.refseq_accession
        else:
            contig_name = contig.name
        if contig_allow_list is not None:
            if contig_name not in contig_allow_list:
                continue
        contigs.append((contig_name, contig.length, genome_build.name))
    contig_lines = get_vcf_header_contig_lines(contigs)
    return reference_lines + contig_lines


def write_cleaned_vcf_header(genome_build, source_vcf_filename: str, output_filename: str,
                             new_info_lines: list[str] = None, standard_contigs_only=True):
    contig_regex = re.compile(r"^##contig=<ID=(.+),length=(\d+)")

    header_lines = []
    with open_handle_gzip(source_vcf_filename, "rt") as in_f:
        for line in in_f:
            if not line.startswith("#"):
                break  # End of header
            header_lines.append(line.strip())

    # These are used to validate contigs in header
    genome_fasta = GenomeFasta.get_for_genome_build(genome_build)
    contig_lengths = dict(genome_build.contigs.values_list("pk", "length"))
    contig_to_fasta_names = genome_fasta.get_contig_id_to_name_mappings()

    with open(output_filename, "w") as f:
        found_column_names_line = False
        for line in header_lines:
            if line.startswith("#CHROM"):
                found_column_names_line = True
                # This is where we dump the new stuff
                if new_info_lines:
                    for new_info_line in new_info_lines:
                        f.write(new_info_line + "\n")
                for contig_line in get_contigs_header_lines(genome_build, standard_only=standard_contigs_only):
                    f.write(contig_line + "\n")

            elif m := contig_regex.match(line):
                # Strip existing contig lines from header - though check they match so we don't get build swaps
                contig_name, provided_contig_length = m.groups()
                if contig := genome_build.chrom_contig_mappings.get(contig_name):
                    if standard_contigs_only:
                        if contig.role != SequenceRole.ASSEMBLED_MOLECULE:
                            continue  # No validation

                    if fasta_chrom := contig_to_fasta_names.get(contig.pk):
                        provided_contig_length = int(provided_contig_length)
                        ref_contig_length = contig_lengths[contig.pk]
                        if provided_contig_length != ref_contig_length:
                            msg = f"VCF header contig '{contig_name}' (length={provided_contig_length}) has " + \
                                f"different length than ref contig {fasta_chrom} (length={ref_contig_length})"
                            raise ValueError(msg)
                # Don't write any old contigs
                continue

            f.write(line + "\n")

        if not found_column_names_line:
            raise ValueError("VCF header was missing line starting with '#CHROM'")
