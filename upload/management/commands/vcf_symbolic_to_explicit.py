import sys

from django.conf import settings
from django.core.management.base import BaseCommand

from library.genomics.vcf_enums import VCFColumns
from snpdb.models import GenomeBuild, VariantCoordinate


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--vcf', help='VCF file, default: - (stdin)', default="-")
        parser.add_argument('--genome-build', help='GenomeBuild name', required=True)

    def handle(self, *args, **options):
        vcf_filename = options["vcf"]
        build_name = options["genome_build"]

        genome_build = GenomeBuild.get_name_or_alias(build_name)
        if vcf_filename == '-':
            f = sys.stdin
        else:
            f = open(vcf_filename)

        vcf_header_lines = []
        for line in f:
            if line[0] == '#':
                vcf_header_lines.append(line)
                sys.stdout.write(line)
            else:
                columns = line.strip().split("\t")
                chrom = columns[VCFColumns.CHROM]
                position = columns[VCFColumns.POS]
                ref = columns[VCFColumns.REF]
                alt = columns[VCFColumns.ALT]
                info = columns[VCFColumns.INFO]
                info_dict = {}
                if info and info != '.':
                    for il in info.split(";"):
                        if "=" in il:
                            k, v = il.split("=")
                            info_dict[k] = v
                        else:
                            info_dict[il] = None
                    if svlen := info_dict.pop("SVLEN", None):
                        if alt in settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES:
                            vc = VariantCoordinate(chrom=chrom, position=position, ref=ref, alt=alt, svlen=svlen)
                            vc_e = vc.as_external_explicit(genome_build)
                            chrom, position, ref, alt, _ = vc_e
                            columns[VCFColumns.CHROM] = chrom
                            columns[VCFColumns.POS] = str(position)
                            columns[VCFColumns.REF] = ref
                            columns[VCFColumns.ALT] = alt
                            info_list = []
                            for k, v in info_dict.items():
                                if v is not None:
                                    info_list.append(f"{k}={v}")
                                else:
                                    info_list.append(k)  # Flag
                            columns[VCFColumns.INFO] = ";".join(info_list)
                sys.stdout.write("\t".join(columns) + "\n")
