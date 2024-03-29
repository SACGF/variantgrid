"""
Test program to run VEP

"""

import logging
import os

from django.conf import settings
from django.core.management.base import BaseCommand

from annotation.models import VariantAnnotationPipelineType
from annotation.vep_annotation import run_vep, VEPConfig
from snpdb.models.models_genome import GenomeBuild

DO_SMALL = False


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--test', action='store_true')
        parser.add_argument('--sv', action='store_true')
        parser.add_argument('--genome-build', required=True)

    def handle(self, *args, **options):
        test = options["test"]
        sv = options["sv"]
        build_name = options["genome_build"]
        genome_build = GenomeBuild.get_name_or_alias(build_name)
        vc = VEPConfig(genome_build)

        if test:
            print("Re-generating VCF for unit test")
            vep_suffix = "vep_annotated"
            test_vcf_filename = f"test_{genome_build.name.lower()}"
            if sv:
                test_vcf_filename += "_symbolic_alt"
            unit_test_dir = os.path.join(settings.BASE_DIR, "annotation/tests/test_data")
            vcf_filename = os.path.join(unit_test_dir, f"{test_vcf_filename}.vcf")
            output_dir = unit_test_dir
            base_name = f"test_columns_version{vc.columns_version}_{genome_build.name.lower()}"
        else:
            vep_suffix = f"vep_annotated_{genome_build.name}"
            output_dir = settings.ANNOTATION_VCF_DUMP_DIR

            if DO_SMALL:
                base_name = "dump_small"
                VEP_EXAMPLES_DIR = os.path.join(settings.ANNOTATION_VEP_BASE_DIR, "examples")
                vcf_filename = os.path.join(VEP_EXAMPLES_DIR, f"{base_name}.vcf")
            else:
                base_name = "test"  # "dump_small"
                vcf_filename = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, f"{base_name}.vcf")

        if sv:
            pipeline_type = VariantAnnotationPipelineType.STRUCTURAL_VARIANT
            base_name += "_sv"
        else:
            pipeline_type = VariantAnnotationPipelineType.STANDARD

        output_filename = os.path.join(output_dir, f"{base_name}.{vep_suffix}.vcf.gz")
        return_code, std_out, std_err = run_vep(vcf_filename, output_filename,
                                                genome_build, genome_build.annotation_consortium,
                                                pipeline_type)
        if return_code != 0:
            logging.info(std_out)
            logging.error(std_err)
        else:
            print(f"wrote: {output_filename}")
