#!/usr/bin/env python3
"""
Score structural-variant conservation (phastCons/phyloP _max) with pyBigWig (#1657).

Standalone entry point around annotation.sv_conservation: reads an SV VCF (the dumped or annotated VCF,
whose INFO carries variant_id/END/SVLEN), computes the 4 conservation _max columns across the build's
bigWig tracks, and writes the sidecar TSV the import path consumes (or an explicit --output). The bigWig
tracks are scored in parallel; --threads defaults to settings.ANNOTATION_VEP_FORK (the setting that also
drives VEP --fork).
"""
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from annotation.sv_conservation import (
    conservation_sidecar_filename,
    get_sv_conservation_tracks,
    score_sv_vcf,
    write_conservation_sidecar,
)
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    help = "Score SV conservation (phastCons/phyloP max) with pyBigWig (#1657)"

    def add_arguments(self, parser):
        parser.add_argument("--vcf", required=True, help="SV VCF (variant_id/END/SVLEN in INFO)")
        parser.add_argument("--genome-build", required=True)
        parser.add_argument("--output", help="Sidecar TSV path (default: <vcf>.conservation.tsv)")
        parser.add_argument("--threads", type=int, default=None,
                            help="Thread pool size (default: settings.ANNOTATION_VEP_FORK)")

    def handle(self, *args, **options):
        genome_build = GenomeBuild.get_name_or_alias(options["genome_build"])
        tracks = get_sv_conservation_tracks(genome_build)
        if not tracks:
            raise CommandError(f"No conservation bigWig tracks configured for {genome_build}")

        threads = options["threads"]
        if threads is None:
            threads = settings.ANNOTATION_VEP_FORK or 1

        vcf_filename = options["vcf"]
        results = score_sv_vcf(vcf_filename, genome_build, threads=threads)
        output = options["output"] or conservation_sidecar_filename(vcf_filename)
        write_conservation_sidecar(output, results, tracks)
        self.stdout.write(f"Scored {len(results)} SVs across {len(tracks)} tracks "
                          f"({threads} threads) -> {output}")
