#!/usr/bin/env python3
"""
Fills the gnomAD columns of T2T variants annotated from gnomAD 4.1.t2t_liftover, which has no record where CHM13's
reference allele differs from GRCh38's, from the 4.1.t2t_liftover_ref_swaps file that adds them (#2788).

Runs on each VariantAnnotationVersion with backfilled_gnomad_ref_swaps False (migration 0018 sets it on the ACTIVE
T2T version): dumps its variants with no gnomAD AF, `bcftools annotate`s them against the settings' gnomAD file,
imports the result and sets the flag.
"""
import os
import subprocess
import tempfile

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from annotation.backfill_columns import (
    dump_annotated_variants,
    import_backfill_vcf,
    mark_backfilled,
    resolve_backfill_columns,
)
from annotation.enums import VariantAnnotationPipelineType
from annotation.models.models_version import VariantAnnotationVersion
from annotation.vep_config import VEPConfig, parse_gnomad_version_from_filename

BACKFILL_FLAG = "backfilled_gnomad_ref_swaps"
GNOMAD_VEP_CONFIG_KEY = "gnomad4"

# VariantGrid column -> INFO field in the gnomAD file
GNOMAD_INFO_FIELDS = {
    "gnomad_ac": "AC",
    "gnomad_af": "AF",
    "gnomad_an": "AN",
    "gnomad_afr_af": "AF_afr",
    "gnomad_amr_af": "AF_amr",
    "gnomad_asj_af": "AF_asj",
    "gnomad_eas_af": "AF_eas",
    "gnomad_fin_af": "AF_fin",
    "gnomad_mid_af": "AF_mid",
    "gnomad_nfe_af": "AF_nfe",
    "gnomad_oth_af": "AF_remaining",
    "gnomad_sas_af": "AF_sas",
    "gnomad_popmax": "grpmax",
    "gnomad_popmax_ac": "AC_grpmax",
    "gnomad_popmax_af": "AF_grpmax",
    "gnomad_popmax_an": "AN_grpmax",
    "gnomad_xy_ac": "AC_XY",
    "gnomad_xy_af": "AF_XY",
    "gnomad_xy_an": "AN_XY",
}


class Command(BaseCommand):
    category = "one-off"
    help = "Fill gnomAD on T2T variants annotated before the 4.1.t2t_liftover_ref_swaps file (#2788)"

    def add_arguments(self, parser):
        parser.add_argument("--dry-run", action="store_true",
                            help="Report how many rows the import would change, changing nothing")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        vavs = VariantAnnotationVersion.objects.filter(**{BACKFILL_FLAG: False}).order_by("pk")
        if not vavs.exists():
            self.stdout.write(f"No VariantAnnotationVersion has {BACKFILL_FLAG}=False")
            return

        for vav in vavs:
            self._backfill(vav, dry_run)

    def _backfill(self, vav: VariantAnnotationVersion, dry_run: bool):
        gnomad_filename = self._gnomad_filename(vav)
        targets = resolve_backfill_columns(vav, list(GNOMAD_INFO_FIELDS), info_overrides=GNOMAD_INFO_FIELDS)
        os.makedirs(settings.ANNOTATION_VCF_DUMP_DIR, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=settings.ANNOTATION_VCF_DUMP_DIR) as tmp_dir:
            dump_filename = os.path.join(tmp_dir, "dump.vcf.gz")
            count = dump_annotated_variants(vav, dump_filename, targets=targets, only_missing=True,
                                            pipeline_type=VariantAnnotationPipelineType.STANDARD)
            self.stdout.write(f"{vav}: {count} variants with no gnomAD values")
            annotated_filename = os.path.join(tmp_dir, "annotated.vcf.gz")
            self.stdout.write(f"bcftools annotate against '{gnomad_filename}' (streams the whole file)")
            self._annotate(vav, dump_filename, gnomad_filename, tmp_dir, annotated_filename)
            report = import_backfill_vcf(vav, annotated_filename, targets, dry_run=dry_run,
                                         emit=self.stdout.write)

        verb = "would change" if dry_run else "changed"
        self.stdout.write(f"{vav}: read {report['records_read']} records, {verb} {report['rows_changed']} rows")
        if not dry_run:
            mark_backfilled(vav, BACKFILL_FLAG)
            self.stdout.write(f"Set {BACKFILL_FLAG} on {vav}")

    @staticmethod
    def _gnomad_filename(vav: VariantAnnotationVersion) -> str:
        """ The file VEP annotates this version's new variants from, which the backfill has to match """
        gnomad_filename = VEPConfig(vav.genome_build)[GNOMAD_VEP_CONFIG_KEY]
        gnomad_version = parse_gnomad_version_from_filename(gnomad_filename)
        if gnomad_version != vav.gnomad:
            raise CommandError(f"{vav} has gnomAD {vav.gnomad} but settings vep_config['{GNOMAD_VEP_CONFIG_KEY}'] "
                               f"is {gnomad_version} ('{gnomad_filename}')")
        if not os.path.exists(gnomad_filename):
            raise CommandError(f"'{gnomad_filename}' does not exist - run 'python3 manage.py download_annotation_data'")
        return gnomad_filename

    @staticmethod
    def _annotate(vav: VariantAnnotationVersion, dump_filename: str, gnomad_filename: str, tmp_dir: str,
                  annotated_filename: str):
        """ The gnomAD file names contigs by accession, the dump by chromosome """
        chrom_mapping = os.path.join(settings.BASE_DIR, "snpdb", "genome",
                                     f"chrom_mapping_{vav.genome_build.name}.map")
        renamed_filename = os.path.join(tmp_dir, "renamed.vcf.gz")
        columns = ",".join(f"INFO/{info_field}" for info_field in GNOMAD_INFO_FIELDS.values())
        commands = [
            [settings.BCFTOOLS_COMMAND, "annotate", "--rename-chrs", chrom_mapping, dump_filename,
             "-Oz", "-o", renamed_filename],
            [settings.BCFTOOLS_COMMAND, "index", "--tbi", renamed_filename],
            [settings.BCFTOOLS_COMMAND, "annotate", "-a", gnomad_filename, "-c", columns, renamed_filename,
             "-Oz", "-o", annotated_filename],
        ]
        for command in commands:
            if (process := subprocess.run(command, check=False, capture_output=True, text=True)).returncode != 0:
                raise CommandError(f"'{' '.join(command)}' failed ({process.returncode}): {process.stderr}")
