# Generated by Django 4.2.10 on 2024-03-18 02:30
import os

from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_missing_fasta(genome_build_name):
    build_annotation = settings.ANNOTATION.get(genome_build_name)
    if build_annotation.get("enabled"):
        if vep_config := build_annotation.get("vep_config"):
            if vep_fasta_filename := vep_config.get("fasta"):
                if not os.path.exists(vep_fasta_filename):
                    return True
    return False


def _test_has_missing_37_fasta(apps):
    return _test_has_missing_fasta("GRCh37")


def _test_has_missing_38_fasta(apps):
    return _test_has_missing_fasta("GRCh37")


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0091_one_off_populate_missing_symbolic_hgvs"),
    ]

    operations = [
        ManualOperation.operation_other(args=[
            "Download Ensembl Fasta file - go to fasta dir and run: ${VARIANTGRID_DIR}/annotation/annotation_data/generate_annotation/vep_fasta_grch37.sh"],
            test=_test_has_missing_37_fasta),
        ManualOperation.operation_other(args=[
            "Download Ensembl Fasta file - go to fasta dir and run: ${VARIANTGRID_DIR}/annotation/annotation_data/generate_annotation/vep_fasta_grch38.sh"],
            test=_test_has_missing_38_fasta),
    ]