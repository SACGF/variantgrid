"""
Rewrite Ensembl's SpliceVault VEP plugin file to MANE RefSeq transcript ids, for a RefSeq build (#2426).
@see annotation/splicevault.py
"""
import os

from django.core.management.base import BaseCommand, CommandError

from annotation.splicevault import (
    convert_splicevault_file,
    default_output_filename,
    mane_ensembl_to_refseq,
)
from genes.enums import AnnotationConsortium
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    category = "import"
    help = "Convert SpliceVault_data_<build>.tsv.gz from Ensembl ENST to MANE RefSeq NM ids, bgzip + tabix"

    def add_arguments(self, parser):
        parser.add_argument("--genome-build", required=True)
        parser.add_argument("--input", required=True,
                            help="Ensembl download, eg SpliceVault_data_GRCh38.tsv.gz")
        parser.add_argument("--output",
                            help="Default: alongside the input, eg SpliceVault_data_GRCh38.refseq_mane.tsv.gz. "
                                 "Its basename is the VariantAnnotationVersion 'splicevault' pin")

    def handle(self, *args, **options):
        genome_build = GenomeBuild.get_name_or_alias(options["genome_build"])
        if genome_build.annotation_consortium != AnnotationConsortium.REFSEQ:
            raise CommandError(f"{genome_build} is annotated with {genome_build.annotation_consortium.label} "
                               f"transcripts - configure the Ensembl file as it is")
        input_filename = options["input"]
        if not os.path.exists(input_filename):
            raise CommandError(f"{input_filename} not found")
        output_filename = options["output"] or default_output_filename(input_filename)
        if os.path.abspath(output_filename) == os.path.abspath(input_filename):
            raise CommandError("--output would overwrite --input")

        ensembl_to_refseq = mane_ensembl_to_refseq()
        if not ensembl_to_refseq:
            raise CommandError("No MANE transcript pairs - load the MANE CachedWebResource first")

        conversion = convert_splicevault_file(input_filename, output_filename, ensembl_to_refseq)
        self.stdout.write(f"Wrote {output_filename} (+ .tbi): {conversion.rows_written:,} of "
                          f"{conversion.rows_read:,} rows, {conversion.transcripts_converted:,} transcripts "
                          f"converted to MANE RefSeq. Dropped {conversion.rows_dropped:,} rows on "
                          f"{conversion.transcripts_dropped:,} transcripts with no MANE pair")
