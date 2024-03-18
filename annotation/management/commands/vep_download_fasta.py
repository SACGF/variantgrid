import os
import subprocess

from django.conf import settings
from django.core.management.base import BaseCommand

from annotation.models import VariantAnnotationVersion
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    """
        We switched to using VEP fastas
        This ensures they have been downloaded
    """

    @staticmethod
    def _download_and_bgzip_fasta(directory, url):
        filename = url.split('/')[-1]
        command = f"wget --quiet -O - {url} | gzip -d | bgzip > {filename}; samtools faidx {filename}"
        print(f"Executing: '{command}'")
        subprocess.check_output(command, shell=True, cwd=directory)

    @staticmethod
    def _get_url(genome_build, vep_version):
        if genome_build == "GRCh37":
            # Always the same
            return "https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
        elif genome_build == "GRCh38":
            return "http://ftp.ensembl.org/pub/release-{vep_version}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
        else:
            raise ValueError(f"Unknown genome build {genome_build}")

    def handle(self, *args, **options):

        for genome_build in GenomeBuild.builds_with_annotation():
            vav = VariantAnnotationVersion.latest(genome_build)

            org_name = "homo_sapiens"
            if genome_build.annotation_consortium == AnnotationConsortium.REFSEQ:
                org_name += "_refseq"

            dir_name = os.path.join(settings.ANNOTATION_VEP_CACHE_DIR, org_name, f"{vav.vep}_{genome_build.name}")
            url = self._get_url(genome_build.name, vav.vep)
            filename = url.split("/")[-1]
            fasta_file = os.path.join(dir_name, filename)
            print(f"{genome_build}/{vav.get_annotation_consortium_display()} - looking for VEP fasta: '{fasta_file}'")
            if os.path.exists(fasta_file):
                print(f"file already exists")
            else:
                print(f"Downloading '{fasta_file}'")
                self._download_and_bgzip_fasta(dir_name, url)
