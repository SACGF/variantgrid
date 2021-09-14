import logging

from django.core.management import BaseCommand

from genes.models import TranscriptVersion, GeneAnnotationImport, TranscriptVersionSequenceInfo
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    def handle(self, *args, **options):

        have_tvi = set()
        for tvi in TranscriptVersionSequenceInfo.objects.all():
            accession = f"{tvi.transcript_id}.{tvi.version}"
            have_tvi.add(accession)

        tv_lengths_37 = {}
        tv_lengths_38 = {}
        for tv in TranscriptVersion.objects.filter(transcript__annotation_consortium='R', alignment_gap=False,
                                                   genome_build=GenomeBuild.grch37()):
            tv_lengths_37[tv.accession] = tv.length

        for tv in TranscriptVersion.objects.filter(transcript__annotation_consortium='R', alignment_gap=False,
                                                   genome_build=GenomeBuild.grch38()):
            tv_lengths_38[tv.accession] = tv.length

        different_lengths = set()
        for accession in set(tv_lengths_37) & set(tv_lengths_38):
            length_37 = tv_lengths_37[accession]
            length_38 = tv_lengths_38[accession]
            if length_37 != length_38:
                different_lengths.add(accession)

        only_one_build = set(tv_lengths_37) ^ set(tv_lengths_38)
        # Imports w/o GFFs (only from genePred) don't have alignment info (gap_count) so we can't detect alignment gaps
        no_gff_imports = GeneAnnotationImport.objects.filter(filename__contains='genePred')
        no_gff_tvs = {tv.accession for tv in TranscriptVersion.objects.filter(import_source__in=no_gff_imports)}

        need_to_retrieve = (different_lengths | only_one_build | no_gff_tvs) - have_tvi
        num_to_retrieve = len(need_to_retrieve)
        logging.info("Retrieving %d Transcript Version Sequence Info records (takes ~1min per 1000)", num_to_retrieve)
        TranscriptVersionSequenceInfo.get_refseq_transcript_versions(need_to_retrieve)
