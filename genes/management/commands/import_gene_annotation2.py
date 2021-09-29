import os
import re
from collections import Counter, defaultdict
from typing import Tuple, Optional, Dict

from django.core.management.base import BaseCommand
from pyhgvs.utils import read_genepred

from genes.cached_web_resource.refseq import retrieve_refseq_gene_summaries
from genes.gene_matching import GeneMatcher
from genes.models import GeneAnnotationImport, HGNC, \
    GeneSymbol, Gene, GeneVersion, Transcript, TranscriptVersion, GeneAnnotationRelease, ReleaseGeneVersion, \
    ReleaseTranscriptVersion
from genes.models_enums import AnnotationConsortium
from library.django_utils import highest_pk
from library.file_utils import open_handle_gzip
from library.utils import invert_dict
from snpdb.models.models_enums import SequenceRole
from snpdb.models.models_genome import GenomeBuild, GenomeFasta


class Command(BaseCommand):
    BATCH_SIZE = 2000

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.annotation_consortium = None
        self.genome_build = None
        self.contig_id_to_fasta = None
        self.hgnc_ids = set(HGNC.objects.values_list("pk", flat=True))
        # Known objects containers are updated with new inserts
        self.known_gene_symbols = set(GeneSymbol.objects.all().values_list("pk", flat=True))
        self.known_gene_versions_by_gene_id = defaultdict(dict)
        self.known_transcript_versions_by_transcript_id = defaultdict(dict)

    def add_arguments(self, parser):
        consortia = [ac[1] for ac in AnnotationConsortium.choices]
        builds = [gb.name for gb in GenomeBuild.builds_with_annotation()]

        parser.add_argument('--dry-run', action='store_true', help="Don't actually modify anything")
        parser.add_argument('--genome-build', choices=builds, required=True)
        parser.add_argument('--annotation-consortium', choices=consortia, required=True)
        parser.add_argument('--release', required=False,
                            help="Make a release (to match VEP) store all gene/transcript versions")
        parser.add_argument('filename', help="PyReference JSON.gz")

    def handle(self, *args, **options):
        filename = options["filename"]


def convert_transcript_pyreference_to_pyhgvs(transcript_data: Dict) -> Dict:
    start = transcript_data["start"]
    end = transcript_data["stop"]
    strand = transcript_data["strand"]
    # PyHGVS has cds_start/cds_end be equal to start/end for non-coding transcripts
    cds_start = transcript_data.get("cds_start", start)
    cds_end = transcript_data.get("cds_end", end)
    # PyHGVS exons are in genomic order, pyhgvs are in stranded
    features = transcript_data["features_by_type"]
    exons = [[ed["start"], ed["stop"]] for ed in features["exon"]]
    if strand == '-':
        exons.reverse()

    return {
        'chrom': transcript_data["chr"],
        'start': start,
        'end': end,
        'strand': strand,
        'cds_start': cds_start,
        'cds_end': cds_end,
        'exons': exons,
    }
