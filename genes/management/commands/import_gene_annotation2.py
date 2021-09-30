import json
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
    """
        * Insert new gene symbols
        * Insert new Gene
        * Insert new GeneVerison
        * Insert new Transcripts
        * Insert new TranscriptVersions

        UPDATES

        GeneVersion - symbol + import
        TranscriptVersion - data, gene_version + import


        -------------


        ReleaseGeneVersion
        ReleaseTranscriptVersion
    """
    BATCH_SIZE = 2000

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.annotation_consortium = None
        self.genome_build = None
        self.contig_id_to_fasta = None
        self.hgnc_ids = set(HGNC.objects.values_list("pk", flat=True))
        # Known objects containers are updated with new inserts
        self.known_gene_symbols = set(GeneSymbol.objects.all().values_list("pk", flat=True))

    def add_arguments(self, parser):
        consortia = [ac[1] for ac in AnnotationConsortium.choices]
        builds = [gb.name for gb in GenomeBuild.builds_with_annotation()]

        parser.add_argument('--dry-run', action='store_true', help="Don't actually modify anything")
        parser.add_argument('--genome-build', choices=builds, required=True)
        parser.add_argument('--annotation-consortium', choices=consortia, required=True)
        parser.add_argument('--release', required=False,
                            help="Make a release (to match VEP) store all gene/transcript versions")
        group = parser.add_mutually_exclusive_group()
        group.add_argument('--pyreference-json', help='PyReference JSON.gz')
        group.add_argument('--merged-json', help='Merged JSON (from multiple PyReference files)')

    def handle(self, *args, **options):
        if pyreference_json := options["merged_json"]:
            with open_handle_gzip(pyreference_json) as f:
                pyreference_data = json.load(f)
            merged_data = self._convert_to_merged_data(pyreference_data)
        else:
            merged_json = options["merged_json"]
            with open_handle_gzip(merged_json) as f:
                merged_data = json.load(f)

        self._import_merged_data(merged_data)

    def _convert_to_merged_data(self, pyreference_data: Dict) -> Dict:
        pass

    def _import_merged_data(self, data: Dict):
        """
            [{
                "gene_annotation_import": {"filename": "", "url": "", "file_md5sum": ""},
                "gene_symbols": [],
                "gene": [],
                "gene_version": [],
                "transcripts": [],
                "transcript_versions": [],
             },
             }
        """
        pass


def convert_transcript_pyreference_to_pyhgvs(transcript_data: Dict) -> Dict:
    start = transcript_data["start"]
    end = transcript_data["stop"]
    strand = transcript_data["strand"]
    # PyHGVS has cds_start/cds_end be equal to start/end for non-coding transcripts
    cds_start = transcript_data.get("cds_start", start)
    cds_end = transcript_data.get("cds_end", end)
    # PyHGVS exons are in genomic order, PyReference are in stranded
    features = transcript_data["features_by_type"]
    exons = [[ed["start"], ed["stop"]] for ed in features["exon"]]
    cdna_match = [cdm.get("gap") for cdm in features.get("cDNA_match", [])]

    if strand == '-':
        exons.reverse()
        cdna_match.reverse()

    pyhgvs_data = {
        'chrom': transcript_data["chr"],
        'start': start,
        'end': end,
        'strand': strand,
        'cds_start': cds_start,
        'cds_end': cds_end,
        'exons': exons,
    }

    # Optional stuff
    if cdna_match:
        pyhgvs_data["cdna_match"] = cdna_match
    if transcript_data.get("partial"):
        pyhgvs_data["partial"] = 1

    return pyhgvs_data

