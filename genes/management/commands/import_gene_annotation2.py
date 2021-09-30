import gzip
import json
import logging
from typing import Dict, List, Set

from django.core.management.base import BaseCommand
from genes.models import HGNC, GeneSymbol, GeneAnnotationImport
from genes.models_enums import AnnotationConsortium
from library.file_utils import open_handle_gzip
from snpdb.models.models_genome import GenomeBuild


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

    def add_arguments(self, parser):
        consortia = [ac[1] for ac in AnnotationConsortium.choices]
        builds = [gb.name for gb in GenomeBuild.builds_with_annotation()]

        parser.add_argument('--dry-run', action='store_true', help="Don't actually modify anything")
        parser.add_argument('--genome-build', choices=builds, required=True)
        parser.add_argument('--annotation-consortium', choices=consortia, required=True)
        parser.add_argument('--release', required=False,
                            help="Make a release (to match VEP) store all gene/transcript versions")
        parser.add_argument('--save-merged-file', help="Write a file (requires pyreference-json)")
        group = parser.add_mutually_exclusive_group()
        group.add_argument('--pyreference-json', nargs="+", action="extend", help='PyReference JSON.gz')
        group.add_argument('--merged-json', help='Merged JSON (from multiple PyReference files)')

    def handle(self, *args, **options):
        if pyreference_json := options["pyreference_json"]:
            pyreference_data = []
            for prj_filename in pyreference_json:
                logging.info("Loading %s", prj_filename)
                with open_handle_gzip(prj_filename) as f:
                    pyreference_data.append(json.load(f))
            merged_data = self._convert_to_merged_data(pyreference_data)
            if save_merged_file := options["save_merged_file"]:
                with gzip.open(save_merged_file, 'w') as outfile:
                    json_str = json.dumps(merged_data)
                    outfile.write(json_str.encode('ascii'))
                    exit(0)
        elif merged_json := options["merged_json"]:
            with open_handle_gzip(merged_json) as f:
                merged_data = json.load(f)
        else:
            raise ValueError("You need to specify at least one of '--pyreference-json' or '--merged-json'")

        self._import_merged_data(merged_data)

    @staticmethod
    def _get_most_recent_transcripts(pyreference_data) -> List[Set]:
        transcripts_in_files = [set(prd["transcripts_by_id"]) for prd in pyreference_data]
        most_recent_transcripts = []
        all_transcripts = set()
        for tif in reversed(transcripts_in_files):
            unique_transcripts = tif - all_transcripts
            most_recent_transcripts.append(unique_transcripts)
            all_transcripts |= unique_transcripts
        most_recent_transcripts.reverse()
        return most_recent_transcripts

    def _convert_to_merged_data(self, pyreference_data: List[Dict]) -> List[Dict]:
        """ We want to make the minimal amount of data to insert - so only keep the last copy of transcripts """
        print("_import_merged_data")
        most_recent_transcripts = self._get_most_recent_transcripts(pyreference_data)

        merged_data = []
        for prd, transcripts in zip(pyreference_data, most_recent_transcripts):
            gene_version = {}
            transcript_versions = {}

            for gene_id, gene in prd["genes_by_id"].items():
                for transcript_accession in gene["transcripts"]:
                    if transcript_accession in transcripts:
                        gene_version[gene_id] = convert_gene_pyreference_to_gene_version_data(gene)

            for transcript_accession in transcripts:
                transcript = prd["transcripts_by_id"][transcript_accession]
                transcript_versions[transcript_accession] = convert_transcript_pyreference_to_pyhgvs(transcript)

            data = {
                "gene_annotation_import": prd["reference_gtf"],
                "gene_version": gene_version,
                "transcript_versions": transcript_versions,
            }
            merged_data.append(data)

        return merged_data


    def _import_merged_data(self, merged_data: List[Dict]):
        """
        """
        print("_import_merged_data")

        gene_symbols = []
        genes = []
        gene_versions = []
        transcripts = []
        transcript_versions = []

        for data in merged_data:
            import_data = data["gene_annotation_import"]
            logging.info("%s has %d transcripts", import_data, len(data["transcript_versions"]))
#            import_source = GeneAnnotationImport.objects.create(annotation_consortium=self.annotation_consortium,
#                                                                genome_build=self.genome_build,
#                                                                filename=import_data["path"],
#                                                                url=import_data["url"],
#                                                                file_md5sum=import_data["md5sum"])

            # Go through and create gene_symbols etc, assigning import_source

        # Now, break up into new / old?

        # bulk_create
        # bulk update (with import source)




def convert_gene_pyreference_to_gene_version_data(gene_data: Dict) -> Dict:
    gene_version_data = {
        'biotype': ",".join(gene_data["biotype"]),
        'description': gene_data.get("description"),
        'gene_symbol': gene_data["name"],
    }

    if hgnc_str := gene_data.get("HGNC"):
        # Has HGNC: (5 characters) at start of it
        gene_version_data["hgnc"] = hgnc_str[5:]
    return gene_version_data


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

