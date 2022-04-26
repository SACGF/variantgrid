import json
import logging
from typing import Dict, Tuple

from django.core.management.base import BaseCommand
from django.db.models.functions import Upper

from genes.cached_web_resource.refseq import retrieve_refseq_gene_summaries
from genes.gene_matching import GeneMatcher
from genes.models import GeneSymbol, GeneAnnotationImport, Gene, GeneVersion, TranscriptVersion, Transcript, HGNC, \
    GeneAnnotationRelease, ReleaseGeneVersion, ReleaseTranscriptVersion
from genes.models_enums import AnnotationConsortium
from library.file_utils import open_handle_gzip
from library.utils import invert_dict
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    BATCH_SIZE = 2000

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.import_source_by_url = {}

    def add_arguments(self, parser):
        consortia = [ac[1] for ac in AnnotationConsortium.choices]
        builds = [gb.name for gb in GenomeBuild.builds_with_annotation()]

        parser.add_argument('--genome-build', choices=builds, required=True)
        parser.add_argument('--annotation-consortium', choices=consortia, required=True)
        parser.add_argument('--release', required=False,
                            help="Make a release (to match VEP) store all gene/transcript versions")
        parser.add_argument('--json-file', required=True,
                            help='cdot JSON.gz')
        parser.add_argument('--clear-obsolete', action='store_true', help='Clear old transcripts')

    def handle(self, *args, **options):
        build_name = options["genome_build"]
        annotation_consortium_name = options["annotation_consortium"]
        genome_build = GenomeBuild.get_name_or_alias(build_name)
        ac_dict = invert_dict(dict(AnnotationConsortium.choices))
        annotation_consortium = ac_dict[annotation_consortium_name]

        for import_source in GeneAnnotationImport.objects.filter(annotation_consortium=annotation_consortium,
                                                                 genome_build=genome_build):
            self.import_source_by_url[import_source.url] = import_source

        json_file = options["json_file"]
        with open_handle_gzip(json_file) as f:
            cdot_data = json.load(f)
            cdot_version = cdot_data.get("cdot_version")
            if cdot_version is None:
                raise ValueError("JSON does not contain 'cdot_version' key")
            print(f"JSON uses cdot version {cdot_version}")

        if release_version := options["release"]:
            self._create_release(genome_build, annotation_consortium, release_version, cdot_data)
        else:
            self._import_merged_data(genome_build, annotation_consortium, cdot_data)

        if options["clear_obsolete"]:
            print("Clearing old Transcript Versions")
            tv_qs = TranscriptVersion.objects.filter(transcript__annotation_consortium=annotation_consortium,
                                                     genome_build=genome_build)
            ret = tv_qs.filter(data__genome_builds__isnull=True).delete()
            print(f"Deleted: {ret}")

    def _get_import_source_by_url(self, genome_build, annotation_consortium, url):
        import_source = self.import_source_by_url.get(url)
        if not import_source:
            import_source = GeneAnnotationImport.objects.create(annotation_consortium=annotation_consortium,
                                                                genome_build=genome_build,
                                                                url=url)
            self.import_source_by_url[url] = import_source
        return import_source

    def _create_release(self, genome_build: GenomeBuild, annotation_consortium, release_version, cdot_data):
        """ A GeneAnnotationRelease doesn't change/store transcript data, but does keep track of eg what
            symbols are used and how things are linked together """

        # A release should be from a single GTF - so all URLs should be the same, so take any one
        random_transcript = next(iter(cdot_data["transcripts"].values()))
        url = random_transcript["genome_builds"][genome_build.name]["url"]
        import_source = self.import_source_by_url[url]  # For a release, this must be there as it was imported before
        release, created = GeneAnnotationRelease.objects.update_or_create(version=release_version,
                                                                          genome_build=genome_build,
                                                                          annotation_consortium=annotation_consortium,
                                                                          defaults={
                                                                              "gene_annotation_import": import_source
                                                                          })
        if not created:
            print("Release exists - clearing existing data")
            release.releasegeneversion_set.all().delete()
            release.releasetranscriptversion_set.all().delete()

        gene_version_ids_by_accession, transcript_version_ids_by_accession = self._get_gene_and_transcript_version_pk_lookups(genome_build, annotation_consortium)

        release_transcript_version_list = []
        gene_versions_used_by_transcripts = set()
        for transcript_accession, tv_data in cdot_data["transcripts"].items():
            transcript_version_id = transcript_version_ids_by_accession[transcript_accession]
            rtv = ReleaseTranscriptVersion(release=release, transcript_version_id=transcript_version_id)
            release_transcript_version_list.append(rtv)

            if gene_version := tv_data.get("gene_version"):
                gene_versions_used_by_transcripts.add(gene_version)

        release_gene_version_list = []
        for gene_accession in cdot_data["genes"]:

            # Gene accession may not have
            # We only store gene versions that are used in the merged files (which is what's used to insert data)
            if gene_accession in gene_versions_used_by_transcripts:
                gene_version_id = gene_version_ids_by_accession[gene_accession]
                release_gene_version_list.append(ReleaseGeneVersion(release=release, gene_version_id=gene_version_id))

        if release_gene_version_list:
            print(f"Inserting {len(release_gene_version_list)} gene versions for release {release}")
            ReleaseGeneVersion.objects.bulk_create(release_gene_version_list)

        if release_transcript_version_list:
            print(f"Inserting {len(release_transcript_version_list)} transcript versions for release {release}")
            ReleaseTranscriptVersion.objects.bulk_create(release_transcript_version_list)

        print("Matching existing gene list symbols to this release...")
        gm = GeneMatcher(release)
        gm.match_unmatched_in_hgnc_and_gene_lists()

    @staticmethod
    def _get_gene_and_transcript_version_pk_lookups(genome_build: GenomeBuild, annotation_consortium) -> Tuple[Dict, Dict]:
        gene_version_qs = GeneVersion.objects.filter(genome_build=genome_build,
                                                     gene__annotation_consortium=annotation_consortium)
        gene_version_ids_by_accession = {}  # Uses version if non-zero
        for (pk, gene_id, version) in gene_version_qs.values_list("pk", "gene_id", "version"):
            if version:
                gene_accession = f"{gene_id}.{version}"
            else:
                gene_accession = gene_id
            gene_version_ids_by_accession[gene_accession] = pk

        transcript_version_qs = TranscriptVersion.objects.filter(genome_build=genome_build,
                                                                 transcript__annotation_consortium=annotation_consortium)
        tv_values = transcript_version_qs.values_list("pk", "transcript_id", "version")
        transcript_version_ids_by_accession = {f"{transcript_id}.{version}": pk
                                                     for (pk, transcript_id, version) in tv_values}
        return gene_version_ids_by_accession, transcript_version_ids_by_accession

    def _import_merged_data(self, genome_build: GenomeBuild, annotation_consortium, cdot_data: Dict):
        print("_import_merged_data")

        known_uc_gene_symbols = set(GeneSymbol.objects.annotate(uc_symbol=Upper("symbol")).values_list("uc_symbol", flat=True))
        genes_qs = Gene.objects.filter(annotation_consortium=annotation_consortium)
        known_genes_ids = set(genes_qs.values_list("identifier", flat=True))
        transcripts_qs = Transcript.objects.filter(annotation_consortium=annotation_consortium)
        known_transcript_ids = set(transcripts_qs.values_list("identifier", flat=True))
        hgnc_ids = set(HGNC.objects.all().values_list("pk", flat=True))
        gene_version_ids_by_accession, transcript_version_ids_by_accession = self._get_gene_and_transcript_version_pk_lookups(genome_build, annotation_consortium)

        new_gene_symbols = set()
        new_genes = []
        new_gene_versions = []
        modified_gene_versions = []

        def fix_accession(_gene_accession: str) -> str:
            if _gene_accession.startswith("_"):  # Fake UTA
                # If there's only 1 - we can use that
                # otherwise rename slightly to w/ Gene.FAKE_GENE_ID_PREFIX
                _gene_accession = Gene.FAKE_GENE_ID_PREFIX + _gene_accession[1:]
            return _gene_accession

        for gene_accession, gv_data in cdot_data["genes"].items():
            gene_accession = fix_accession(gene_accession)
            gene_id, version = GeneVersion.get_gene_id_and_version(gene_accession)
            if version is None:
                version = 0  # RefSeq genes have no version, store as 0

            if gene_id not in known_genes_ids:
                new_genes.append(Gene(identifier=gene_id,
                                      annotation_consortium=annotation_consortium))
                known_genes_ids.add(gene_id)

            if symbol := gv_data["gene_symbol"]:
                uc_symbol = symbol.upper()
                if uc_symbol not in known_uc_gene_symbols:
                    new_gene_symbols.add(symbol)
                    known_uc_gene_symbols.add(uc_symbol)

            hgnc_id = None  # Foreign Key - only set if we have it (ie not withdrawn)
            if hgnc_identifier := gv_data.get("hgnc"):
                hgnc_identifier = int(hgnc_identifier)
                if hgnc_identifier in hgnc_ids:
                    hgnc_id = hgnc_identifier

            import_source = self._get_import_source_by_url(genome_build, annotation_consortium, gv_data["url"])
            gene_version = GeneVersion(gene_id=gene_id,
                                       version=version,
                                       gene_symbol_id=symbol,
                                       hgnc_identifier=hgnc_identifier,
                                       hgnc_id=hgnc_id,
                                       description=gv_data.get("description"),
                                       biotype=gv_data.get("biotype"),
                                       genome_build=genome_build,
                                       import_source=import_source)

            if pk := gene_version_ids_by_accession.get(gene_accession):
                gene_version.pk = pk
                modified_gene_versions.append(gene_version)
            else:
                new_gene_versions.append(gene_version)

        if new_gene_symbols:
            logging.info("Creating %d new gene symbols", len(new_gene_symbols))
            GeneSymbol.objects.bulk_create([GeneSymbol(symbol=symbol) for symbol in new_gene_symbols],
                                           batch_size=self.BATCH_SIZE)

        if new_genes:
            logging.info("Creating %d new genes", len(new_genes))
            Gene.objects.bulk_create(new_genes, batch_size=self.BATCH_SIZE)

        if new_gene_versions:
            logging.info("Creating %d new gene versions", len(new_gene_versions))
            GeneVersion.objects.bulk_create(new_gene_versions, batch_size=self.BATCH_SIZE)
            # Update with newly inserted records - so that we have a PK to use below
            gene_version_ids_by_accession.update({gv.accession: gv.pk for gv in new_gene_versions})

        # Could potentially be duplicate gene versions (diff transcript versions from diff GFFs w/same GeneVersion)
        if modified_gene_versions:
            logging.info("Updating %d gene versions", len(modified_gene_versions))
            gv_fields = ["gene_symbol_id", "hgnc_identifier", "hgnc_id", "description", "biotype", "import_source"]
            GeneVersion.objects.bulk_update(modified_gene_versions,
                                            gv_fields,
                                            batch_size=self.BATCH_SIZE)

        new_transcript_ids = set()
        new_transcript_versions = []
        modified_transcript_versions = []

        for transcript_accession, tv_data in cdot_data["transcripts"].items():
            transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
            if transcript_id not in known_transcript_ids:
                new_transcript_ids.add(transcript_id)

            build_data = tv_data["genome_builds"][genome_build.name]  # Should always be there as single build file
            gene_accession = fix_accession(tv_data.pop("gene_version"))
            gene_version_id = gene_version_ids_by_accession[gene_accession]
            import_source = self._get_import_source_by_url(genome_build, annotation_consortium, build_data["url"])
            contig = genome_build.chrom_contig_mappings[build_data["contig"]]
            transcript_version = TranscriptVersion(transcript_id=transcript_id,
                                                   version=version,
                                                   gene_version_id=gene_version_id,
                                                   genome_build=genome_build,
                                                   contig=contig,
                                                   import_source=import_source,
                                                   biotype=tv_data.get("biotype"),
                                                   data=tv_data)
            if pk := transcript_version_ids_by_accession.get(transcript_accession):
                transcript_version.pk = pk
                modified_transcript_versions.append(transcript_version)
            else:
                new_transcript_versions.append(transcript_version)

        if new_transcript_ids:
            logging.info("Creating %d new transcripts", len(new_transcript_ids))
            new_transcripts = [Transcript(identifier=transcript_id, annotation_consortium=annotation_consortium)
                               for transcript_id in new_transcript_ids]

            Transcript.objects.bulk_create(new_transcripts, batch_size=self.BATCH_SIZE)
            known_transcript_ids.update(new_transcript_ids)

        # No need to update known after insert as there won't be duplicate transcript versions in the merged data
        if new_transcript_versions:
            logging.info("Creating %d new transcript versions", len(new_transcript_versions))
            TranscriptVersion.objects.bulk_create(new_transcript_versions, batch_size=self.BATCH_SIZE)

        if modified_transcript_versions:
            logging.info("Updating %d transcript versions", len(modified_transcript_versions))
            TranscriptVersion.objects.bulk_update(modified_transcript_versions,
                                                  ["gene_version_id", "import_source", "biotype", "data", "contig"],
                                                  batch_size=self.BATCH_SIZE)

        if new_genes and annotation_consortium == AnnotationConsortium.REFSEQ:
            print("Created new RefSeq genes - retrieving gene summaries via API")
            retrieve_refseq_gene_summaries()
