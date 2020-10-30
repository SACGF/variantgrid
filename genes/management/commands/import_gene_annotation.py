"""

See [annotation] page for instructions

GenePred is way easier to parse, but also want GFF3 as then we can have gene ID, biotype and description etc

"""
from collections import Counter, defaultdict
from django.core.management.base import BaseCommand
from pyhgvs.utils import read_genepred
import os
import re

from genes.gene_matching import GeneMatcher
from genes.models import GeneAnnotationImport, HGNCGeneNames, \
    GeneSymbol, Gene, GeneVersion, Transcript, TranscriptVersion, GeneAnnotationRelease, ReleaseGeneVersion, \
    ReleaseTranscriptVersion
from genes.models_enums import AnnotationConsortium
from library.django_utils import highest_pk
from library.file_utils import open_handle_gzip
from library.utils import invert_dict
from snpdb.models.models_genome import GenomeBuild, GenomeFasta
from snpdb.models.models_enums import SequenceRole


class Command(BaseCommand):
    BATCH_SIZE = 2000
    FAKE_GENE_ID_PREFIX = "unknown_"  # Legacy from when we allowed inserting GenePred w/o GFF3

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.annotation_consortium = None
        self.genome_build = None
        self.contig_id_to_fasta = None
        self.hgnc_ids = set(HGNCGeneNames.objects.values_list("pk", flat=True))
        # Known objects containers are updated with new inserts
        self.known_gene_symbols = set(GeneSymbol.objects.all().values_list("pk", flat=True))
        self.known_gene_versions_by_gene_id = defaultdict(dict)
        self.known_transcript_versions_by_transcript_id = defaultdict(dict)

    def add_arguments(self, parser):
        consortia = [ac[1] for ac in AnnotationConsortium.CHOICES]
        builds = [gb.name for gb in GenomeBuild.builds_with_annotation()]

        parser.add_argument('--genome-build', choices=builds, required=True)
        parser.add_argument('--annotation-consortium', choices=consortia, required=True)
        parser.add_argument('--replace', action='store_true', help="Replace gene symbols and relations")
        parser.add_argument('--release', type=int, required=False,
                            help="Make a release (to match VEP) store all gene/transcript versions")

        gff3_help = 'Use "None" if using only genePred (not recommended)'
        parser.add_argument('--gff3', nargs="+", required=True, help=gff3_help)
        parser.add_argument('--genePred', nargs="+", required=True)

    def handle(self, *args, **options):
        build_name = options["genome_build"]
        annotation_consortium_name = options["annotation_consortium"]
        replace = options["replace"]
        release_version = options["release"]
        gff3_filenames = options["gff3"]
        genepred_filenames = options["genePred"]

        genome_build = GenomeBuild.get_name_or_alias(build_name)
        ac_dict = invert_dict(dict(AnnotationConsortium.CHOICES))
        annotation_consortium = ac_dict[annotation_consortium_name]

        # gff/genePred sanity checks
        if len(gff3_filenames) != len(genepred_filenames):
            raise ValueError("Must supply exactly the same amount of gff3 and genePred files!")
        if len(gff3_filenames) != 1:
            if release_version:
                raise ValueError("You can only specify 1 GFF in a release!")
            if replace:
                raise ValueError("You can only specify 1 GFF when using replace!")

        for gff3_filename, genepred_filename in zip(gff3_filenames, genepred_filenames):
            if not os.path.exists(genepred_filename):
                raise FileNotFoundError(genepred_filename)
            name, extension = os.path.splitext(genepred_filename)
            if extension != '.genePred':
                raise ValueError(f"genePred files must have .genePred extension, was: '{extension}'")
            if not gff3_filename.startswith(name):
                raise ValueError(f"Start of GFF3 file '{gff3_filename}' must match start of genePred '{name}'")
            if not os.path.exists(gff3_filename):
                raise FileNotFoundError(gff3_filename)

        self.annotation_consortium = annotation_consortium
        self.genome_build = genome_build
        genome_fasta = GenomeFasta.get_for_genome_build(genome_build)
        self.contig_id_to_fasta = genome_fasta.get_contig_id_to_name_mappings()

        # Known objects containers are updated with new inserts
        gv_qs = GeneVersion.objects.filter(gene__annotation_consortium=annotation_consortium, genome_build=genome_build)
        self.update_known_gene_versions_by_gene_id(gv_qs)

        tv_qs = TranscriptVersion.objects.filter(transcript__annotation_consortium=annotation_consortium, genome_build=genome_build)
        self.update_known_transcript_versions_by_transcript_id(tv_qs)

        num_remaining_file_pairs = len(gff3_filenames)
        for gff3_filename, genepred_filename in zip(gff3_filenames, genepred_filenames):
            num_remaining_file_pairs -= 1
            update_known_objects = num_remaining_file_pairs > 0
            self.insert_gene_annotations(gff3_filename, genepred_filename, update_known_objects,
                                         replace=replace, release_version=release_version)

        # Remove orphaned fake genes
        used_genes = TranscriptVersion.objects.filter(gene_version__gene__identifier__startswith=self.FAKE_GENE_ID_PREFIX).values_list("gene_version__gene")
        qs = Gene.objects.filter(identifier__startswith=self.FAKE_GENE_ID_PREFIX).exclude(identifier__in=used_genes)
        ret = qs.delete()
        if ret:
             print(f"Deleted orphaned {self.FAKE_GENE_ID_PREFIX} records:")
             print(ret)

    def update_known_gene_versions_by_gene_id(self, gv_qs):
        for gv in gv_qs:
            self.known_gene_versions_by_gene_id[gv.gene_id][gv.version] = gv

    def update_known_transcript_versions_by_transcript_id(self, tv_qs):
        for tv in tv_qs:
            self.known_transcript_versions_by_transcript_id[tv.transcript_id][tv.version] = tv

    @staticmethod
    def get_gene_accession(gene_id, version):
        """ GeneVersion.accession doesn't display version for RefSeq """
        return f"{gene_id}.{version}"

    def insert_gene_annotations(self, gff3_filename, genepred_filename, update_known_objects,
                                replace=False, release_version=None):
        """ update_known_objects - set false as optimisation if you know it's the last loop """

        print(f"insert_gene_annotations('{gff3_filename}', '{genepred_filename}')")
        import_source = GeneAnnotationImport.objects.create(genome_build=self.genome_build,
                                                            annotation_consortium=self.annotation_consortium,
                                                            filename=gff3_filename)
        release = None
        if release_version:
            release, created = GeneAnnotationRelease.objects.update_or_create(version=release_version,
                                                                              genome_build=self.genome_build,
                                                                              annotation_consortium=self.annotation_consortium,
                                                                              defaults={
                                                                                  "gene_annotation_import": import_source
                                                                              })
            if not created:
                print("Release exists - clearing existing data")
                release.releasegeneversion_set.all().delete()
                release.releasetranscriptversion_set.all().delete()

        parser = gff_parser_factory(self.genome_build, import_source, self.hgnc_ids)

        # Ensembl geneProd doesn't have version, so we have to store what was used in the GFF3
        gff3_gene_versions = {}
        gff3_transcript_versions = {}
        gff3_existing_transcript_versions = TranscriptVersionContainer()
        unknown_gene_symbols = []
        unknown_gene_ids = []
        unknown_gene_versions = []
        unknown_transcripts = []
        unknown_transcript_versions_by_gene_accession = defaultdict(list)
        # When you have genePred only - we don't know the GeneID, so create fake ones with "unknown_"
        # If we import the same transcript ID later via a GFF3, we can change the gene version to the real one
        known_transcript_versions_to_update_with_gene_accession = defaultdict(list)
        known_gene_versions_to_update = []

        with open_handle_gzip(gff3_filename, "rt") as f:
            for line in f:
                try:
                    if line.startswith("#"):
                        continue
                    cols = line.strip().split("\t")
                    feature_type = cols[2]
                    if feature_type in {"CDS", "exon"}:
                        continue

                    attributes = dict((a.split("=") for a in cols[8].split(";")))
                    if parser.is_gene(feature_type, attributes):  # GeneVersion
                        gene_version = parser.get_gene_version(attributes)
                        gene_id = gene_version.gene_id
                        gff3_gene_versions[gene_id] = gene_version.version

                        if gene_version.gene_symbol_id not in self.known_gene_symbols:
                            unknown_gene_symbols.append(GeneSymbol(symbol=gene_version.gene_symbol_id))

                        gene_versions_dict = self.known_gene_versions_by_gene_id.get(gene_id, {})
                        if not gene_versions_dict:
                            unknown_gene_ids.append(gene_id)

                        known_gene_version = gene_versions_dict.get(gene_version.version)
                        if known_gene_version:
                            if replace:
                                if known_gene_version.gene_symbol_id != gene_version.gene_symbol_id:
                                    known_gene_version.import_source = import_source
                                    known_gene_version.gene_symbol_id = gene_version.gene_symbol_id
                                    known_gene_versions_to_update.append(known_gene_version)
                        else:
                            unknown_gene_versions.append(gene_version)

                    elif parser.is_transcript(feature_type, attributes):
                        gene_id, transcript_version = parser.get_gene_id_and_transcript(attributes)
                        transcript_id = transcript_version.transcript_id
                        gff3_transcript_versions[transcript_id] = transcript_version.version
                        transcript_versions_dict = self.known_transcript_versions_by_transcript_id.get(transcript_id, {})
                        if not transcript_versions_dict:
                            unknown_transcripts.append(Transcript(identifier=transcript_id,
                                                                  annotation_consortium=self.annotation_consortium))

                        gene_accession = self.get_gene_accession(gene_id, gff3_gene_versions[gene_id])
                        if known_transcript_version := transcript_versions_dict.get(transcript_version.version):
                            gff3_existing_transcript_versions.add(known_transcript_version)

                            # Always replace if starts with "unknown_" (or replace and different)
                            current_gene_accession = known_transcript_version.gene_version.accession
                            if current_gene_accession.startswith(self.FAKE_GENE_ID_PREFIX) or \
                                    (replace and current_gene_accession != gene_accession):
                                known_transcript_versions_to_update_with_gene_accession[gene_accession].append(known_transcript_version)
                                # print(f"Updating {known_transcript_version} => {gene_accession}")
                        else:
                            # print(f"{transcript_version.accession} is unknown - have {transcript_versions_dict}")
                            unknown_transcript_versions_by_gene_accession[gene_accession].append(transcript_version)
                except:
                    print(f"Couldn't handle line:")
                    print(line)
                    raise

        genepred_data_by_transcript_id = self.get_genepred_data_by_transcript_id(parser, gff3_existing_transcript_versions, genepred_filename)

        # Insert new data
        if unknown_gene_symbols:
            print(f"Inserting {len(unknown_gene_symbols)} gene symbols")
            GeneSymbol.objects.bulk_create(unknown_gene_symbols, batch_size=Command.BATCH_SIZE, ignore_conflicts=True)
            self.known_gene_symbols.update(unknown_gene_symbols)

        if unknown_gene_ids:
            print(f"Inserting {len(unknown_gene_ids)} genes")
            unknown_genes = (Gene(identifier=gene_id, annotation_consortium=self.annotation_consortium) for gene_id in unknown_gene_ids)
            Gene.objects.bulk_create(unknown_genes, batch_size=Command.BATCH_SIZE, ignore_conflicts=True)

        if unknown_gene_versions:
            print(f"Inserting {len(unknown_gene_versions)} gene versions")
            old_max_gene_id = highest_pk(GeneVersion)
            GeneVersion.objects.bulk_create(unknown_gene_versions, batch_size=Command.BATCH_SIZE, ignore_conflicts=True)
            gene_version_qs = GeneVersion.objects.filter(gene__annotation_consortium=self.annotation_consortium,
                                                         genome_build=self.genome_build)
            if old_max_gene_id:
                gene_version_qs = gene_version_qs.filter(pk__gt=old_max_gene_id)
            self.update_known_gene_versions_by_gene_id(gene_version_qs)

        if unknown_transcripts:
            print(f"Inserting {len(unknown_transcripts)} transcripts")
            Transcript.objects.bulk_create(unknown_transcripts, batch_size=Command.BATCH_SIZE, ignore_conflicts=True)

        unknown_transcript_versions = []
        unknown_transcript_ids = set()
        for gene_accession, transcript_versions in unknown_transcript_versions_by_gene_accession.items():
            gene_id, version = gene_accession.rsplit(".", 1)
            gene_versions_dict = self.known_gene_versions_by_gene_id[gene_id]
            try:
                gene_version = gene_versions_dict[int(version)]
            except KeyError:
                print(f"gene_accession: {gene_accession} No {version} - had {gene_versions_dict}")
                raise

            for tv in transcript_versions:
                tv.gene_version = gene_version
                # Add genePred data now so it can be inserted via bulk_create
                # Since there is only 1 version of a transcript in a file, matching by ID (w/o version) is ok
                unknown_transcript_ids.add(tv.transcript_id)
                data = genepred_data_by_transcript_id.get(tv.transcript_id, {})
                self.set_transcript_data(tv, data)
                unknown_transcript_versions.append(tv)

        if unknown_transcript_versions:
            print(f"Inserting {len(unknown_transcript_versions)} transcript versions")
            TranscriptVersion.objects.bulk_create(unknown_transcript_versions, batch_size=Command.BATCH_SIZE, ignore_conflicts=True)

        # UPDATES
        if known_gene_versions_to_update:
            print(f"Updating {len(known_gene_versions_to_update)} GeneVersions to new symbols")
            GeneVersion.objects.bulk_update(known_gene_versions_to_update,
                                            fields=["import_source", "gene_symbol_id"], batch_size=Command.BATCH_SIZE)
            # No need to update gene versions as only ever 1 GFF when using "replace-symbols"

        # Remove any genePred that were inserted in unknown_transcript_versions
        # print("unknown_transcript_ids:")
        # print(unknown_transcript_ids)
        genepred_data_by_transcript_id = {k: v for k, v in genepred_data_by_transcript_id.items() if k not in unknown_transcript_ids}

        transcript_versions_to_update = []
        for transcript_name, data in genepred_data_by_transcript_id.items():
            transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_name)
            existing_transcript_version = gff3_existing_transcript_versions.get(transcript_id, version)
            if existing_transcript_version is None:
                raise ValueError(f"Could not obtain {transcript_name} from GFF data!")
            existing_transcript_version.import_source = import_source
            self.set_transcript_data(existing_transcript_version, data)
            transcript_versions_to_update.append(existing_transcript_version)

        if transcript_versions_to_update:
            print(f"Updating {len(transcript_versions_to_update)} transcript versions")
            fields = ['import_source', 'data']
            TranscriptVersion.objects.bulk_update(transcript_versions_to_update, fields, batch_size=Command.BATCH_SIZE)

        if known_transcript_versions_to_update_with_gene_accession:
            transcript_versions_to_update = []
            for gene_accession, transcript_versions in known_transcript_versions_to_update_with_gene_accession.items():
                gene_id, version = gene_accession.rsplit(".", 1)
                gene_versions_dict = self.known_gene_versions_by_gene_id[gene_id]
                try:
                    gene_version = gene_versions_dict[int(version)]
                except KeyError:
                    print(f"gene_accession: {gene_accession} No {version} - had {gene_versions_dict}")
                    raise

                for tv in transcript_versions:
                    tv.import_source = import_source
                    tv.gene_version = gene_version
                    transcript_versions_to_update.append(tv)

            print(f"Updating {len(transcript_versions_to_update)} TranscriptVersions with proper gene version")
            TranscriptVersion.objects.bulk_update(transcript_versions_to_update,
                                                  ["import_source", "gene_version"], batch_size=Command.BATCH_SIZE)

        if release or (update_known_objects and (unknown_transcript_versions or transcript_versions_to_update)):
            tv_qs = import_source.transcriptversion_set.all()
            self.update_known_transcript_versions_by_transcript_id(tv_qs)

        if release:
            release_gene_version_list = []
            for gene_id, version in gff3_gene_versions.items():
                gv = self.known_gene_versions_by_gene_id[gene_id][version]
                release_gene_version_list.append(ReleaseGeneVersion(release=release, gene_version=gv))
            if release_gene_version_list:
                print(f"Inserting {len(release_gene_version_list)} gene versions for release {release}")
                ReleaseGeneVersion.objects.bulk_create(release_gene_version_list)

            release_transcript_version_list = []
            for transcript_id, version in gff3_transcript_versions.items():
                tv = self.known_transcript_versions_by_transcript_id[transcript_id][version]
                release_transcript_version_list.append(ReleaseTranscriptVersion(release=release, transcript_version=tv))
            if release_transcript_version_list:
                print(f"Inserting {len(release_transcript_version_list)} transcript versions for release {release}")
                ReleaseTranscriptVersion.objects.bulk_create(release_transcript_version_list)

            print("Matching existing gene list symbols to this release...")
            gm = GeneMatcher(release)
            gm.match_unmatched_in_gene_lists()

    @staticmethod
    def set_transcript_data(existing_transcript_version, data):
        data["id"] = existing_transcript_version.accession  # Ensembl genePred data does not contain versions
        existing_transcript_version.data = data

    def get_genepred_data_by_transcript_id(self, parser, gff3_existing_transcripts, genepred_filename):
        """ Ensembl genePred:

            command       gff3ToGenePred Homo_sapiens.GRCh38.97.gff3.gz Homo_sapiens.GRCh38.97.genePred
            transcript    transcript:ENST00000440066
            gene          gene:ENSG00000205916

            RefSeq genePred:

            command        gff3ToGenePred -geneNameAttr=Dbxref -rnaNameAttr=transcript_id ref_GRCh38.p12_top_level.gff3.gz ref_GRCh38.p12_top_level.genePred
            transcript    NM_153443.4
            gene        GeneID:115653
        """
        genepred_data_by_transcript_id = {}

        skipped_chroms = Counter()
        num_existing_skipped = 0
        num_existing_replaced_standard_contigs = 0
        num_existing_replaced_y_pseudoautosomal = 0
        num_no_contig_match = 0
        with open(genepred_filename) as f:
            for transcript_json in read_genepred(f):
                # No longer store gene_name in transcript JSON (will use TranscriptVersion/GeneVersion relation)
                del transcript_json["gene_name"]
                # NCBI file uses refseq accession IDs
                # Map everything to our contigs, then convert to our Fasta
                chrom = transcript_json['chrom']
                contig = self.genome_build.chrom_contig_mappings.get(chrom)
                if contig:
                    fasta_chrom = self.contig_id_to_fasta.get(contig.pk)
                else:
                    fasta_chrom = None

                if fasta_chrom is None:
                    num_no_contig_match += 1
                    skipped_chroms[chrom] += 1
                    continue
                transcript_json['chrom'] = fasta_chrom
                try:
                    transcript_json['id'] = parser.clean_and_validate_genepred_id(transcript_json['id'])
                except:
                    continue  # Skip

                # We only have 1 versioned transcript per build (unique_together constraint)
                # as this simplifies resolving HGVS to variant loci, but because there can be
                # multiple we choose based on contig (ie standard chroms over alt/fixes etc)
                transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_json['id'])
                existing_data = genepred_data_by_transcript_id.get(transcript_id)
                if not existing_data:
                    existing = gff3_existing_transcripts.get(transcript_id, version)
                    if existing and existing.data:
                        existing_data = existing.data

                if existing_data:
                    existing_chrom = existing_data.get("chrom")
                    if existing_chrom:
                        existing_contig = self.genome_build.chrom_contig_mappings[existing_chrom]
                        if existing_contig.role != contig.role and contig.role == SequenceRole.ASSEMBLED_MOLECULE:
                            num_existing_replaced_standard_contigs += 1
                            # print(f"Replaced {transcript_name} was: {existing_contig} now {contig}")
                        elif existing_contig.name == 'Y' and contig.name == 'X':
                            # issue #2047 - Y pseudoautosomal region - take X over existing Y
                            num_existing_replaced_y_pseudoautosomal += 1
                        else:
                            num_existing_skipped += 1
                            continue  # skip (don't update existing)

                genepred_data_by_transcript_id[transcript_id] = transcript_json

        print("=" * 40)
        print(f"genepred_filename: {genepred_filename}")
        print(f"num_existing_skipped = {num_existing_skipped}")
        print(f"num_existing_replaced_standard_contigs = {num_existing_replaced_standard_contigs}")
        print(f"num_existing_replaced_y_pseudoautosomal = {num_existing_replaced_y_pseudoautosomal}")
        print(f"num_no_contig_match = {num_no_contig_match}")
        return genepred_data_by_transcript_id


class TranscriptVersionContainer:
    """ Store the transcripts used in a GFF file (retrieving by version if available) """

    def __init__(self):
        self.transcript_version_by_accession = {}
        self.transcript_version_by_id = {}

    def add(self, transcript_version):
        self.transcript_version_by_accession[transcript_version.accession] = transcript_version
        self.transcript_version_by_id[transcript_version.transcript_id] = transcript_version

    def get(self, transcript_id, version):
        if version:
            accession = f"{transcript_id}.{version}"
            transcript = self.transcript_version_by_accession.get(accession)
        else:
            transcript = self.transcript_version_by_id.get(transcript_id)
        return transcript


def gff_parser_factory(genome_build: GenomeBuild, import_source, *args, **kwargs):
    ac = import_source.annotation_consortium
    if ac == AnnotationConsortium.REFSEQ:
        return RefSeqParser(genome_build, import_source, *args, **kwargs)
    elif ac == AnnotationConsortium.ENSEMBL:
        return EnsemblParser(genome_build, import_source, *args, **kwargs)


class GFFParser:

    def __init__(self, genome_build: GenomeBuild, import_source, hgnc_ids):
        self.genome_build = genome_build
        self.import_source = import_source
        self.hgnc_ids = hgnc_ids

    def is_transcript(self, feature_type, attributes):
        return "transcript_id" in attributes

    def clean_and_validate_genepred_id(self, transcript):
        return transcript


class RefSeqParser(GFFParser):
    GENE_TYPES = {"gene", "pseudogene"}

    @staticmethod
    def get_dbxref(attributes: dict) -> dict:
        dbxref = {}
        dbxref_str = attributes.get("Dbxref")
        if dbxref_str:
            dbxref = dict(d.split(":", 1) for d in dbxref_str.split(","))
        return dbxref

    def is_gene(self, feature_type, attributes):
        return feature_type in RefSeqParser.GENE_TYPES
        #dbxref = self.get_dbxref(attributes)
        #return dbxref.get("GeneID")

    def is_transcript(self, feature_type, attributes):
        valid = super().is_transcript(feature_type, attributes)
        if valid:
            transcript_id = attributes["transcript_id"]
            valid = RefSeqParser.valid_transcript_id(transcript_id)
        return valid

    def get_gene_version(self, attributes):
        """ Looks like:
            {'ID': 'gene10861',
             'Dbxref': 'GeneID:2624,HGNC:HGNC:4171,MIM:137295',
             'Name': 'GATA2',
             'description': 'GATA binding protein 2',
             'gbkey': 'GeneVersion',
             'gene': 'GATA2',
             'gene_biotype': 'protein_coding',
             'gene_synonym': 'DCML,IMD21,MONOMAC,NFE1B'} """

        dbxref = self.get_dbxref(attributes)
        gene_id = dbxref["GeneID"]
        hgnc_id = None
        hgnc = dbxref.get("HGNC")
        if hgnc:
            hgnc = hgnc.replace("HGNC:", "")
            if hgnc in self.hgnc_ids:
                hgnc_id = hgnc

        version_1 = 1  # RefSeq gene IDs don't have a version - always be 1 for unique_together constraint
        return GeneVersion(gene_id=gene_id,
                           version=version_1,
                           gene_symbol_id=attributes["gene"],
                           hgnc_id=hgnc_id,
                           description=attributes.get("description"),
                           biotype=attributes.get("gene_biotype"),  # Only in GRCh38
                           genome_build=self.genome_build,
                           import_source=self.import_source)

    def get_gene_id_and_transcript(self, attributes):
        """ Looks like:
            {'ID': 'rna34632',
             'Parent': 'gene10861',
             'Dbxref': 'GeneID:2624,Genbank:NM_032638.4,HGNC:HGNC:4171,MIM:137295',
             'Name': 'NM_032638.4',
             'gbkey': 'mRNA',
             'gene': 'GATA2',
             'product': 'GATA binding protein 2%2C transcript variant 2',
             'transcript_id': 'NM_032638.4'} """
        dbxref = self.get_dbxref(attributes)
        gene_id = dbxref["GeneID"]
        transcript_name = attributes["transcript_id"]
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_name)

        transcript = TranscriptVersion(transcript_id=transcript_id,
                                       version=version,
                                       genome_build=self.genome_build,
                                       import_source=self.import_source,
                                       biotype=None)  # No RefSeq transcript biotype

        return gene_id, transcript

    @staticmethod
    def valid_transcript_id(transcript_id):
        VALID_TRANSCRIPT_TYPES = ["NM_", "NR_"]
        return any(map(transcript_id.startswith, VALID_TRANSCRIPT_TYPES))

    def clean_and_validate_genepred_id(self, transcript):
        if not self.valid_transcript_id(transcript):
            raise RuntimeError
        return transcript


class EnsemblParser(GFFParser):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hgnc_pattern = re.compile(r"(.*) \[Source:HGNC.*Acc:HGNC:(\d+)\]")

    def is_gene(self, feature_type, attributes):
        return 'gene_id' in attributes

    def get_gene_version(self, attributes):
        """ Looks like:
            {'ID': 'gene:ENSG00000179348',
             'Name': 'GATA2',
             'biotype': 'protein_coding',
             'description': 'GATA binding protein 2 [Source:HGNC Symbol%3BAcc:HGNC:4171]',
             'gene_id': 'ENSG00000179348',
             'logic_name': 'ensembl_havana_gene',
             'version': '11'} """

        hgnc_id = None
        description = attributes.get("description")
        if description:
            if m := self.hgnc_pattern.match(description):
                description, potential_hgnc_id = m.groups()
                if potential_hgnc_id in self.hgnc_ids:
                    hgnc_id = potential_hgnc_id

        return GeneVersion(gene_id=attributes["gene_id"],
                           version=int(attributes["version"]),
                           gene_symbol_id=attributes["Name"],
                           hgnc_id=hgnc_id,
                           description=description,
                           biotype=attributes["biotype"],
                           genome_build=self.genome_build,
                           import_source=self.import_source)

    def get_gene_id_and_transcript(self, attributes):
        """ Looks like:
            {'ID': 'transcript:ENST00000468377',
             'Parent': 'gene:ENSG00000244300',
             'Name': 'GATA2-AS1-202',
             'biotype': 'lncRNA',
             'transcript_id': 'ENST00000468377',
             'transcript_support_level': '2',
             'version': '1'} """
        gene_id = attributes["Parent"]
        if gene_id.startswith("gene:"):
            gene_id = gene_id[5:]
        if not gene_id.startswith("ENSG"):
            raise ValueError(f"Expected gene_id to start with ENSG was: '{gene_id}'")

        transcript = TranscriptVersion(transcript_id=attributes["transcript_id"],
                                       version=int(attributes["version"]),
                                       genome_build=self.genome_build,
                                       import_source=self.import_source,
                                       biotype=attributes["biotype"])
        return gene_id, transcript

    def clean_and_validate_genepred_id(self, transcript):
        if transcript.startswith("transcript:"):
            transcript = transcript[11:]
        else:
            raise ValueError("not a valid Ensembl genePred ID")
        return transcript
