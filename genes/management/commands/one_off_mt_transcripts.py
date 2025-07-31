#!/usr/bin/env python3

"""
    MT genes and transcripts, see https://github.com/SACGF/variantgrid_sapath/issues/330

"""
import gzip
import json
import logging

from collections import defaultdict
from django.core.management.base import BaseCommand
from intervaltree import IntervalTree

from annotation.models import VariantAnnotationVersion, VariantAnnotation, VariantTranscriptAnnotation
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from genes.models import GeneAnnotationImport, GeneSymbol, Gene, GeneVersion, TranscriptVersion, HGNC, Transcript
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild, Contig
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv


class Command(BaseCommand):

    def handle(self, *args, **options):
        self._insert_transcripts()
        self._assign_mt_transcripts()

    @staticmethod
    def _insert_transcripts():
        FILES = {
            GenomeBuild.grch37(): "/tau/references/VariantGrid/variantgrid_setup_data/cdot-0.2.28-Homo_sapiens_GRCh37_RefSeq_105.20190906.gff.mt_only.json.gz",
            GenomeBuild.grch38(): "/tau/references/VariantGrid/variantgrid_setup_data/cdot-0.2.28-Homo_sapiens_GRCh38_RefSeq_109.20190607.gff.mt_only.json.gz",
        }

        for genome_build in GenomeBuild.builds_with_annotation():
            print("-" * 40)
            print(f"{genome_build=}")
            cdot_filename = FILES[genome_build]
            data = json.load(gzip.open(cdot_filename))

            gai, _ = GeneAnnotationImport.objects.get_or_create(genome_build=genome_build,
                                                                annotation_consortium=AnnotationConsortium.REFSEQ,
                                                                filename=cdot_filename)
            gene_symbols = []
            genes = []
            gene_versions = []
            for gene_id, gene_data in data["genes"].items():
                gene_symbol_id = gene_data["gene_symbol"]
                gene_symbols.append(GeneSymbol(symbol=gene_symbol_id))
                genes.append(Gene(identifier=gene_id,
                                  annotation_consortium=AnnotationConsortium.REFSEQ,
                                  summary=gene_data["summary"]))

                hgnc = None
                if hgnc_id := gene_data.get("hgnc_id"):
                    hgnc = HGNC.objects.filter(pk=hgnc_id).first()

                gv = GeneVersion(gene_id=gene_id,
                                 version=1, # Always for RefSeq,
                                 gene_symbol_id=gene_symbol_id,
                                 hgnc=hgnc,
                                 biotype=",".join(gene_data["biotype"]),
                                 genome_build=genome_build,
                                 import_source=gai)

                gene_versions.append(gv)


            if gene_symbols:
                print(f"Creating {len(gene_symbols)} gene_symbols")
                GeneSymbol.objects.bulk_create(gene_symbols, ignore_conflicts=True)

            if genes:
                print(f"Creating {len(genes)} genes")
                Gene.objects.bulk_create(genes, ignore_conflicts=True)

            if gene_versions:
                print(f"Creating {len(gene_versions)} gene versions")
                GeneVersion.objects.bulk_create(gene_versions, ignore_conflicts=True)


            gene_versions_by_gene_id = {}
            for gv in GeneVersion.objects.filter(genome_build=genome_build,
                                                 gene__annotation_consortium=AnnotationConsortium.REFSEQ,
                                                 gene_symbol__in=[gs.symbol for gs in gene_symbols]):
                gene_versions_by_gene_id[gv.gene_id] = gv

            transcript_versions = []
            for transcript_id, transcript_data in data["transcripts"].items():
                transcript_id = transcript_id.replace("fake-rna-", "")
                gv = gene_versions_by_gene_id[transcript_data["gene_version"]]
                transcript, _ = Transcript.objects.get_or_create(identifier=transcript_id,
                                                                 annotation_consortium=AnnotationConsortium.REFSEQ)
                td = transcript_data["genome_builds"][genome_build.name]
                exons = [[e[0], e[1]] for e in td["exons"]]
                pyhgvs_data = {
                    "chrom": td["contig"],
                    "start": exons[0][0],
                    "end": exons[-1][1],
                    "exons": exons,
                    "strand": td["strand"],
                    "cds_start": td["cds_start"],
                    "cds_end": td["cds_end"],
                }

                tv = TranscriptVersion(transcript=transcript,
                                       version=1,
                                       gene_version=gv,
                                       genome_build=genome_build,
                                       import_source=gai,
                                       biotype=",".join(transcript_data["biotype"]),
                                       data=pyhgvs_data)
                transcript_versions.append(tv)

            if transcript_versions:
                print(f"Creating {len(transcript_versions)} transcript versions")
                TranscriptVersion.objects.bulk_create(transcript_versions, ignore_conflicts=True)


    @staticmethod
    def _assign_mt_transcripts():
        # Variants are inserted against contig - thus there will only be 1 for both GRCh37/GRCh38
        # However there will be diff VariantAnnotation for each build
        mt_contig_accession = 'NC_012920.1'
        for genome_build in GenomeBuild.builds_with_annotation():
            print("-" * 40)
            print(f"{genome_build=}")

            mt_tx_transcripts = IntervalTree()  # Will only be chrMT
            tv_qs = TranscriptVersion.objects.filter(data__icontains=mt_contig_accession, genome_build=genome_build)
            for tv in tv_qs:
                contig_str = tv.data["chrom"]
                start = tv.data["start"]
                end = tv.data["end"]
                mt_tx_transcripts[start:end] = tv

            vav = VariantAnnotationVersion.latest(genome_build)
            mt_contig = Contig.objects.get(refseq_accession=mt_contig_accession)

            variant_annotation = []  # These will be updated
            variant_transcript_annotation = []  # These will be updated
            variant_gene_overlaps = []

            annotation_classes = {
                VariantAnnotation: variant_annotation,
                VariantTranscriptAnnotation: variant_transcript_annotation,
            }

            for klass, array in annotation_classes.items():
                va_qs = klass.objects.filter(version=vav, variant__locus__contig=mt_contig)

                for va in va_qs.select_related("variant", "variant__locus"):
                    v = va.variant
                    tv_set = mt_tx_transcripts[v.start:v.end]
                    if tv_set:
                        array.append(va)  # Will be modified

                        # Just get the 1st one out
                        tv = next(iter(tv_set)).data
                        va.consequence = "unknown MT exon change"
                        va.gene = tv.gene_version.gene
                        va.transcript = tv.transcript
                        va.transcript_version = tv
                        va.symbol = tv.gene_version.gene_symbol.symbol

                        if klass == VariantAnnotation:
                            # There is no overlapping_symbols column in VG3
                            # va.overlapping_symbols = ",".join(sorted([tv_i.data.gene_version.gene_symbol.symbol for tv_i in tv_set]))

                            # This needs to go in the right partition - if this is run twice it will insert dupes and fail
                            gene = tv.gene_version.gene
                            vgo = {
                                "version_id": va.version_id,
                                "annotation_run_id": va.annotation_run_id,
                                "variant_id": v.pk,
                                "gene_id": gene.pk,
                            }
                            variant_gene_overlaps.append(vgo)

            ava_fields = [
                "consequence",
                "gene",
                "transcript",
                "transcript_version",
                "symbol",
            ]

            for klass, array in annotation_classes.items():
                if array:
                    print(f"Updating {klass}: {len(array)}")
                    klass.objects.bulk_update(array, ava_fields)

            if variant_gene_overlaps:
                # This needs to go in own partition, so use SQL command line
                print(f"Creating {len(variant_gene_overlaps)} variant_gene_overlaps")
                base_table_name = VariantAnnotationVersion.VARIANT_GENE_OVERLAP
                header = [
                    "version_id",
                    "annotation_run_id",
                    "variant_id",
                    "gene_id",
                ]
                DELIMITER = '\t'

                data_filename = f"/tmp/mt_gene_overlap_{genome_build.name}_vav_{vav.pk}.tsv"
                row_data = BulkVEPVCFAnnotationInserter._annotations_list_to_row_data(header, variant_gene_overlaps)
                write_sql_copy_csv(row_data, data_filename, delimiter=DELIMITER)
                partition_table = vav.get_partition_table(base_table_name=base_table_name)

                logging.info("Inserting file '%s' into partition %s", data_filename, partition_table)
                sql_copy_csv(data_filename, partition_table, header, delimiter=DELIMITER)

                # VariantGeneOverlap.objects.bulk_create(variant_gene_overlaps, ignore_conflicts=True)
