#!/usr/bin/env python3

"""
    MT genes and transcripts, see https://github.com/SACGF/variantgrid_sapath/issues/330

"""
import gzip
import json

from django.core.management.base import BaseCommand

from genes.models import GeneAnnotationImport, GeneSymbol, Gene, GeneVersion, TranscriptVersion, HGNC, Transcript
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild


class Command(BaseCommand):

    def handle(self, *args, **options):
        FILES = {
            GenomeBuild.grch37: "/tau/references/VariantGrid/variantgrid_setup_data/cdot-0.2.28-Homo_sapiens_GRCh37_RefSeq_105.20190906.gff.mt_only.json.gz",
            GenomeBuild.grch38: "/tau/references/VariantGrid/variantgrid_setup_data/cdot-0.2.28-Homo_sapiens_GRCh38_RefSeq_109.20190607.gff.mt_only.json.gz",
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
                gene_symbol = gene_data["gene_symbol"]
                gene_symbols.append(GeneSymbol(symbol=gene_symbol))
                genes.append(Gene(identifier=gene_id,
                                  annotation_consortium=AnnotationConsortium.REFSEQ,
                                  summary=gene_data["summary"]))

                hgnc = None
                if hgnc_id := gene_data.get("hgnc_id"):
                    hgnc = HGNC.objects.filter(pk=hgnc_id).first()

                gv = GeneVersion(gene_id=gene_id,
                                 version=1, # Always for RefSeq,
                                 gene_symbol=gene_symbol,
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
                                                 annotation_consortium=AnnotationConsortium.REFSEQ,
                                                 gene_symbols__in=[gs.symbol for gs in gene_symbols]):
                gene_versions_by_gene_id[gv.gene_id] = gv

            transcript_versions = []
            for transcript_id, transcript_data in data["transcript"].items():
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


