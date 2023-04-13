#!/usr/bin/env python3

"""
We inserted gene annotation (from e.g. Universal Transcript Archive) that was missing the Gene ID
We thus inserted it as "unknown_GENE_SYMBOL"

This command finds transcripts in the same build, and if all the genes are the same, switch to that

"""

from django.core.management.base import BaseCommand

from genes.models import Gene, TranscriptVersion


class Command(BaseCommand):

    def handle(self, *args, **options):
        fake_genes = Gene.objects.filter(pk__startswith=Gene.FAKE_GENE_ID_PREFIX)
        print(f"Started with {fake_genes.count()} fake genes")

        num_changed = 0
        num_no_other_transcripts = 0
        num_multiple_gene_versions = 0
        for tv in TranscriptVersion.objects.filter(gene_version__gene__in=fake_genes):
            gene_symbol = tv.gene_version.gene_symbol
            other_transcript_versions = tv.transcript.transcriptversion_set.exclude(pk=tv.pk)
            other_transcript_versions = other_transcript_versions.exclude(gene_version__gene__in=fake_genes)
            other_gene_versions = set()
            for other_tv in other_transcript_versions.filter(genome_build=tv.genome_build,
                                                             gene_version__gene_symbol=gene_symbol):
                other_gene_versions.add(other_tv.gene_version)

            num_gene_versions = len(other_gene_versions)
            if num_gene_versions:
                if num_gene_versions == 1:
                    tv.gene_version = other_gene_versions.pop()
                    tv.save()
                    num_changed += 1
                else:
                    print(f"{tv} had {num_gene_versions} other gene versions. Skipping")
                    num_multiple_gene_versions += 1
            else:
                num_no_other_transcripts += 1

        print(f"Changed: {num_changed}. Skipped - No other trans: {num_no_other_transcripts}, multiple GVs: {num_multiple_gene_versions}")
        Gene.delete_orphaned_fake_genes()
