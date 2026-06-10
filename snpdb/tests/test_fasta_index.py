from django.test import TestCase

from snpdb.models import GenomeBuild, GenomeFasta


class FastaIndexTest(TestCase):
    def test_load_genome_fasta_index(self):
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        GenomeFasta.objects.all().delete()  # Force get_for_genome_build to rebuild the index
        genome_fasta = GenomeFasta.get_for_genome_build(genome_build)
        contig_qs = genome_fasta.genomefastacontig_set.all()
        self.assertTrue(contig_qs.exists())
        # Every contig must map back to the genome build
        build_contig_ids = set(genome_build.contigs.values_list("pk", flat=True))
        for genome_fasta_contig in contig_qs:
            self.assertIn(genome_fasta_contig.contig_id, build_contig_ids)
