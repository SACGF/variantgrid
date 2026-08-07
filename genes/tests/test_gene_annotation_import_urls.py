from django.test import TestCase

from genes.gene_annotation_import_urls import canonical_import_url


class CanonicalImportUrlTest(TestCase):
    def test_scheme_becomes_https(self):
        self.assertEqual(
            "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz",
            canonical_import_url(
                "ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"))

    def test_ensembl_gff3_becomes_gtf(self):
        self.assertEqual(
            "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz",
            canonical_import_url(
                "ftp://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.gff3.gz"))

    def test_refseq_archive_gff3_is_left_alone(self):
        """ RefSeq archive files are named '.gff3.gz' but have no GTF equivalent, so rewriting the
            extension would invent a url that doesn't exist """
        archive_url = ("ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/"
                       "ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz")
        self.assertEqual(archive_url.replace("ftp://", "https://"), canonical_import_url(archive_url))

    def test_uta_connection_string_untouched(self):
        """ Not a file - a database connection """
        uta = "postgresql://uta.biocommons.org/uta_20241220"
        self.assertEqual(uta, canonical_import_url(uta))
