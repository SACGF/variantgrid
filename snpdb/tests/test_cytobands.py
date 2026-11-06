from django.test import TestCase

from snpdb.cytobands import CytobandRange, parse_cytoband
from snpdb.models import GenomeBuild


class ParseCytobandTest(TestCase):
    def test_parse_cytoband(self):
        for text, expected in [
            ("11p13", CytobandRange("11", "p13", "p13")),
            ("chrXq28", CytobandRange("X", "q28", "q28")),
            ("17P", CytobandRange("17", "p", "p")),
            ("11p15.5 - p15.4", CytobandRange("11", "p15.5", "p15.4")),
            ("11p11.2-11q12", CytobandRange("11", "p11.2", "q12")),
            ("11p13-12q12", None),  # across chromosomes
            ("21", None),
            ("chr1:100-200", None),
        ]:
            with self.subTest(text=text):
                self.assertEqual(parse_cytoband(text), expected)


class CytobandIntervalTest(TestCase):
    """ Against the cytoband files loaded by the migration - coordinates as in the file, start made 1-based """

    def _interval(self, genome_build, text):
        if interval := parse_cytoband(text).get_interval(genome_build):
            return interval.contig.name, interval.start, interval.end
        return None

    def test_get_interval(self):
        grch37 = GenomeBuild.grch37()
        grch38 = GenomeBuild.grch38()
        for genome_build, text, expected in [
            (grch37, "11p13", ("11", 31000001, 36400000)),
            (grch37, "11p15", ("11", 1, 21700000)),  # every p15.x sub-band
            (grch37, "11p15.5-p15.4", ("11", 1, 10700000)),
            (grch37, "11p12-p13", ("11", 31000001, 43500000)),  # either order
            (grch37, "17p", ("17", 1, 24000000)),
            (grch38, "11p15.5-p15.4", ("11", 1, 11700000)),
            (grch37, "11p99", None),
            (grch37, "11p13-p99", None),
        ]:
            with self.subTest(genome_build=genome_build.name, text=text):
                self.assertEqual(self._interval(genome_build, text), expected)
