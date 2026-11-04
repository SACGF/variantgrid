import random
import threading
from pathlib import Path
from unittest.mock import patch

from django.test import SimpleTestCase
from pysam.libcfaidx import FastaFile

from library.genomics.fasta_wrapper import ThreadLocalFastaFile, make_contig_fastas_thread_local

GRCH37_SPARSE_FASTA = str(Path(__file__).parents[2] / "variantgrid/data/reference/sparse_test_fastas/"
                          "GCF_000001405.25_GRCh37.p13_genomic.fna.gz")


class TestThreadLocalFastaFile(SimpleTestCase):
    """ #2624 - a pysam handle read from two threads or across a fork breaks for the life of the process """

    def test_concurrent_threads_fetch_without_errors(self):
        fasta = ThreadLocalFastaFile(GRCH37_SPARSE_FASTA)
        contigs = [c for c in fasta.references if c.startswith("NC_")][:3]
        lengths = {c: fasta.get_reference_length(c) for c in contigs}
        errors = []

        def fetch_randomly():
            for _ in range(500):
                contig = random.choice(contigs)
                start = random.randrange(lengths[contig] - 1)
                try:
                    fasta.fetch(contig, start, start + 1)
                except Exception as e:
                    errors.append(e)

        threads = [threading.Thread(target=fetch_randomly) for _ in range(4)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual([], errors)

    def test_forked_child_opens_its_own(self):
        fasta = ThreadLocalFastaFile(GRCH37_SPARSE_FASTA)
        parent_handle = fasta.fasta_file
        self.assertIs(parent_handle, fasta.fasta_file)
        with patch("library.genomics.fasta_wrapper.os.getpid", return_value=-1):
            self.assertIsNot(parent_handle, fasta.fasta_file)

    def test_make_contig_fastas_thread_local_one_proxy_per_file(self):
        fasta_file = FastaFile(GRCH37_SPARSE_FASTA)
        contig_fastas = dict.fromkeys(fasta_file.references, fasta_file)
        make_contig_fastas_thread_local(contig_fastas)
        proxies = set(map(id, contig_fastas.values()))
        self.assertEqual(1, len(proxies))
        proxy = next(iter(contig_fastas.values()))
        self.assertIsInstance(proxy, ThreadLocalFastaFile)
        self.assertEqual(fasta_file.filename, proxy.filename)
