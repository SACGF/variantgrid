"""
pysam FASTA access: FastaFileWrapper (fasta[chrom][start:end], contig names converted on lookup) and
ThreadLocalFastaFile / make_contig_fastas_thread_local, which give each thread and process its own handle under a
cdot seqfetcher that is shared (cached) across them (#2624).
"""
import os
import threading

from pysam.libcfaidx import FastaFile


class FastaFileContigWrapper:
    """ A single contig from a pysam FastaFile, sliceable by
        genomic coordinate ie contig[start:end] """

    def __init__(self, fasta_file, contig):
        self.fasta_file = fasta_file
        self.contig = contig

    def __getitem__(self, _slice):
        return self.fasta_file.fetch(reference=self.contig,
                                     start=_slice.start,
                                     end=_slice.stop)


class FastaFileWrapper:
    """ A reference genome fasta, indexed by contig name ie fasta[chrom][start:end]
        Contig names can be converted on lookup via convert_chrom_func """

    def __init__(self, fasta_filename, convert_chrom_func=None):
        self.fasta_filename = fasta_filename
        self.fasta_file = FastaFile(fasta_filename)
        self.convert_chrom_func = convert_chrom_func

    def __contains__(self, rname):
        return self.fasta_file.__contains__(rname)

    def __getitem__(self, contig):
        original_contig = contig
        if self.convert_chrom_func:
            contig = self.convert_chrom_func(contig)

        if contig not in self.fasta_file:
            if contig != original_contig:
                failed_lookup = f"'{contig}' (originally: '{original_contig}')"
            else:
                failed_lookup = f"'{original_contig}'"
            msg = f"Error loading HGVS: contig {failed_lookup} not in the ref genome '{self.fasta_filename}'"
            raise KeyError(msg)
        return FastaFileContigWrapper(self.fasta_file, contig)


class ThreadLocalFastaFile:
    """ A FastaFile opened once per thread and per process. One pysam handle read from two threads, or from both sides
        of a fork, breaks for the life of the process ("failure when retrieving sequence" on every fetch) """

    def __init__(self, filename):
        self.filename = filename
        self._local = threading.local()

    @property
    def fasta_file(self) -> FastaFile:
        local = self._local
        # A forked child keeps the forking thread's locals, so the pid says whose handle it is
        if getattr(local, "pid", None) != os.getpid():
            local.fasta_file = FastaFile(self.filename)
            local.pid = os.getpid()
        return local.fasta_file

    def fetch(self, *args, **kwargs) -> str:
        return self.fasta_file.fetch(*args, **kwargs)

    def __getattr__(self, name):
        return getattr(self.fasta_file, name)


def make_contig_fastas_thread_local(contig_fastas: dict):
    """ Swaps the FastaFiles in a cdot GenomeFastaSeqFetcher's contig_fastas ({contig: FastaFile}) for
        ThreadLocalFastaFiles, so the seqfetcher can be cached and shared. Until cdot does this itself (SACGF/cdot#234) """
    proxies = {}
    for contig, fasta_file in contig_fastas.items():
        if (proxy := proxies.get(id(fasta_file))) is None:
            proxy = proxies[id(fasta_file)] = ThreadLocalFastaFile(fasta_file.filename)
            fasta_file.close()
        contig_fastas[contig] = proxy
