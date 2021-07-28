from pysam.libcfaidx import FastaFile


class FastaFileContigWrapper:
    """ PyHGVS expects FastaSeqFileDB to return an object
        representing the contig - fake this for pysam """

    def __init__(self, fasta_file, contig):
        self.fasta_file = fasta_file
        self.contig = contig

    def __getitem__(self, _slice):
        return self.fasta_file.fetch(reference=self.contig,
                                     start=_slice.start,
                                     end=_slice.stop)


class FastaFileWrapper:
    """ Implements pygr.SeqFileDB interface
        as pygr (pyghvs dependency) doesn't support Python3 """

    def __init__(self, fasta_filename, convert_chrom_func=None):
        self.fasta_filename = fasta_filename
        self.fasta_file = FastaFile(fasta_filename)
        self.convert_chrom_func = convert_chrom_func

    def get_seq(self, chrom, start, end):
        return self[chrom][start:end]

    def __contains__(self, rname):
        return self.fasta_file.__contains__(rname)

    def __getitem__(self, contig):
        original_contig = contig
        if self.convert_chrom_func:
            contig = self.convert_chrom_func(contig)

        try:
            return FastaFileContigWrapper(self.fasta_file, contig)
        except KeyError:
            if contig != original_contig:
                failed_lookup = f"'{contig}' (originally: '{original_contig}')"
            else:
                failed_lookup = f"'{original_contig}'"
            msg = f"Error loading HGVS: contig {failed_lookup} not in the ref genome '{self.fasta_filename}'"
            raise KeyError(msg)
