"""
Biocommons hgvs Babelfish over a genome build's reference FASTA, for the VCF-coordinate work that needs no
transcripts: telling whether an explicit allele is a del/dup/inv (vcf_coordinate_as_symbolic) and expanding a
symbolic one (vcf_coordinate_as_explicit). Entry points: get_genome_babelfish, babelfish_assembly_name.
The HGVS converter builds its own Babelfish over the transcript data provider, with the same assembly name.
"""
from functools import cache

from cdot.hgvs.dataproviders.fasta_seqfetcher import GenomeFastaSeqFetcher
from hgvs.extras.babelfish import Babelfish

from snpdb.models.models_genome import GenomeBuild


def babelfish_assembly_name(genome_build: GenomeBuild) -> str:
    # GRCh37 needs the patch name to get the MT chromosome mapping
    if genome_build.name == 'GRCh37':
        return genome_build.get_build_with_patch()
    return genome_build.name


class GenomeFastaSequenceProvider(GenomeFastaSeqFetcher):
    """ The one data provider method Babelfish's VCF coordinate conversions call. Looks contigs up by
        accession, as the HGVS converter's fasta seqfetcher does """

    def get_seq(self, ac, start_i=None, end_i=None) -> str:
        return self.fetch_seq(ac, start_i, end_i)


@cache
def get_genome_babelfish(genome_build: GenomeBuild) -> Babelfish:
    """ No symbolic_alt_min_length - callers pass the length rule they need """
    return Babelfish(GenomeFastaSequenceProvider(genome_build.reference_fasta), babelfish_assembly_name(genome_build))
