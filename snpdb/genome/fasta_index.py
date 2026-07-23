import logging

import pandas as pd

from library.utils import file_sha256sum


def load_genome_fasta_index(genome_fasta: 'GenomeFasta', genome_build: 'GenomeBuild'):
    # Delete any older GenomeFasta entries for this build/annotation consortium (eg before we changed settings)
    gf_qs = genome_build.genomefasta_set.filter(annotation_consortium=genome_build.annotation_consortium)
    gf_qs.exclude(pk=genome_fasta.pk).delete()

    index_filename = genome_fasta.filename + ".fai"
    genome_fasta.index_filename = index_filename
    genome_fasta.index_sha256sum = file_sha256sum(index_filename)
    genome_fasta.genome_build = genome_build
    genome_fasta.annotation_consortium = genome_build.annotation_consortium
    genome_fasta.save()

    genome_fasta.genomefastacontig_set.all().delete()

    NAMES = ["NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"]
    df = pd.read_csv(index_filename, sep='\t', header=None, names=NAMES)
    # Can't import GenomeFastaContig at module level (models_genome imports this module)
    GenomeFastaContig = genome_fasta.genomefastacontig_set.model
    genome_fasta_contigs = []
    for _, row in df[["NAME", "LENGTH"]].iterrows():
        name = row["NAME"]
        length = row["LENGTH"]

        contig = genome_build.chrom_contig_mappings.get(name)
        if contig:
            if contig.length != length:
                msg = f"Fasta contig {name} length = {length} but {genome_build} contig {contig} length={contig.length}"
                raise ValueError(msg)

            genome_fasta_contigs.append(GenomeFastaContig(genome_fasta=genome_fasta,
                                                          name=name,
                                                          length=length,
                                                          contig=contig))
        else:
            logging.warning("Could not find contig in %s for fasta name '%s'", genome_build, name)

    GenomeFastaContig.objects.bulk_create(genome_fasta_contigs, batch_size=2000)
    return genome_fasta
