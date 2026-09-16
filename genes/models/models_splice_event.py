"""
The names a report gives a recurrent splice junction - AR-V7, EGFRvIII, MET exon 14 skipping.

A splice call arrives as two breakpoints in one gene and becomes a gene-level Variant whose alt
carries the junction's label (@see genes.gene_splice, and snpdb.gene_level_variants for why it is a
Variant at all). This table is what turns those coordinates into the label: it is a naming table,
not a per-variant record, so a junction we have no row for still imports - under a label made from
its own coordinates.

Seeded with the junctions the TSO 500 panel reports; a lab adds rows for the ones its own panel
reports that we have not named.
"""
from django.db import models
from django.db.models.deletion import CASCADE

from genes.models.models_gene import GeneSymbol


class SpliceEvent(models.Model):
    """ A recurrent splice junction and the name a report gives it.

        Identity has two halves, each unique: the coordinates a caller writes (what an import looks
        up), and the gene plus label (what a classification and the alt meet on). The breakpoints
        are build-specific, so a junction is one row per build. """

    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    label = models.TextField()    # what the alt carries: V7, vIII, ex14skip
    display = models.TextField()  # what the report writes: "AR-V7 splice variant"
    genome_build = models.ForeignKey('snpdb.GenomeBuild', on_delete=CASCADE)
    contig = models.ForeignKey('snpdb.Contig', on_delete=CASCADE)
    donor = models.IntegerField()     # Breakpoint 1: last base of the 5' exon
    acceptor = models.IntegerField()  # Breakpoint 2, as the caller writes it

    class Meta:
        unique_together = (("genome_build", "contig", "donor", "acceptor"),
                           ("gene_symbol", "label", "genome_build"))

    def __str__(self):
        return self.display
