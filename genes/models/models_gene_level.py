"""
The gene identity every gene-level Variant is built on, and the whole-gene copy number event.

@see snpdb.gene_level_variants for why a gene-level event is stored as a Variant at all, and
genes.models.models_gene_fusion for the other event that uses these ids.
"""
from typing import Optional

from django.core.exceptions import ValidationError
from django.db import IntegrityError, models, transaction
from django.db.models import TextField
from django.db.models.deletion import CASCADE, PROTECT, SET_NULL

from genes.models.models_gene import HGNC, Gene, GeneSymbol
from library.genomics.vcf_enums import GeneIdNamespace, GeneLevelSymbolicAlt


class GeneLevelId(models.Model):
    """ A stable number for one gene of a gene-level event - the pk is what the Variant carries as its
        Locus.position and inside its symbolic alt (@see GeneLevelSymbolicAlt), so this table is the
        identifier space gene-level identity is built on. @see snpdb.gene_level_variants for why such
        an event is stored as a Variant in the first place.

        pk is the HGNC ID where the gene has one, so the identity of an ordinary event is the same
        number on every deployment. Symbols HGNC doesn't carry - clone-based identifiers like
        RP11-458D21.5, which turn up as fusion partners routinely - get a pk allocated above
        CUSTOM_ID_START, so every call the caller made can still become a Variant.

        Custom numbers are local to this deployment, which is why the alt namespaces them as GENE
        rather than HGNC. Anything leaving the system - display, export, a classification sent to
        another instance - uses symbol_str rather than the pk, and the receiver resolves it back
        through its own table. @see GeneLevelSymbolicAlt for where the numbers are used, and
        GeneFusion.canonical_str / GeneCopyNumberEvent.canonical_str for the string forms.

        Identity, once handed out, is fixed - the same way a Variant's coordinates are - so a name
        that only later becomes resolvable keeps the custom pk its variants were created with. What
        moves it onto the HGNC number is re-loading the caller's file: resolution mints the right
        identity from the start, and the old rows stay as they are. """

    CUSTOM_ID_START = 1_000_000
    CUSTOM_ID_RETRIES = 5

    # The name to use whenever this leaves the system - the approved symbol where we have one,
    # otherwise the name exactly as the caller wrote it
    symbol_str = TextField(unique=True, db_collation='case_insensitive')
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=SET_NULL)
    hgnc = models.ForeignKey(HGNC, null=True, on_delete=SET_NULL)
    # What the caller's breakpoint landed in - the same gene as an Entrez id and an ENSG, since a
    # release is one consortium's. Annotation reads these before falling back to the symbol, which
    # is what makes a side found by position reachable from a gene list @see gene_level_annotation
    genes = models.ManyToManyField(Gene, blank=True)

    def __str__(self):
        return self.symbol_str

    @property
    def is_custom(self) -> bool:
        return self.pk >= GeneLevelId.CUSTOM_ID_START

    @property
    def alt_namespace(self) -> str:
        return GeneIdNamespace.GENE if self.is_custom else GeneIdNamespace.HGNC

    @staticmethod
    def get_or_create_for_symbol(symbol_str: str, gene_symbol_id: Optional[str],
                                 hgnc: Optional[HGNC]) -> 'GeneLevelId':
        """ symbol_str is what to show and send - callers resolve aliases first, so SEPT14 arrives
            here as SEPTIN14 and collapses onto the one row """

        if hgnc is not None:
            defaults = {"symbol_str": symbol_str, "gene_symbol_id": gene_symbol_id, "hgnc": hgnc}
            gene_level_id, _ = GeneLevelId.objects.get_or_create(pk=hgnc.pk, defaults=defaults)
            return gene_level_id

        if gene_level_id := GeneLevelId.objects.filter(symbol_str=symbol_str).first():
            return gene_level_id

        # Allocating our own pk, so another worker can take the number between the max() and the insert
        for _ in range(GeneLevelId.CUSTOM_ID_RETRIES):
            try:
                with transaction.atomic():
                    last = GeneLevelId.objects.filter(pk__gte=GeneLevelId.CUSTOM_ID_START).order_by("-pk").first()
                    pk = last.pk + 1 if last else GeneLevelId.CUSTOM_ID_START
                    return GeneLevelId.objects.create(pk=pk, symbol_str=symbol_str,
                                                      gene_symbol_id=gene_symbol_id)
            except IntegrityError:
                pass
            if gene_level_id := GeneLevelId.objects.filter(symbol_str=symbol_str).first():
                return gene_level_id
        raise IntegrityError(f"Could not allocate a GeneLevelId id for '{symbol_str}'")


class GeneCopyNumberEventKind(models.TextChoices):
    """ Which way a whole-gene copy number call went. GAIN rather than amplification because the
        threshold that makes a gain an amplification is the lab's, not ours - identity never encodes
        one, and the copy ratio stays per-sample on CohortGenotype. """

    GAIN = "G", "Gain"
    LOSS = "L", "Loss"

    @property
    def alt_kind(self) -> str:
        """ The GeneLevelSymbolicAlt this kind's Variant alt is written with """
        return GENE_COPY_NUMBER_ALT_KINDS[self]

    @property
    def canonical_word(self) -> str:
        """ 'amplification' is what the reports say; 'gain' is what we store """
        return GENE_COPY_NUMBER_CANONICAL_WORDS[self]

    @staticmethod
    def from_alt_kind(alt_kind: str) -> Optional['GeneCopyNumberEventKind']:
        return GENE_COPY_NUMBER_KINDS_BY_ALT.get(alt_kind)


GENE_COPY_NUMBER_ALT_KINDS = {
    GeneCopyNumberEventKind.GAIN: GeneLevelSymbolicAlt.GAIN,
    GeneCopyNumberEventKind.LOSS: GeneLevelSymbolicAlt.LOSS,
}
GENE_COPY_NUMBER_KINDS_BY_ALT = {alt_kind: kind for kind, alt_kind in GENE_COPY_NUMBER_ALT_KINDS.items()}
GENE_COPY_NUMBER_CANONICAL_WORDS = {
    GeneCopyNumberEventKind.GAIN: "amplification",
    GeneCopyNumberEventKind.LOSS: "loss",
}


def gene_copy_number_canonical_str(gene: GeneLevelId, kind: str) -> str:
    """ 'EGFR amplification' / 'EGFR loss' - the words on the reports and in the caller's combined
        output, and the form to send anywhere off this deployment (@see GeneLevelId). 'deletion' is
        accepted on input but never written, since as output it reads as a coordinate event. """
    return f"{gene.symbol_str} {GeneCopyNumberEventKind(kind).canonical_word}"


class GeneCopyNumberEvent(models.Model):
    """ A whole-gene copy number call, one-to-one with the Variant carrying it - the twin of
        GeneFusion, and @see snpdb.gene_level_variants for why it is stored as a Variant.

        Identity is the gene plus the direction, not the coordinates: the segment a caller reports is
        the panel's target window rather than the event, so it moves with the manifest and differs
        between assays for what everyone calls "EGFR amplification". The copy ratio is per
        observation and lives on each sample's CohortGenotype, where the copy number column reads it
        (@see VCF.copy_number_field). """

    variant = models.OneToOneField('snpdb.Variant', on_delete=CASCADE)
    gene = models.ForeignKey(GeneLevelId, related_name='copy_number_events', on_delete=PROTECT)
    kind = models.CharField(max_length=1, choices=GeneCopyNumberEventKind.choices)

    def __str__(self):
        return self.canonical_str

    def get_absolute_url(self):
        return self.variant.get_absolute_url()

    @property
    def canonical_str(self) -> str:
        """ The form to display and to send anywhere off this deployment - @see GeneLevelId """
        return gene_copy_number_canonical_str(self.gene, self.kind)

    @property
    def gene_level_ids(self) -> list[GeneLevelId]:
        """ The genes this event is about, so gene lists and annotation read it the way they read a
            fusion's partners """
        return [self.gene]

    def clean(self):
        super().clean()
        parsed = GeneLevelSymbolicAlt.parse(self.variant.alt.seq)
        if parsed is None:
            raise ValidationError(f"Variant alt '{self.variant.alt.seq}' is not a gene-level alt")

        alt_kind, _namespace, gene_id = parsed
        if alt_kind != GeneCopyNumberEventKind(self.kind).alt_kind:
            raise ValidationError(f"Variant alt '{self.variant.alt.seq}' does not match {self.kind=}")
        if gene_id != self.gene_id:
            raise ValidationError(f"Variant alt '{self.variant.alt.seq}' does not match {self.gene_id=}")
        if self.variant.locus.position != self.gene_id:
            raise ValidationError(f"Variant position {self.variant.locus.position} is not {self.gene_id=}")
