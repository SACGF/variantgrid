from typing import Optional

from django.core.exceptions import ValidationError
from django.db import IntegrityError, models, transaction
from django.db.models import TextField
from django.db.models.deletion import CASCADE, PROTECT, SET_NULL

from genes.models.models_gene import HGNC, GeneSymbol
from library.genomics.vcf_enums import GeneIdNamespace, GeneLevelSymbolicAlt


class FusionGeneId(models.Model):
    """ A stable number for one side of a gene fusion - the pk is what a fusion Variant carries as its
        Locus.position and inside its symbolic alt (@see GeneLevelSymbolicAlt), so this table is the
        identifier space fusion identity is built on. @see snpdb.gene_level_variants for why a fusion
        is stored as a Variant in the first place.

        pk is the HGNC ID where the gene has one, so the identity of an ordinary fusion is the same
        number on every deployment. Symbols HGNC doesn't carry - clone-based identifiers like
        RP11-458D21.5, which turn up as fusion partners routinely - get a pk allocated above
        CUSTOM_ID_START, so every call the caller made can still become a Variant.

        Custom numbers are local to this deployment, which is why the alt namespaces them as GENE
        rather than HGNC. Anything leaving the system - display, export, a classification sent to
        another instance - uses symbol_str rather than the pk, and the receiver resolves it back
        through its own table. @see GeneLevelSymbolicAlt for where the numbers are used, and
        GeneFusion.canonical_str for the string form.

        A symbol that later gains an HGNC entry has its hgnc FK filled in, which improves annotation.
        Existing variants keep the custom pk they were created with - identity, once handed out, is
        fixed, the same way a Variant's coordinates are. """

    CUSTOM_ID_START = 1_000_000
    CUSTOM_ID_RETRIES = 5

    # The name to use whenever this leaves the system - the approved symbol where we have one,
    # otherwise the name exactly as the caller wrote it
    symbol_str = TextField(unique=True, db_collation='case_insensitive')
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=SET_NULL)
    hgnc = models.ForeignKey(HGNC, null=True, on_delete=SET_NULL)

    def __str__(self):
        return self.symbol_str

    @property
    def is_custom(self) -> bool:
        return self.pk >= FusionGeneId.CUSTOM_ID_START

    @property
    def alt_namespace(self) -> str:
        return GeneIdNamespace.GENE if self.is_custom else GeneIdNamespace.HGNC

    @staticmethod
    def get_or_create_for_symbol(symbol_str: str, gene_symbol_id: Optional[str],
                                 hgnc: Optional[HGNC]) -> 'FusionGeneId':
        """ symbol_str is what to show and send - callers resolve aliases first, so SEPT14 arrives
            here as SEPTIN14 and collapses onto the one row """

        if hgnc is not None:
            defaults = {"symbol_str": symbol_str, "gene_symbol_id": gene_symbol_id, "hgnc": hgnc}
            fusion_gene_id, _ = FusionGeneId.objects.get_or_create(pk=hgnc.pk, defaults=defaults)
            return fusion_gene_id

        if fusion_gene_id := FusionGeneId.objects.filter(symbol_str=symbol_str).first():
            return fusion_gene_id

        # Allocating our own pk, so another worker can take the number between the max() and the insert
        for _ in range(FusionGeneId.CUSTOM_ID_RETRIES):
            try:
                with transaction.atomic():
                    last = FusionGeneId.objects.filter(pk__gte=FusionGeneId.CUSTOM_ID_START).order_by("-pk").first()
                    pk = last.pk + 1 if last else FusionGeneId.CUSTOM_ID_START
                    return FusionGeneId.objects.create(pk=pk, symbol_str=symbol_str,
                                                        gene_symbol_id=gene_symbol_id)
            except IntegrityError:
                pass
            if fusion_gene_id := FusionGeneId.objects.filter(symbol_str=symbol_str).first():
                return fusion_gene_id
        raise IntegrityError(f"Could not allocate a FusionGeneId id for '{symbol_str}'")

    CUSTOM_ID_RETRIES = 5


def fusion_canonical_str(anchor: FusionGeneId, partner: Optional[FusionGeneId]) -> str:
    """ 'BCR::ABL1' - the VICC gene-level form, which HGNC and HGVS both point at for fusions.
        '::' is the fusion separator; a single hyphen means a read-through transcript, which is a
        different event, so it is never written here. An unnamed partner shows as '?'. """
    partner_str = partner.symbol_str if partner else "?"
    return f"{anchor.symbol_str}::{partner_str}"


class GeneFusion(models.Model):
    """ A gene fusion, one-to-one with the Variant carrying it. The Variant is what makes a fusion
        reachable from gene lists, comp-het, grids and classifications - @see snpdb.gene_level_variants,
        which explains why that is how these are stored.

        Identity is the gene pair, not the breakpoints: one caller reports ENTPD3-RPL14 three times
        with three different 5' breakpoints, so coordinates are per-observation and live in
        CohortGenotype.info. Reciprocal fusions are separate rows - the 5' promoter drives a
        different protein, so BCR-ABL1 and ABL1-BCR are not the same event.

        The anchor is the 5' partner when is_ordered, otherwise the lower-id gene of the pair.
        partner is null only where the caller named one gene and left the other unspecified. """

    variant = models.OneToOneField('snpdb.Variant', on_delete=CASCADE)
    anchor = models.ForeignKey(FusionGeneId, related_name='fusions_as_anchor', on_delete=PROTECT)
    partner = models.ForeignKey(FusionGeneId, null=True, related_name='fusions_as_partner', on_delete=PROTECT)
    is_ordered = models.BooleanField(default=False)

    @property
    def canonical_str(self) -> str:
        """ The form to display and to send anywhere off this deployment - @see FusionGeneId """
        return fusion_canonical_str(self.anchor, self.partner)

    def __str__(self):
        return self.canonical_str

    def get_absolute_url(self):
        return self.variant.get_absolute_url()

    @property
    def fusion_gene_ids(self) -> list[FusionGeneId]:
        """ Both partners, so gene lists and comp-het find the fusion from either side """
        genes = [self.anchor]
        if self.partner:
            genes.append(self.partner)
        return genes

    def clean(self):
        super().clean()
        parsed = GeneLevelSymbolicAlt.parse(self.variant.alt.seq)
        if parsed is None:
            raise ValidationError(f"Variant alt '{self.variant.alt.seq}' is not a gene-level alt")

        kind, _namespace, partner_id = parsed
        expected_kind = GeneLevelSymbolicAlt.FUSION if self.is_ordered else GeneLevelSymbolicAlt.FUSION_UNORDERED
        if kind != expected_kind:
            raise ValidationError(f"Variant alt '{self.variant.alt.seq}' does not match {self.is_ordered=}")
        if partner_id != self.partner_id:
            raise ValidationError(f"Variant alt '{self.variant.alt.seq}' does not match {self.partner_id=}")
        if self.variant.locus.position != self.anchor_id:
            raise ValidationError(f"Variant position {self.variant.locus.position} is not {self.anchor_id=}")
