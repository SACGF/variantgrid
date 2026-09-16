from typing import Optional

from django.core.exceptions import ValidationError
from django.db import models
from django.db.models.deletion import CASCADE, PROTECT

from genes.models.models_gene_level import GeneLevelId
from library.genomics.vcf_enums import GeneLevelSymbolicAlt


def fusion_canonical_str(anchor: GeneLevelId, partner: Optional[GeneLevelId]) -> str:
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
    anchor = models.ForeignKey(GeneLevelId, related_name='fusions_as_anchor', on_delete=PROTECT)
    partner = models.ForeignKey(GeneLevelId, null=True, related_name='fusions_as_partner', on_delete=PROTECT)
    is_ordered = models.BooleanField(default=False)

    def __str__(self):
        return self.canonical_str

    def get_absolute_url(self):
        return self.variant.get_absolute_url()

    @property
    def canonical_str(self) -> str:
        """ The form to display and to send anywhere off this deployment - @see GeneLevelId """
        return fusion_canonical_str(self.anchor, self.partner)

    @property
    def gene_level_ids(self) -> list[GeneLevelId]:
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
