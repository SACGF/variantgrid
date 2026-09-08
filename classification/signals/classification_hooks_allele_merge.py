from django.dispatch import receiver

from classification.models import Classification
from classification.models.clinical_context_utils import (
    classifications_needing_rehoming,
    rehome_classifications,
)
from snpdb.models import Allele, allele_merged_signal


@receiver(allele_merged_signal, sender=Allele)
def allele_merged_handler(sender, old_allele: Allele, new_allele: Allele, **kwargs):  # pylint: disable=unused-argument
    # ALLELES MERGED
    # merge() moves Classification.allele across in bulk, so the clinical context and grouping they were
    # filed under are still those of the allele that was merged away
    moved_qs = Classification.objects.filter(allele=new_allele)
    rehome_classifications(classifications_needing_rehoming(moved_qs),
                           force_recalc_text=f"{old_allele} merged into {new_allele}")
