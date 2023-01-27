from django.contrib.auth.models import User
from django.dispatch import receiver

from classification.models import classification_post_publish_signal, Classification, ClassificationModification, \
    variants_classification_changed_signal, allele_info_changed_signal, ImportedAlleleInfo
from library.utils import DebugTimer
from snpdb.models import Allele

"""
Converts various events to variants_classification_changed_signal

"""


def _send_variants_classification_changed_for(allele: Allele):
    for variant_allele in allele.variant_alleles():
        variants_classification_changed_signal.send(sender=Classification,
                                                    variants=[variant_allele.variant],
                                                    genome_build=variant_allele.genome_build)


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              user: User,
              debug_timer: DebugTimer,
              **kwargs):  # pylint: disable=unused-argument
    if not previously_published or (previously_published.clinical_significance != newly_published.clinical_significance):
        if allele_info := classification.allele_info:
            if allele := allele_info.allele:
                _send_variants_classification_changed_for(allele)


@receiver(allele_info_changed_signal, sender=ImportedAlleleInfo)
def allele_info_changed_update_classifications(sender, allele_info: ImportedAlleleInfo, **kwargs):
    """
    An AlleleInfo has had its matching to a variant changed, apply the change to classifications linked
    to the allele info, as well as calling variants_classification_changed_signal
    """
    for classification in allele_info.classification_set.all():
        classification.apply_allele_info_to_classification()
        classification.save()

    # update variants_classification_changed_signal
    if allele := allele_info.allele:
        _send_variants_classification_changed_for(allele)
