""" Builds the PhenotypeDescription a Patient or Cohort asks for when its phenotype text changes
    (HasPhenotypeDescriptionMixin.process_phenotype_if_changed) """
from django.dispatch import receiver

from patients.models.has_phenotype_description_mixin import phenotype_description_wanted_signal
from patients.phenotype_matching import create_phenotype_description


@receiver(phenotype_description_wanted_signal)
def phenotype_description_wanted_handler(sender, owner, text, approved_by, phenotype_matcher, defer_processing,
                                         **kwargs):  # pylint: disable=unused-argument
    create_phenotype_description(text, owner=owner, approved_by=approved_by, phenotype_matcher=phenotype_matcher,
                                 defer_processing=defer_processing)
