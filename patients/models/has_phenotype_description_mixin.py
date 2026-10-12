"""
HasPhenotypeDescriptionMixin: a model with a `phenotype` TextField whose text is split into sentences and
matched to ontology terms. Its PhenotypeDescription points back at it through the one-to-one named after the
model (PhenotypeDescription.patient, PhenotypeDescription.cohort), so deleting the owner deletes it.

Entry points: save_phenotype (from the owner's save), process_phenotype_if_changed (bulk matching).
New text goes out on phenotype_description_wanted_signal, and patients/signals/phenotype_description.py builds
the description: the matcher needs ontology, which sits above snpdb, and snpdb imports this module.
"""
from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from django.db.models import QuerySet
from django.dispatch import Signal

# The reverse accessor PhenotypeDescription's owner fields give the owner
PHENOTYPE_DESCRIPTION = "phenotype_description"

# Sent (owner, text, approved_by, phenotype_matcher, defer_processing) when owner's text needs a new description
phenotype_description_wanted_signal = Signal()


def phenotype_text_is_excluded(text: str) -> bool:
    """ True if text contains settings.PATIENT_PHENOTYPE_EXCLUDE_STRING — caller should skip
        persisting auto-matched phenotype terms (human review required first). """
    exclude_string = getattr(settings, "PATIENT_PHENOTYPE_EXCLUDE_STRING", None)
    return bool(exclude_string and text and exclude_string in text)


class HasPhenotypeDescriptionMixin:
    """ For a model with a `phenotype` TextField that PhenotypeDescription has an owner field for, named after the
        model (self._meta.model_name) """

    def get_phenotype_description(self):
        """ The PhenotypeDescription of the current text, or None - cached on the instance either way """
        try:
            return getattr(self, PHENOTYPE_DESCRIPTION)
        except ObjectDoesNotExist:
            return None

    def get_ontology_term_ids(self) -> list[str]:
        if phenotype_description := self.get_phenotype_description():
            return phenotype_description.get_ontology_term_ids()
        return []

    def get_gene_symbols(self, ontology_version) -> QuerySet:
        return ontology_version.cached_gene_symbols_for_terms_tuple(tuple(self.get_ontology_term_ids()))

    def process_phenotype_if_changed(self, phenotype_matcher=None, phenotype_approval_user=None,
                                     defer_processing=False) -> bool:
        """ pass in phenotype_matcher to save re-loading
            if you don't pass in phenotype_approval_user assumed it is done automatically and thus needs user approval
            if defer_processing is True the PhenotypeDescription/TextPhenotype rows are created but NLP matching is
            skipped (used by bulk_patient_phenotype_matching to batch the heavy work for parallel execution)
            returns whether a new description was made """

        phenotype_description = self.get_phenotype_description()
        if phenotype_description and phenotype_description.original_text != self.phenotype:
            phenotype_description.delete()  # Its sentences and approval go with it
            self._state.fields_cache.pop(PHENOTYPE_DESCRIPTION, None)
            phenotype_description = None

        if phenotype_description or not self.phenotype:
            return False
        if phenotype_text_is_excluded(self.phenotype):
            return False  # Marker present — defer persistence until a human cleans up the text.

        phenotype_description_wanted_signal.send(sender=type(self), owner=self, text=self.phenotype,
                                                 approved_by=phenotype_approval_user,
                                                 phenotype_matcher=phenotype_matcher,
                                                 defer_processing=defer_processing)
        return True

    @staticmethod
    def pop_kwargs(kwargs_dict):
        """ remove kwargs (so save for other model doesn't fail """
        defaults = {"check_patient_text_phenotype": True,
                    "phenotype_approval_user": None,
                    "phenotype_matcher": None}
        return {k: kwargs_dict.pop(k, v) for k, v in defaults.items()}

    def save_phenotype(self, kwargs_dict):
        """ Pass kwargs_dict as dict - will pop fields it uses:
            "check_patient_text_phenotype" and "phenotype_approval_user" """

        # Some browsers send Text inputs with \r\n - while AJAX sends it as \n
        # strip \r to keep it consistent so that highlighting offsets line up
        if self.phenotype:
            self.phenotype = self.phenotype.replace('\r', '')

        kwargs = HasPhenotypeDescriptionMixin.pop_kwargs(kwargs_dict)
        if kwargs.pop("check_patient_text_phenotype", False):
            self.process_phenotype_if_changed(**kwargs)
