from typing import List

from django.db.models import QuerySet
from lazy import lazy


class HasPhenotypeDescriptionMixin:

    def _get_phenotype_input_text_field(self):
        """ Each subclass needs to implement - returns a string """
        raise NotImplementedError()

    def _get_phenotype_description_relation_class_and_kwargs(self):
        """ Each subclass needs to implement the way to get their PhenotypeDescription object """
        raise NotImplementedError()

    @property
    def phenotype_input_text(self):
        return getattr(self, self._get_phenotype_input_text_field())

    @phenotype_input_text.setter
    def phenotype_input_text(self, value):
        setattr(self, self._get_phenotype_input_text_field(), value)

    @lazy
    def phenotype_description_relation(self):
        try:
            klass, kwargs = self._get_phenotype_description_relation_class_and_kwargs()
            _phenotype_description_relation = klass.objects.get(**kwargs)
        except:
            _phenotype_description_relation = None
        return _phenotype_description_relation

    @property
    def phenotype_description(self):
        if self.phenotype_description_relation:
            return self.phenotype_description_relation.phenotype_description
        return None

    def get_ontology_term_ids(self) -> List[str]:
        if self.phenotype_description:
            terms = self.phenotype_description.get_ontology_term_ids()
        else:
            terms = []
        return terms

    def get_gene_symbols(self) -> QuerySet:
        from genes.models import Gene

        if self.phenotype_description:
            gene_qs = self.phenotype_description.get_gene_symbols()
        else:
            gene_qs = Gene.objects.none()
        return gene_qs

    def process_phenotype_if_changed(self, phenotype_matcher=None, phenotype_approval_user=None):
        """ pass in phenotype_matcher to save re-loading
            if you don't pass in phenotype_approval_user assumed it is done automatically and thus needs user approval
            returns whether phenotype changed """

        # Stop circular import
        from annotation.phenotype_matcher import PhenotypeMatcher
        from annotation.phenotype_matching import create_phenotype_description

        if phenotype_matcher is None:
            phenotype_matcher = PhenotypeMatcher()

        phenotype_input_text = self.phenotype_input_text
        phenotype_description_relation = self.phenotype_description_relation
        phenotype_description = None

        if phenotype_description_relation:
            phenotype_description = phenotype_description_relation.phenotype_description
            changed = phenotype_description and phenotype_description.original_text != phenotype_input_text
            if changed:
                phenotype_description.delete()
                phenotype_description = None
        else:
            klass, kwargs = self._get_phenotype_description_relation_class_and_kwargs()
            phenotype_description_relation = klass(**kwargs)  # Create a new one

        parsed_phenotypes = False
        if phenotype_input_text:
            if phenotype_description:  # Not deleted - so no change needed
                pass
            else:
                # TODO: Do as async job??
                phenotype_description_relation.phenotype_description = create_phenotype_description(phenotype_input_text, phenotype_matcher)
                phenotype_description_relation.approved_by = phenotype_approval_user
                phenotype_description_relation.save()

                parsed_phenotypes = True

        return parsed_phenotypes

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
        if self.phenotype_input_text:
            self.phenotype_input_text = self.phenotype_input_text.replace('\r', '')

        kwargs = HasPhenotypeDescriptionMixin.pop_kwargs(kwargs_dict)
        if kwargs.pop("check_patient_text_phenotype", False):
            self.process_phenotype_if_changed(**kwargs)
