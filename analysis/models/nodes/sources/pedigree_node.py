from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models import Q
from django.db.models.deletion import SET_NULL

from analysis.models.nodes.sources import AbstractCohortBasedNode
from analysis.models.nodes.sources._stats_cache import (
    get_cached_label_count_for_cohort, get_handler_for_node, UNCACHEABLE,
)
from patients.models_enums import Zygosity
from pedigree.models import Pedigree, PedigreeInheritance, CohortSamplePedFileRecord


class PedigreeNode(AbstractCohortBasedNode):
    pedigree = models.ForeignKey(Pedigree, null=True, on_delete=SET_NULL)
    inheritance_model = models.CharField(max_length=2, choices=PedigreeInheritance.choices)
    require_zygosity = models.BooleanField(default=True)
    min_inputs = 0
    max_inputs = 0

    def _get_cohort(self):
        cohort = None
        if self.pedigree:
            cohort = self.pedigree.cohort
        return cohort

    def _has_filters_that_affect_label_counts(self) -> bool:
        # Inheritance model is a per-pedigree shape we don't precompute today,
        # so any inheritance setting defeats the cache. Quality filters
        # (min_ad/dp/gq, max_pl) still defeat as well.
        if super()._has_filters_that_affect_label_counts():
            return True
        return bool(self.inheritance_model)

    def _get_cached_label_count(self, label):
        if self.pedigree is None:
            return None
        if self._has_filters_that_affect_label_counts():
            return None
        filter_code = self.get_filter_code()
        if filter_code not in (0, 1):
            return None
        handler = get_handler_for_node(self)
        filter_key = handler.filter_key_for_node(self)
        if filter_key is UNCACHEABLE:
            return None
        return get_cached_label_count_for_cohort(
            cohort=self.pedigree.cohort,
            sample=None,
            filter_key=filter_key,
            annotation_version=self.analysis.annotation_version,
            passing_filter=bool(filter_code),
            zygosities=self._cached_label_count_zygosities(),
            label=label,
        )

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        cohort, arg_q_dict = self.get_cohort_and_arg_q_dict()
        if cohort:
            q = None
            if self.inheritance_model == PedigreeInheritance.AUTOSOMAL_RECESSIVE:
                q = self.get_recessive_q(cohort.cohort_genotype_collection)
            elif self.inheritance_model == PedigreeInheritance.AUTOSOMAL_DOMINANT:
                q = self.get_dominant_q(cohort.cohort_genotype_collection)

            if q:
                cohort_arg_q_dict = {
                    self.cohort_genotype_collection.cohortgenotype_alias: {str(q): q}
                }
                self.merge_arg_q_dicts(arg_q_dict, cohort_arg_q_dict)
        return arg_q_dict

    def get_affected_unaffected_sample_zygosities_dict(self, unaffected_zygosities, affected_zygosities):
        sample_zygosities_dict = {}
        sample_require_zygosity_dict = {}

        for cspfr in CohortSamplePedFileRecord.objects.filter(pedigree=self.pedigree):
            if cspfr.ped_file_record.affection:
                zyg = affected_zygosities
            else:
                zyg = unaffected_zygosities
            sample = cspfr.cohort_sample.sample
            sample_zygosities_dict[sample] = zyg
            sample_require_zygosity_dict[sample] = self.require_zygosity
        return sample_zygosities_dict, sample_require_zygosity_dict

    def get_recessive_q(self, cohort_genotype_collection):
        unaffected = {Zygosity.HET}
        affected = {Zygosity.HOM_ALT}
        sample_zygosities_dict, sample_require_zygosity_dict = self.get_affected_unaffected_sample_zygosities_dict(unaffected, affected)
        q = cohort_genotype_collection.get_zygosity_q(sample_zygosities_dict, sample_require_zygosity_dict)
        return q

    def get_dominant_q(self, cohort_genotype_collection):
        unaffected = {r'\.'}
        affected = {Zygosity.HET, Zygosity.HOM_ALT}
        sample_zygosities_dict, sample_require_zygosity_dict = self.get_affected_unaffected_sample_zygosities_dict(unaffected, affected)
        return cohort_genotype_collection.get_zygosity_q(sample_zygosities_dict, sample_require_zygosity_dict)

    def _get_method_summary(self):
        method_summary = ''
        if self.pedigree:
            pass
        return method_summary

    def get_node_name(self):
        name = ''
        if self.pedigree:
            return f"{self.pedigree.name} ({self.inheritance_model})"
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Variants from a family of samples, filtered by genotype according to inheritance models"

    @staticmethod
    def get_node_class_label():
        return "Pedigree"

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if not self.pedigree:
            errors.append("No pedigree selected.")
        else:
            errors.extend(self._get_genome_build_errors("pedigree", self.pedigree.genome_build))

        return errors

    def __str__(self):
        return self.get_node_name()


auditlog.register(Pedigree)
