from typing import Optional

from auditlog.registry import auditlog
from django.db import models
from django.db.models import Q
from django.db.models.deletion import SET_NULL

from analysis.models.nodes.cohort_mixin import merged_cohort_warnings
from analysis.models.nodes.family_inheritance import AbstractFamilyInheritance
from analysis.models.nodes.node_display import NodeIcon
from analysis.models.nodes.sources.cohort_node import AbstractCohortBasedNode
from patients.enums import Zygosity, ZygosityRequirement
from pedigree.enums import PedigreeInheritance
from pedigree.models import CohortSamplePedFileRecord, Pedigree


class PedigreeNode(AbstractCohortBasedNode):
    pedigree = models.ForeignKey(Pedigree, null=True, on_delete=SET_NULL)
    inheritance_model = models.CharField(max_length=2, choices=PedigreeInheritance.choices)
    zygosity_requirement = models.CharField(max_length=1, choices=ZygosityRequirement.choices,
                                            default=ZygosityRequirement.CALL_REQUIRED)
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
        """ A record with unknown affection (None) is left unconstrained - the same rule as
            Pedigree.get_samples(affected=False), which excludes it """
        sample_zygosities_dict = {}
        sample_zygosity_requirements = {}

        for cspfr in CohortSamplePedFileRecord.objects.filter(pedigree=self.pedigree):
            affection = cspfr.ped_file_record.affection
            if affection is None:
                continue
            if affection:
                zyg = affected_zygosities
            else:
                zyg = unaffected_zygosities
            sample = cspfr.cohort_sample.sample
            sample_zygosities_dict[sample] = zyg
            sample_zygosity_requirements[sample] = ZygosityRequirement(self.zygosity_requirement)
        return sample_zygosities_dict, sample_zygosity_requirements

    def get_recessive_q(self, cohort_genotype_collection):
        unaffected = AbstractFamilyInheritance.NOT_HOM_ALT
        affected = {Zygosity.HOM_ALT}
        sample_zygosities_dict, sample_zygosity_requirements = self.get_affected_unaffected_sample_zygosities_dict(unaffected, affected)
        return cohort_genotype_collection.get_zygosity_q(sample_zygosities_dict, sample_zygosity_requirements)

    def get_dominant_q(self, cohort_genotype_collection):
        unaffected = AbstractFamilyInheritance.NO_VARIANT
        affected = AbstractFamilyInheritance.HAS_VARIANT
        sample_zygosities_dict, sample_zygosity_requirements = self.get_affected_unaffected_sample_zygosities_dict(unaffected, affected)
        return cohort_genotype_collection.get_zygosity_q(sample_zygosities_dict, sample_zygosity_requirements)

    def _get_method_summary(self):
        method_summary = ''
        if self.pedigree:
            pass
        return method_summary

    def get_warnings(self) -> list[str]:
        warnings = super().get_warnings()
        if self.pedigree:
            warnings.extend(merged_cohort_warnings(self.pedigree.cohort, [ZygosityRequirement(self.zygosity_requirement)]))
        return warnings

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

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        return NodeIcon(symbol="node-icon-pedigree")

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
