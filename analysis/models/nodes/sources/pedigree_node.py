from typing import Optional

from django.db import models
from django.db.models import Q
from django.db.models.deletion import SET_NULL

from analysis.models.nodes.sources import AbstractCohortBasedNode
from patients.models_enums import Zygosity
from pedigree.models import Pedigree, PedigreeInheritance, CohortSamplePedFileRecord


class PedigreeNode(AbstractCohortBasedNode):
    pedigree = models.ForeignKey(Pedigree, null=True, on_delete=SET_NULL)
    inheritance_model = models.CharField(max_length=2, choices=PedigreeInheritance.CHOICES)
    require_zygosity = models.BooleanField(default=True)
    min_inputs = 0
    max_inputs = 0

    def _get_cohort(self):
        cohort = None
        if self.pedigree:
            cohort = self.pedigree.cohort
        return cohort

    def _get_node_q(self) -> Optional[Q]:
        cohort, q = self.get_cohort_and_q()
        if cohort:
            if self.inheritance_model == PedigreeInheritance.AUTOSOMAL_RECESSIVE:
                q &= self.get_recessive_q(cohort.cohort_genotype_collection)
            elif self.inheritance_model == PedigreeInheritance.AUTOSOMAL_DOMINANT:
                q &= self.get_dominant_q(cohort.cohort_genotype_collection)
        return q

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
    def get_node_class_label():
        return "Pedigree"

    def _get_configuration_errors(self):
        errors = []
        if not self.pedigree:
            errors.append("No pedigree selected.")

        return errors

    def __str__(self):
        return self.get_node_name()
