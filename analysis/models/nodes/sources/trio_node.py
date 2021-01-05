from typing import Optional, List

from django.db import models
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.enums import TrioInheritance, NodeErrorSource, AnalysisTemplateType
from analysis.models.nodes.sources import AbstractCohortBasedNode
from annotation.models.models import VariantTranscriptAnnotation
from patients.models_enums import Zygosity, Sex
from snpdb.models import Trio, Sample


class TrioNode(AbstractCohortBasedNode):
    trio = models.ForeignKey(Trio, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=TrioInheritance.choices, default=TrioInheritance.RECESSIVE)
    require_zygosity = models.BooleanField(default=True)

    NO_VARIANT = {Zygosity.MISSING, Zygosity.HOM_REF}  # 2 different het would be "missing" (as has no ref)
    HAS_VARIANT = {Zygosity.HET, Zygosity.HOM_ALT}
    UNAFFECTED_AND_AFFECTED_ZYGOSITIES = [NO_VARIANT, HAS_VARIANT]

    @property
    def min_inputs(self):
        return self.max_inputs

    @property
    def max_inputs(self):
        if self.inheritance == TrioInheritance.COMPOUND_HET:
            return 1
        return 0

    @staticmethod
    def get_trio_inheritance_errors(trio: Trio, inheritance) -> List[str]:
        errors = []
        if trio:
            if inheritance == TrioInheritance.DOMINANT:
                if not trio.mother_affected or trio.father_affected:
                    errors.append("Dominant inheritance requires an affected parent")
            elif inheritance == TrioInheritance.XLINKED_RECESSIVE:
                proband_sample = trio.proband.sample
                proband_is_female = proband_sample.patient and proband_sample.patient.sex == Sex.FEMALE
                if proband_is_female:
                    errors.append("X-linked recessive inheritance doesn't currently work with female proband")
                elif trio.mother_affected:
                    errors.append("X-linked recessive inheritance requires an unaffected mother")
        return errors

    def get_errors(self, include_parent_errors=True, flat=False):
        errors = super().get_errors(include_parent_errors=include_parent_errors)
        # Allow template to configure anything w/o checks
        if self.analysis.template_type != AnalysisTemplateType.TEMPLATE:
            if trio_errors := self.get_trio_inheritance_errors(self.trio, self.inheritance):
                errors.extend(((NodeErrorSource.CONFIGURATION, e) for e in trio_errors))
        if flat:
            errors = self.flatten_errors(errors)
        return errors

    def _get_cohort(self):
        cohort = None
        if self.trio:
            cohort = self.trio.cohort
        return cohort

    def modifies_parents(self):
        return self.trio is not None

    def _get_node_q(self) -> Optional[Q]:
        cohort, q = self.get_cohort_and_q()
        if cohort:
            cohort_genotype_collection = cohort.cohort_genotype_collection
            if self.inheritance == TrioInheritance.COMPOUND_HET:
                parent = self.get_single_parent()
                q &= self.get_compound_het_q(cohort_genotype_collection, parent)
            else:
                if self.inheritance == TrioInheritance.RECESSIVE:
                    q &= self.get_recessive_q(cohort_genotype_collection)
                elif self.inheritance == TrioInheritance.DOMINANT:
                    q &= self.get_dominant_q(cohort_genotype_collection)
                elif self.inheritance == TrioInheritance.DENOVO:
                    q &= self.get_denovo_q(cohort_genotype_collection)
                elif self.inheritance == TrioInheritance.PROBAND_HET:
                    q &= self.get_proband_het_q(cohort_genotype_collection)
                elif self.inheritance == TrioInheritance.XLINKED_RECESSIVE:
                    q &= self.get_xlinked_recessive_q(cohort_genotype_collection)

            q &= self.get_vcf_locus_filters_q()
        return q

    def get_zyg_q(self, cohort_genotype_collection, trio_zyg_data):
        """ trio_zyg_data = tuple of mum_zyg_set, dad_zyg_set, proband_zyg_set """
        mother_sample = self.trio.mother.sample
        father_sample = self.trio.father.sample
        proband_sample = self.trio.proband.sample

        mum_dad_proband = [mother_sample, father_sample, proband_sample]
        sample_zygosities_dict = dict(zip(mum_dad_proband, trio_zyg_data))
        sample_require_zygosity_dict = {mother_sample: self.require_zygosity,
                                        father_sample: self.require_zygosity,
                                        proband_sample: True}  # 947 - Always require zygosity for Proband

        return cohort_genotype_collection.get_zygosity_q(sample_zygosities_dict, sample_require_zygosity_dict)

    def get_recessive_q(self, cohort_genotype_collection):
        zyg_data = ({Zygosity.HET}, {Zygosity.HET}, {Zygosity.HOM_ALT})
        return self.get_zyg_q(cohort_genotype_collection, zyg_data)

    def get_compound_het_q(self, cohort_genotype_collection, parent):

        mum_but_not_dad = self.get_zyg_q(cohort_genotype_collection, ({Zygosity.HET}, self.NO_VARIANT, {Zygosity.HET}))
        dad_but_not_mum = self.get_zyg_q(cohort_genotype_collection, (self.NO_VARIANT, {Zygosity.HET}, {Zygosity.HET}))
        comp_het_q = mum_but_not_dad | dad_but_not_mum
        # Need to pass in kwargs in case we have parent (eg Venn node) that doesn't have same cohort annotation kwargs
        annotation_kwargs = self.get_annotation_kwargs()

        def get_parent_genes(q):
            qs = parent.get_queryset(q, extra_annotation_kwargs=annotation_kwargs)
            return qs.values_list("varianttranscriptannotation__gene", flat=True).distinct()

        common_genes = set(get_parent_genes(mum_but_not_dad)) & set(get_parent_genes(dad_but_not_mum))
        comp_het_genes = VariantTranscriptAnnotation.get_overlapping_genes_q(common_genes)
        return parent.get_q() & comp_het_q & comp_het_genes

    def get_dominant_q(self, cohort_genotype_collection):
        mother_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(self.trio.mother_affected)]
        father_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(self.trio.father_affected)]
        zyg_data = (mother_zyg, father_zyg, self.HAS_VARIANT)

        return self.get_zyg_q(cohort_genotype_collection, zyg_data)

    def get_denovo_q(self, cohort_genotype_collection):
        zyg_data = (self.NO_VARIANT, self.NO_VARIANT, self.HAS_VARIANT)
        return self.get_zyg_q(cohort_genotype_collection, zyg_data)

    def get_proband_het_q(self, cohort_genotype_collection):
        """ Equivalent of SampleNode, but in TrioNode so analysis templates only need to have Trio as argument """
        zyg_data = ({}, {}, {Zygosity.HET})
        return self.get_zyg_q(cohort_genotype_collection, zyg_data)

    def get_xlinked_recessive_q(self, cohort_genotype_collection):
        q_x = Q(locus__contig__name='X')  # will work for hg19 and GRCh38
        zyg_data = ({Zygosity.HET}, {}, {Zygosity.HOM_ALT})
        return q_x & self.get_zyg_q(cohort_genotype_collection, zyg_data)

    def _get_method_summary(self):
        return "TODO: method summary"

    def get_node_name(self):
        name_parts = [TrioInheritance(self.inheritance).label]
        if not self.require_zygosity:
            name_parts.append(' (non strict)')

        filter_description = self.get_filter_description()
        if filter_description:
            name_parts.append(f"({filter_description})")

        return "\n".join(name_parts)

    def get_css_classes(self):
        css_classes = super().get_css_classes()
        if self.trio:
            if self.trio.mother_affected:
                css_classes.append("mother-affected")
            if self.trio.father_affected:
                css_classes.append("father-affected")
        return css_classes

    @staticmethod
    def get_node_class_label():
        return 'Trio'

    def _get_configuration_errors(self):
        errors = []
        if not self.trio:
            errors.append("No trio selected")
        return errors

    def _get_cohorts_and_sample_visibility_for_node(self):
        cohorts = []
        visibility = {}
        if self.trio:
            cohort = self.trio.cohort
            cohorts = [cohort]
            visibility = {s: cohort.has_genotype for s in self.trio.get_samples()}
        return cohorts, visibility

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        proband_sample = None
        if self.trio:
            proband_sample = self.trio.proband.sample
        return proband_sample

    def __str__(self):
        return f"TrioNode: {self.pk}"
