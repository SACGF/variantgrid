import operator
from functools import reduce
from typing import Optional, List

from django.core.exceptions import ObjectDoesNotExist
from django.db import models
from django.db.models import SET_NULL
from django.db.models.query_utils import Q

from analysis.models import GeneCoverageMixin
from analysis.models.nodes.analysis_node import AnalysisNode, NodeAlleleFrequencyFilter
from analysis.models.nodes.cohort_mixin import SampleMixin
from annotation.models import SampleClinVarAnnotationStats, SampleClinVarAnnotationStatsPassingFilter, \
    SampleGeneAnnotationStats, SampleGeneAnnotationStatsPassingFilter, \
    SampleVariantAnnotationStats, SampleVariantAnnotationStatsPassingFilter
from genes.models import SampleGeneList
from patients.models_enums import Zygosity
from snpdb.models import SampleStats, SampleStatsPassingFilter, Sample
from snpdb.models.models_enums import BuiltInFilters


class SampleNode(SampleMixin, GeneCoverageMixin, AnalysisNode):
    """ Use restrict_to_qc_gene_list to keep track of that as sample_gene_list is cleared when a sample changes
        including in an AnalysisTemplates """
    sample = models.ForeignKey(Sample, null=True, on_delete=SET_NULL)
    # When setting sample, if restrict_to_qc_gene_list = True, sample_gene_list is set to active sample gene list
    sample_gene_list = models.ForeignKey(SampleGeneList, null=True, blank=True, on_delete=SET_NULL)
    has_gene_coverage = models.BooleanField(null=True)
    min_ad = models.IntegerField(default=0)
    min_dp = models.IntegerField(default=0)
    min_gq = models.IntegerField(default=0)
    max_pl = models.IntegerField(null=True, blank=True)
    zygosity_ref = models.BooleanField(default=False)
    zygosity_het = models.BooleanField(default=True)
    zygosity_hom = models.BooleanField(default=True)
    zygosity_unk = models.BooleanField(default=False)
    restrict_to_qc_gene_list = models.BooleanField(default=False)

    FIELDS_THAT_CHANGE_QUERYSET = ("min_ad", "min_dp", "min_gq", "max_pl", "restrict_to_qc_gene_list")
    SAMPLE_FIELD_MAPPINGS = [("ad", "allele_depth"),
                             ("dp", "read_depth"),
                             ("gq", "genotype_quality")]
    min_inputs = 0
    max_inputs = 0

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        proband_sample = None
        if self.sample:
            proband_sample = self.sample
        return proband_sample

    def _get_zygosities(self):
        # No genotype - all will be unknown
        if self.sample and self.sample.has_genotype is False:
            return [Zygosity.UNKNOWN_ZYGOSITY]

        ZYGOSITY_LOOKUP = [
            (self.zygosity_ref, Zygosity.HOM_REF),
            (self.zygosity_het, Zygosity.HET),
            (self.zygosity_hom, Zygosity.HOM_ALT),
            (self.zygosity_unk, Zygosity.UNKNOWN_ZYGOSITY),
        ]
        return [z for s, z in ZYGOSITY_LOOKUP if s]

    def get_zygosity_description(self):
        zygosities = self._get_zygosities()
        if len(zygosities) != len(Zygosity.CHOICES):  # Subset
            zyg_dict = dict(Zygosity.CHOICES)
            zyg = ','.join([zyg_dict[z] for z in zygosities])
        else:
            zyg = ''  # Any
        return zyg

    def _get_zygosity_q(self):
        zygosity = self._get_zygosities()
        zygosity_q = None
        if zygosity:
            field = self.sample.get_cohort_genotype_field("zygosity")
            zygosity_q = Q(**{f"{field}__in": zygosity})

        return zygosity_q

    def _get_node_q(self) -> Optional[Q]:
        if not self.sample:
            return self.q_none()

        q_and = []
        # _get_zygosity_q handles no genotype (UNKNOWN)
        if zygosity_q := self._get_zygosity_q():
            q_and.append(zygosity_q)

        if self.sample.has_genotype:
            for node_field, ov_field in self.SAMPLE_FIELD_MAPPINGS:
                min_value = getattr(self, f"min_{node_field}")
                if min_value:
                    ov_path = self.sample.get_cohort_genotype_field(ov_field)
                    q_and.append(Q(**{f"{ov_path}__gte": min_value}))

            if self.max_pl is not None:
                pl_path = self.sample.get_cohort_genotype_field("phred_likelihood")
                q_and.append(Q(**{f"{pl_path}__lte": self.max_pl}))

            if af_q := NodeAlleleFrequencyFilter.get_sample_q(self, self.sample):
                q_and.append(af_q)

        if self.restrict_to_qc_gene_list:
            if self.sample_gene_list:
                q = self.sample_gene_list.gene_list.get_q(self.analysis.annotation_version.variant_annotation_version)
            else:
                q = self.q_none()  # Safety - don't show anything if missing
            q_and.append(q)

        q_and.append(self.get_vcf_locus_filters_q())
        return reduce(operator.and_, q_and)

    def _get_method_summary(self):
        if self.sample:
            sample_description = f"Sample {self.sample.name}"
            for node_field, _ in self.SAMPLE_FIELD_MAPPINGS:
                min_value = getattr(self, f"min_{node_field}")
                if min_value:
                    sample_description += f" {node_field.upper()}>={min_value}"
            if self.max_pl is not None:
                sample_description += f" PL<={self.max_pl}"

            if self.sample.has_genotype:
                zyg = self.get_zygosity_description()
                if zyg:
                    sample_description += f" ({zyg})"
            method_summary = f"{sample_description}, from VCF {self.sample.vcf}"
        else:
            method_summary = 'No sample selected'
        return method_summary

    def get_node_name(self):
        name_parts = []
        if self.sample:
            name_parts.append(self.sample.name)
            if self.sample.has_genotype:
                zyg = self.get_zygosity_description()
                if zyg:
                    name_parts.append(f"\n({zyg})")

            filter_description = self.get_filter_description()
            if filter_description:
                name_parts.append(f"\n({filter_description})")

        return "\n".join(name_parts)

    @staticmethod
    def get_help_text() -> str:
        return "Variants from a VCF sample, usually one genotype (patient, cell or organism)"

    def get_gene_lists(self) -> List:
        """ Used for gene coverage """
        gene_lists = []
        if self.sample_gene_list and self.restrict_to_qc_gene_list:
            gene_lists.append(self.sample_gene_list.gene_list)
        return gene_lists

    def _has_filters_that_affect_label_counts(self):
        # This returns None if no filters to apply
        if NodeAlleleFrequencyFilter.get_sample_q(self, self.sample):
            return True

        return any([getattr(self, f) for f in self.FIELDS_THAT_CHANGE_QUERYSET])

    def _get_cached_label_count(self, label):
        """ Input counts can be static, so use cached AnnotationStats if we can """
        if self._has_filters_that_affect_label_counts():
            return None  # Have to do counts

        CLASSES = {
            BuiltInFilters.TOTAL: (SampleStats, SampleStatsPassingFilter),
            BuiltInFilters.CLINVAR: (SampleClinVarAnnotationStats, SampleClinVarAnnotationStatsPassingFilter),
            BuiltInFilters.OMIM: (SampleGeneAnnotationStats, SampleGeneAnnotationStatsPassingFilter),
            BuiltInFilters.IMPACT_HIGH_OR_MODERATE: (SampleVariantAnnotationStats, SampleVariantAnnotationStatsPassingFilter),
        }
        annotation_version = self.analysis.annotation_version
        count = None
        try:
            filter_code = self.get_filter_code()
            klazz = CLASSES[label][filter_code]  # Maybe exception and out of range
            obj = klazz.load_version(self.sample, annotation_version)
            zygosities = [self.zygosity_ref, self.zygosity_het, self.zygosity_hom, self.zygosity_unk]
            if self.sample and not self.sample.has_genotype:
                zygosities = [True] * len(zygosities)  # Show everything

            count = obj.count_for_zygosity(*zygosities, label=label)
        except (IndexError, KeyError, ObjectDoesNotExist):
            pass  # OK, will just calculate it

        return count

    @staticmethod
    def get_node_class_label():
        return "Sample"

    def _get_configuration_errors(self) -> List:
        errors = super()._get_configuration_errors()
        if not self.sample:
            errors.append("No sample selected.")
        else:
            errors.extend(self._get_genome_build_errors("sample", self.sample.genome_build))
        if self.restrict_to_qc_gene_list and self.sample_gene_list is None:
            errors.append("Restricted to Sample Gene List, but none specified!")
        return errors

    def get_rendering_args(self):
        patient_args = {}
        if self.sample and self.sample.patient:
            patient_args = self.sample.patient.get_json_dict()

        return {"patient": patient_args}
