import time
from abc import ABC, abstractmethod
from typing import Optional, List, Set, Tuple

from cache_memoize import cache_memoize
from django.db import models
from django.db.models import Count
from django.db.models.deletion import SET_NULL
from django.db.models.query_utils import Q

from analysis.models.enums import TrioInheritance, NodeErrorSource, AnalysisTemplateType
from analysis.models.nodes.sources import AbstractCohortBasedNode
from annotation.models.models import VariantTranscriptAnnotation
from library.constants import DAY_SECS
from patients.models_enums import Zygosity, Sex
from snpdb.models import Trio, Sample, Contig


class AbstractTrioInheritance(ABC):
    """ Do inheritance filtering in subclasses to keep filters/methods consistent """
    NO_VARIANT = {Zygosity.MISSING, Zygosity.HOM_REF}  # 2 different het would be "missing" (as has no ref)
    HAS_VARIANT = {Zygosity.HET, Zygosity.HOM_ALT}
    UNAFFECTED_AND_AFFECTED_ZYGOSITIES = [NO_VARIANT, HAS_VARIANT]

    def __init__(self, node: 'TrioNode'):
        self.node = node

    def _get_zyg_q(self, cohort_genotype_collection, trio_zyg_data) -> Q:
        """ trio_zyg_data = tuple of mum_zyg_set, dad_zyg_set, proband_zyg_set """
        mother_sample = self.node.trio.mother.sample
        father_sample = self.node.trio.father.sample
        proband_sample = self.node.trio.proband.sample

        mum_dad_proband = [mother_sample, father_sample, proband_sample]
        sample_zygosities_dict = dict(zip(mum_dad_proband, trio_zyg_data))
        sample_require_zygosity_dict = {mother_sample: self.node.require_zygosity,
                                        father_sample: self.node.require_zygosity,
                                        proband_sample: True}  # 947 - Always require zygosity for Proband
        return cohort_genotype_collection.get_zygosity_q(sample_zygosities_dict, sample_require_zygosity_dict)

    @staticmethod
    def _zygosity_options(zyg: Set, allow_unknown=False):
        if allow_unknown:
            # Make a new set so as to not alter passed in value
            zyg = zyg | {Zygosity.UNKNOWN_ZYGOSITY}
        zyg = zyg - {Zygosity.MISSING}  # Implementation detail - don't show to user
        return ", ".join(sorted([Zygosity.display(z) for z in zyg]))

    def get_zygosities_method(self, mum_z: Set, dad_z: Set, proband_z: Set):
        proband = self._zygosity_options(proband_z)
        mum = self._zygosity_options(mum_z, not self.node.require_zygosity)
        dad = self._zygosity_options(dad_z, not self.node.require_zygosity)
        filters = {"Proband": proband, "Mother": mum, "Father": dad}
        return ", ".join([f"{k}: {v}" for k, v in filters.items() if v])

    @abstractmethod
    def get_q(self) -> Q:
        pass

    @abstractmethod
    def get_method(self) -> str:
        pass

    def get_contigs(self) -> Optional[Set[Contig]]:
        """ None means we don't know """
        return None


class SimpleTrioInheritance(AbstractTrioInheritance):
    @abstractmethod
    def _get_mum_dad_proband_zygosities(self) -> Tuple[Set, Set, Set]:
        pass

    def get_q(self) -> Q:
        cohort_genotype_collection = self.node.trio.cohort.cohort_genotype_collection
        return self._get_zyg_q(cohort_genotype_collection, self._get_mum_dad_proband_zygosities())

    def get_method(self) -> str:
        return self.get_zygosities_method(*self._get_mum_dad_proband_zygosities())


class Recessive(SimpleTrioInheritance):
    def _get_mum_dad_proband_zygosities(self) -> Tuple[Set, Set, Set]:
        return {Zygosity.HET}, {Zygosity.HET}, {Zygosity.HOM_ALT}


class Dominant(SimpleTrioInheritance):
    def _get_mum_dad_proband_zygosities(self) -> Tuple[Set, Set, Set]:
        mother_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(self.node.trio.mother_affected)]
        father_zyg = self.UNAFFECTED_AND_AFFECTED_ZYGOSITIES[int(self.node.trio.father_affected)]
        return mother_zyg, father_zyg, self.HAS_VARIANT


class Denovo(SimpleTrioInheritance):
    def _get_mum_dad_proband_zygosities(self) -> Tuple[Set, Set, Set]:
        return self.NO_VARIANT, self.NO_VARIANT, self.HAS_VARIANT


class ProbandHet(SimpleTrioInheritance):
    """ Equivalent of SampleNode, but in TrioNode so analysis templates only need to have Trio as argument """
    def _get_mum_dad_proband_zygosities(self) -> Tuple[Set, Set, Set]:
        return set(), set(), {Zygosity.HET}


class XLinkedRecessive(SimpleTrioInheritance):
    def _get_mum_dad_proband_zygosities(self) -> Tuple[Set, Set, Set]:
        return {Zygosity.HET}, set(), {Zygosity.HOM_ALT}

    def get_q(self) -> Q:
        return super().get_q() & Q(locus__contig__name='X')  # will work for hg19 and GRCh38

    def get_method(self) -> str:
        return super().get_method() + " and contig name = 'X'"

    def get_contigs(self) -> Optional[Set[Contig]]:
        return set(self.node.trio.genome_build.contigs.filter(name='X'))


class CompHet(AbstractTrioInheritance):
    def _mum_but_not_dad(self):
        return {Zygosity.HET}, self.NO_VARIANT, {Zygosity.HET}

    def _dad_but_not_mum(self):
        return self.NO_VARIANT, {Zygosity.HET}, {Zygosity.HET}

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.node.pk, s.node.version))
    def _get_parent_comp_het_q_and_two_hit_genes(self):
        cohort_genotype_collection = self.node.trio.cohort.cohort_genotype_collection

        parent = self.node.get_single_parent()
        mum_but_not_dad = self._get_zyg_q(cohort_genotype_collection, self._mum_but_not_dad())
        dad_but_not_mum = self._get_zyg_q(cohort_genotype_collection, self._dad_but_not_mum())
        comp_het_q = mum_but_not_dad | dad_but_not_mum

        # Need to pass in kwargs in case we have parent (eg Venn node) that doesn't have same cohort annotation kwargs
        annotation_kwargs = self.node.get_annotation_kwargs()

        def get_parent_genes(q):
            qs = parent.get_queryset(q, extra_annotation_kwargs=annotation_kwargs)
            return qs.values_list("varianttranscriptannotation__gene", flat=True).distinct()

        # This ends up doing 3 queries (where we call set() - to work out what Q we need to return)
        common_genes = set(get_parent_genes(mum_but_not_dad)) & set(get_parent_genes(dad_but_not_mum))
        q_in_genes = Q(varianttranscriptannotation__gene__in=common_genes)
        parent_genes_qs = parent.get_queryset(q_in_genes, extra_annotation_kwargs=annotation_kwargs)
        parent_genes_qs = parent_genes_qs.values_list("varianttranscriptannotation__gene")
        two_hits = parent_genes_qs.annotate(gene_count=Count("pk")).filter(gene_count__gte=2)
        two_hit_genes = set(two_hits.values_list("varianttranscriptannotation__gene", flat=True).distinct())
        return parent, comp_het_q, two_hit_genes

    def get_q(self) -> Q:
        parent, comp_het_q, two_hit_genes = self._get_parent_comp_het_q_and_two_hit_genes()
        comp_het_genes = VariantTranscriptAnnotation.get_overlapping_genes_q(two_hit_genes)
        return parent.get_q() & comp_het_q & comp_het_genes

    def get_method(self) -> str:
        mum1, dad1, _ = self._mum_but_not_dad()
        mum_but_not_dad = self.get_zygosities_method(mum1, dad1, set())
        mum2, dad2, _ = self._dad_but_not_mum()
        dad_but_not_mum = self.get_zygosities_method(mum2, dad2, set())
        return f"Proband: HET, and >=2 hits from genes where ({mum_but_not_dad}) OR ({dad_but_not_mum})"

    def get_contigs(self) -> Optional[Set[Contig]]:
        _, _, two_hit_genes = self._get_parent_comp_het_q_and_two_hit_genes()
        contig_qs = Contig.objects.filter(transcriptversion__genome_build=self.node.trio.genome_build,
                                          transcriptversion__gene_version__gene__in=two_hit_genes)
        return set(contig_qs.distinct())


class TrioNode(AbstractCohortBasedNode):
    trio = models.ForeignKey(Trio, null=True, on_delete=SET_NULL)
    inheritance = models.CharField(max_length=1, choices=TrioInheritance.choices, default=TrioInheritance.RECESSIVE)
    require_zygosity = models.BooleanField(default=True)

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
                if not (trio.mother_affected or trio.father_affected):
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

    def _inheritance_factory(self):
        inhertiance_classes = {
            TrioInheritance.COMPOUND_HET: CompHet,
            TrioInheritance.RECESSIVE: Recessive,
            TrioInheritance.DOMINANT: Dominant,
            TrioInheritance.DENOVO: Denovo,
            TrioInheritance.PROBAND_HET: ProbandHet,
            TrioInheritance.XLINKED_RECESSIVE: XLinkedRecessive
        }
        klass = inhertiance_classes[TrioInheritance(self.inheritance)]
        return klass(self)

    def _get_node_q(self) -> Optional[Q]:
        cohort, q = self.get_cohort_and_q()
        if cohort:
            inheritance = self._inheritance_factory()
            q &= inheritance.get_q()
            q &= self.get_vcf_locus_filters_q()
        return q

    def _get_node_contigs(self) -> Optional[Set[Contig]]:
        node_contigs = None
        if self.trio:
            inheritance = self._inheritance_factory()
            node_contigs = inheritance.get_contigs()
        return node_contigs

    def _get_method_summary(self):
        if self._get_cohort():
            inheritance = self._inheritance_factory()
            method = inheritance.get_method()
        else:
            method = "No cohort selected"
        return method

    def get_node_name(self):
        name_parts = [TrioInheritance(self.inheritance).label]
        if not self.require_zygosity:
            name_parts.append(' (non strict)')

        filter_description = self.get_filter_description()
        if filter_description:
            name_parts.append(f"({filter_description})")

        return "\n".join(name_parts)

    @staticmethod
    def get_help_text() -> str:
        return "Mother/Father/Proband - filter for recessive/dominant/denovo inheritance"

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

    def _get_configuration_errors(self) -> List:
        errors = super()._get_configuration_errors()
        if not self.trio:
            errors.append("No trio selected")
        else:
            errors.extend(self._get_genome_build_errors("trio", self.trio.genome_build))
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
