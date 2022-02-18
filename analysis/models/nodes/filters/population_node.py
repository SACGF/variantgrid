import operator
from functools import reduce
from typing import Optional, Dict

from django.conf import settings
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q
from lazy import lazy

from analysis.models.enums import GroupOperation
from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.models.models import VariantAnnotation
from classification.enums import ClinicalSignificance
from classification.models.classification import Classification
from patients.models_enums import SimpleZygosity, GnomADPopulation
from snpdb.models import Sample, VariantZygosityCountCollection


class PopulationNode(AnalysisNode):
    EVERYTHING = 100.0  # Percent
    percent = models.FloatField(default=EVERYTHING)
    group_operation = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ANY)
    # highest in gnomAD - for diff groups see PopulationNodeGnomADPopulation below
    gnomad_af = models.BooleanField(default=True)
    gnomad_popmax_af = models.BooleanField(default=False)
    gnomad_hom_alt_max = models.IntegerField(null=True, blank=True)
    show_gnomad_filtered = models.BooleanField(default=True)
    af_1kg = models.BooleanField(default=True)
    af_uk10k = models.BooleanField(default=True)
    topmed_af = models.BooleanField(default=False)
    zygosity = models.CharField(max_length=1, choices=SimpleZygosity.choices, default=SimpleZygosity.ANY_GERMLINE)
    use_internal_counts = models.BooleanField(default=False)
    max_samples = models.IntegerField(null=True, blank=True)
    internal_percent = models.FloatField(default=EVERYTHING)
    keep_internally_classified_pathogenic = models.BooleanField(default=True)

    POPULATION_DATABASE_FIELDS = ["gnomad_af", "gnomad_popmax_af", "af_1kg", "af_uk10k", "topmed_af"]

    @lazy
    def num_samples_for_build(self) -> int:
        return Sample.objects.filter(vcf__genome_build=self.analysis.genome_build).count()

    @property
    def filtering_by_population(self):
        return self.percent != self.EVERYTHING

    def modifies_parents(self):
        return any([self.filtering_by_population, self.use_internal_counts,
                    self.gnomad_hom_alt_max, not self.show_gnomad_filtered])

    def _get_kwargs_for_parent_annotation_kwargs(self):
        kwargs = {}
        if (self.gnomad_af or self.gnomad_popmax_af) and self.percent <= 5:
            kwargs["common_variants"] = False
        return kwargs

    def _get_annotation_kwargs_for_node(self, **kwargs) -> Dict:
        annotation_kwargs = super()._get_annotation_kwargs_for_node(**kwargs)
        if self.use_internal_counts:
            vzcc = VariantZygosityCountCollection.objects.get(name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)
            annotation_kwargs.update(vzcc.get_annotation_kwargs(**kwargs))
        return annotation_kwargs

    def _get_node_q(self) -> Optional[Q]:
        and_q = []
        if self.filtering_by_population:
            population_databases = set()
            for field in self.POPULATION_DATABASE_FIELDS:
                if getattr(self, field):
                    population_databases.add(field)

            for gnomad_pop in self.populationnodegnomadpopulation_set.all():
                field = VariantAnnotation.get_gnomad_population_field(gnomad_pop.population)
                population_databases.add(field)

            if population_databases:
                # The group operation is backwards from what you may expect, as the widget takes MAX
                # but we are filtering for NULL or less than
                # ANY - remove if frequency is above cutoff in any database (much more strict).
                # ALL - Only remove if frequency is above cutoff in all databases.
                # Note: ALL retains variants with a missing entry from any database
                OPERATIONS = {  # These are opposite of normal GroupOperation
                    GroupOperation.ALL: operator.or_,
                    GroupOperation.ANY: operator.and_,
                }
                group_operation = OPERATIONS[GroupOperation(self.group_operation)]
                max_allele_frequency = self.percent / 100
                filters = []
                for field in population_databases:
                    q_isnull = Q(**{f"variantannotation__{field}__isnull": True})
                    q_max_value = Q(**{f"variantannotation__{field}__lte": max_allele_frequency})
                    filters.append(q_isnull | q_max_value)

                and_q.append(reduce(group_operation, filters))

        if self.gnomad_hom_alt_max is not None:
            q_hom_alt_lt = Q(variantannotation__gnomad_hom_alt__lte=self.gnomad_hom_alt_max)
            q_hom_alt_null = Q(variantannotation__gnomad_hom_alt__isnull=True)
            and_q.append(q_hom_alt_lt | q_hom_alt_null)

        if not self.show_gnomad_filtered:
            q_gnomad_filtered = Q(variantannotation__gnomad_filtered=False)
            q_gnomad_filtered_null = Q(variantannotation__gnomad_filtered__isnull=True)
            and_q.append(q_gnomad_filtered | q_gnomad_filtered_null)

        # Internal filters
        if self.use_internal_counts:
            # Global counts is added to every get_annotation_kwargs()
            vzcc = VariantZygosityCountCollection.objects.get(name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)

            if self.internal_percent != self.EVERYTHING:
                max_samples_from_percent = int(self.num_samples_for_build * (self.internal_percent / 100))
                # Min of 1 - so don't filter out self if only copy!
                max_samples_from_percent = max(max_samples_from_percent, 1)
                less_than = Q(**{vzcc.germline_counts_alias + "__lte": max_samples_from_percent})
                is_null = Q(**{vzcc.germline_counts_alias + "__isnull": True})
                and_q.append(less_than | is_null)

            if self.max_samples is not None:
                if self.zygosity == SimpleZygosity.ANY_GERMLINE:
                    column = vzcc.germline_counts_alias
                elif self.zygosity == SimpleZygosity.ANY_ZYGOSITY:
                    column = vzcc.all_zygosity_counts_alias
                elif self.zygosity == SimpleZygosity.HET:
                    column = vzcc.het_alias
                elif self.zygosity == SimpleZygosity.HOM_ALT:
                    column = vzcc.hom_alias
                else:
                    msg = f"Unknown SimpleZygosity {self.zygosity}"
                    raise ValueError(msg)

                less_than = Q(**{column + "__lte": self.max_samples})
                is_null = Q(**{column + "__isnull": True})

                and_q.append(less_than | is_null)

        if and_q:
            or_q = [reduce(operator.and_, and_q)]
            if self.keep_internally_classified_pathogenic:
                path_and_likely_path = [ClinicalSignificance.LIKELY_PATHOGENIC, ClinicalSignificance.PATHOGENIC]
                q_classified = Classification.get_variant_q(self.analysis.user, self.analysis.genome_build,
                                                            clinical_significance_list=path_and_likely_path)
                or_q.append(q_classified)

            q_node = reduce(operator.or_, or_q)
            q = reduce(operator.and_, [q_node])
        else:
            q = None
        return q

    def _get_method_summary(self):
        if self.modifies_parents():
            filters = []
            max_allele_frequency = self.percent / 100
            for field in self.POPULATION_DATABASE_FIELDS:
                if getattr(self, field):
                    filters.append('%s <= %g' % (field, max_allele_frequency))

            for gnomad_pop in self.populationnodegnomadpopulation_set.all():
                field = VariantAnnotation.get_gnomad_population_field(gnomad_pop.population)
                filters.append('%s <= %g' % (field, max_allele_frequency))

            if self.gnomad_hom_alt_max is not None:
                filters.append(f"gnomad_hom_alt_max <= {self.gnomad_hom_alt_max}")

            if not self.show_gnomad_filtered:
                filters.append("Removing gnomAD filtered")

            method_summary = ', '.join(filters)
        else:
            method_summary = "Set to 100% - No filtering applied"

        return method_summary

    def get_node_name(self):
        ops = []
        if self.filtering_by_population:
            ops.append("<= %g%%" % self.percent)
        if self.gnomad_hom_alt_max is not None:
            ops.append(f"<= {self.gnomad_hom_alt_max} gnomad hom alt")

        if self.use_internal_counts:
            internal_filter_ops = []
            if self.internal_percent != self.EVERYTHING:
                internal_filter_ops.append("%g%%" % self.internal_percent)
            if self.max_samples is not None:
                zygosity = self.get_zygosity_display()
                max_samples_desc = f"<= {self.max_samples} samples ({zygosity})"
                internal_filter_ops.append(max_samples_desc)
            if internal_filter_ops:
                ops.append("Internal DB:")
                ops.extend(internal_filter_ops)

        if ops:
            name = '\n'.join(ops)
        else:
            name = ''
        return name

    @staticmethod
    def get_help_text() -> str:
        return "Filter on population frequencies in public databases (gnomAD/1KG/UK10K/TopMed) " \
               "or number of samples in this database"

    def save_clone(self):
        gnomad_pops = list(self.populationnodegnomadpopulation_set.all())
        copy = super().save_clone()
        for g_pop in gnomad_pops:
            copy.populationnodegnomadpopulation_set.create(population=g_pop.population)
        return copy

    @staticmethod
    def get_node_class_label():
        return "Population"


class PopulationNodeGnomADPopulation(models.Model):
    population_node = models.ForeignKey(PopulationNode, on_delete=CASCADE)
    population = models.CharField(max_length=3, choices=GnomADPopulation.choices)
