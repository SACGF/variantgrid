import operator
from functools import reduce
from typing import Optional

from auditlog.registry import auditlog
from cache_memoize import cache_memoize
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q

from analysis.models.enums import GroupOperation
from analysis.models.nodes.analysis_node import AnalysisNode, NodeAuditLogMixin
from annotation.models.models import VariantAnnotation
from classification.enums import ClinicalSignificance
from classification.models.classification import Classification
from library.constants import MINUTE_SECS
from patients.models_enums import SimpleZygosity, GnomADPopulation
from snpdb.models import VariantZygosityCountCollection


class PopulationNode(AnalysisNode):
    EVERYTHING = 100.0  # Percent
    percent = models.FloatField(default=EVERYTHING)
    group_operation = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ANY)
    # highest in gnomAD - for diff groups see PopulationNodeGnomADPopulation below
    gnomad_af = models.BooleanField(default=True)
    gnomad_popmax_af = models.BooleanField(default=False)
    gnomad_fafmax_faf95_max = models.BooleanField(default=False)
    gnomad_fafmax_faf99_max = models.BooleanField(default=False)
    gnomad_hom_alt_max = models.IntegerField(null=True, blank=True)
    show_gnomad_filtered = models.BooleanField(default=True)
    af_1kg = models.BooleanField(default=True)
    af_uk10k = models.BooleanField(default=True)
    topmed_af = models.BooleanField(default=False)
    zygosity = models.CharField(max_length=1, choices=SimpleZygosity.choices, default=SimpleZygosity.NON_REF_CALL)
    use_internal_counts = models.BooleanField(default=False)
    max_samples = models.IntegerField(null=True, blank=True)
    internal_percent = models.FloatField(default=EVERYTHING)
    keep_internally_classified_pathogenic = models.BooleanField(default=True)

    @property
    def filtering_by_population(self):
        return self.percent != self.EVERYTHING

    @property
    def population_database_fields(self) -> list[str]:
        fields = ["gnomad_af", "gnomad_popmax_af", "af_1kg", "af_uk10k", "topmed_af"]
        if self.has_filtering_allele_frequency:
            fields += ["gnomad_fafmax_faf95_max", "gnomad_fafmax_faf99_max"]
        return fields

    @property
    def has_filtering_allele_frequency(self) -> bool:
        try:
            return self.analysis.annotation_version.variant_annotation_version.gnomad_major_version >= 4
        except AttributeError as e:
            return False

    def modifies_parents(self):
        return any([self.filtering_by_population, self.use_internal_counts,
                    self.gnomad_hom_alt_max, not self.show_gnomad_filtered])

    def _has_common_variants(self) -> bool:
        rare_only = (self.gnomad_af or self.gnomad_popmax_af) and self.percent <= 5
        return not rare_only

    def _get_kwargs_for_parent_annotation_kwargs(self, common_variants=None, **kwargs):
        node_kwargs = {}
        if not (common_variants or self._has_common_variants()):
            node_kwargs["common_variants"] = False
        return node_kwargs

    def _get_annotation_kwargs_for_node(self, **kwargs) -> dict:
        annotation_kwargs = super()._get_annotation_kwargs_for_node(**kwargs)
        if self.use_internal_counts:
            vzcc = VariantZygosityCountCollection.get_global_germline_counts()
            annotation_kwargs.update(vzcc.get_annotation_kwargs(**kwargs))
        return annotation_kwargs

    def _get_node_q(self) -> Optional[Q]:
        and_q = []
        if self.filtering_by_population:
            population_databases = []
            for field in self.population_database_fields:
                if getattr(self, field):
                    population_databases.append(field)

            for gnomad_pop in self.populationnodegnomadpopulation_set.all():
                field = VariantAnnotation.get_gnomad_population_field(gnomad_pop.population)
                population_databases.append(field)

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

        if and_q:
            q = reduce(operator.and_, and_q)
        else:
            q = None
        return q

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        node_arg_q_dict = {}

        and_q = []

        if self.use_internal_counts:
            # Internal count filters use filtered relations - so we ideally want to apply against them to reduce the
            # number of joins. However, if we have keep_internally_classified_pathogenic - we need to apply them with
            # an OR - so need to put into the standard node_arg_q_dict[None] later

            vzcc = VariantZygosityCountCollection.get_global_germline_counts()
            # Max germline can be applied either by percent or from max_samples and zygosity = NON_REF
            max_germline_list = []
            if self.max_samples is not None:
                if self.zygosity == SimpleZygosity.NON_REF_CALL:
                    max_germline_list.append(self.max_samples)
                else:
                    if self.zygosity == SimpleZygosity.ANY_CALL:
                        column = vzcc.any_call_counts_alias
                    elif self.zygosity == SimpleZygosity.HET:
                        column = vzcc.het_alias
                    elif self.zygosity == SimpleZygosity.HOM_ALT:
                        column = vzcc.hom_alias
                    else:
                        msg = f"Unknown SimpleZygosity {self.zygosity}"
                        raise ValueError(msg)

                    less_than = Q(**{column + "__lte": self.max_samples})
                    is_null = Q(**{column + "__isnull": True})

                    q = less_than | is_null
                    if self.keep_internally_classified_pathogenic:
                        and_q.append(q)
                    else:
                        node_arg_q_dict[column] = {str(q): q}

            if self.internal_percent != self.EVERYTHING:
                max_samples_from_percent = int(self.num_samples_for_build * (self.internal_percent / 100))
                # Min of 1 - so don't filter out self if only copy!
                max_samples_from_percent = max(max_samples_from_percent, 1)
                max_germline_list.append(max_samples_from_percent)

            if max_germline_list:
                max_germline = min(max_germline_list)
                less_than = Q(**{vzcc.non_ref_call_alias + "__lte": max_germline})
                is_null = Q(**{vzcc.non_ref_call_alias + "__isnull": True})
                q = less_than | is_null
                if self.keep_internally_classified_pathogenic:
                    and_q.append(q)
                else:
                    node_arg_q_dict[vzcc.non_ref_call_alias] = {str(q): q}

        if node_q := self._get_node_q():
            and_q.append(node_q)

        if and_q:
            or_q = [reduce(operator.and_, and_q)]
            if self.keep_internally_classified_pathogenic:
                parent = self.get_single_parent()
                classified_variant_ids = self._get_parent_classified_variant_ids(parent)
                or_q.append(Q(pk__in=classified_variant_ids))

            q = reduce(operator.or_, or_q)
            node_arg_q_dict[None] = {str(q): q}

        print(f"{node_arg_q_dict=}")
        return node_arg_q_dict

    @staticmethod
    @cache_memoize(5 * MINUTE_SECS, args_rewrite=lambda p: (p.pk, p.version))
    def _get_parent_classified_variant_ids(parent) -> list:
        """ Uses parent pk/version so can be shared with siblings. Doesn't use exclude as that ended up being slow """
        qs = parent.get_queryset()
        path_and_likely_path = [ClinicalSignificance.LIKELY_PATHOGENIC, ClinicalSignificance.PATHOGENIC]
        q_classified = Classification.get_variant_q(parent.analysis.user, parent.analysis.genome_build,
                                                    clinical_significance_list=path_and_likely_path)
        qs = qs.filter(q_classified)
        return list(qs.values_list("pk", flat=True))

    def _get_method_summary(self):
        if self.modifies_parents():
            filters = []
            max_allele_frequency = self.percent / 100
            for field in self.population_database_fields:
                if getattr(self, field):
                    filters.append(f'{field} <= {max_allele_frequency:g}')

            for gnomad_pop in self.populationnodegnomadpopulation_set.all():
                field = VariantAnnotation.get_gnomad_population_field(gnomad_pop.population)
                filters.append(f'{field} <= {max_allele_frequency:g}')

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
            ops.append(f"<= {self.percent:g}%")
        if self.gnomad_hom_alt_max is not None:
            ops.append(f"<= {self.gnomad_hom_alt_max} gnomad hom alt")

        if self.use_internal_counts:
            internal_filter_ops = []
            if self.internal_percent != self.EVERYTHING:
                internal_filter_ops.append(f"{self.internal_percent:g}%")
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


class PopulationNodeGnomADPopulation(NodeAuditLogMixin, models.Model):
    population_node = models.ForeignKey(PopulationNode, on_delete=CASCADE)
    population = models.CharField(max_length=3, choices=GnomADPopulation.choices)

    def _get_node(self):
        return self.population_node

    def __str__(self):
        return f"PopNode {self.population_node_id}: {self.get_population_display()}"


auditlog.register(PopulationNode)
auditlog.register(PopulationNodeGnomADPopulation)
