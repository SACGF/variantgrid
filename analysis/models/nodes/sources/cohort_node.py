from typing import Dict, Optional, List

from django.db import models
from django.db.models import Q
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.expressions import Value, F
from django.db.models.functions import Concat, Substr, Length, Replace
from django.urls.base import reverse
from functools import reduce
import operator

from analysis.models import GroupOperation, AnalysisNode
from analysis.models.nodes.cohort_mixin import CohortMixin
from analysis.models.nodes.zygosity_count_node import AbstractZygosityCountNode
from patients.models_enums import Zygosity, SimpleZygosity
from snpdb.models import Cohort, CohortSample, VariantsType


class AbstractCohortBasedNode(CohortMixin, AnalysisNode):
    min_ad = models.IntegerField(default=0)
    min_dp = models.IntegerField(default=0)
    min_gq = models.IntegerField(default=0)
    max_pl = models.IntegerField(null=True, blank=True)
    min_ad_op = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ANY)
    min_dp_op = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ANY)
    min_gq_op = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ANY)
    max_pl_op = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ANY)

    COHORT_GENOTYPE_FIELD_MAPPINGS = [("min_ad", "samples_allele_depth", "gte"),
                                      ("min_dp", "samples_read_depth", "gte"),
                                      ("min_gq", "samples_genotype_quality", "gte"),
                                      ("max_pl", "samples_phred_likelihood", "lte")]

    class Meta:
        abstract = True

    def _get_q_and_list(self) -> List[Q]:
        q_and = super()._get_q_and_list()

        cohort = self._get_cohort()
        cgc = cohort.cohort_genotype_collection

        array_indicies = [cgc.get_array_index_for_sample_id(sample_id) for sample_id in cohort.get_sample_ids()]
        for field, cg_path, q_op in self.COHORT_GENOTYPE_FIELD_MAPPINGS:
            value = getattr(self, field)
            if value:
                filters = []
                for i in array_indicies:
                    filters.append(Q(**{f"{cgc.cohortgenotype_alias}__{cg_path}__{i}__{q_op}": value}))
                group_operation = getattr(self, field + "_op")
                q_and.append(GroupOperation.reduce(filters, group_operation))
        return q_and


class CohortNode(AbstractCohortBasedNode, AbstractZygosityCountNode):
    COUNT = 0
    SIMPLE_ZYGOSITY = 1
    PER_SAMPLE_ZYGOSITY = 2

    cohort = models.ForeignKey(Cohort, null=True, on_delete=SET_NULL)
    zygosity = models.CharField(max_length=1, choices=SimpleZygosity.choices, default=SimpleZygosity.ANY_GERMLINE)
    zygosity_op = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ALL)
    accordion_panel = models.IntegerField(default=0)
    # /zygosity count

    min_inputs = 0
    max_inputs = 0

    def _get_cohort(self):
        return self.cohort

    def _get_cohorts_and_sample_visibility_for_node(self):
        cohorts, visibility = super()._get_cohorts_and_sample_visibility()
        if self.accordion_panel == self.PER_SAMPLE_ZYGOSITY:
            try:
                cnzf_collection = CohortNodeZygosityFiltersCollection.objects.get(cohort_node=self, cohort=self.cohort)
                for cnzf in cnzf_collection.cohortnodezygosityfilter_set.all():
                    sample = cnzf.cohort_sample.sample
                    visibility[sample] = cnzf.show_in_grid and visibility[sample]
            except CohortNodeZygosityFiltersCollection.DoesNotExist:
                pass
        return cohorts, visibility

    def _get_annotation_kwargs_for_node(self) -> Dict:
        annotation_kwargs = super()._get_annotation_kwargs_for_node()

        if self.cohort:
            # We need to join to our cohort genotype before annotate, or the counts etc will be for the whole table
            if self.cohort.is_sub_cohort():
                cgc = self.cohort.cohort_genotype_collection
                sample_substrings = []
                for sample in self.cohort.get_samples():
                    i = cgc.get_sql_index_for_sample_id(sample.pk)
                    sample_substrings.append(Substr(f"{cgc.cohortgenotype_alias}__samples_zygosity", i, length=1))

                sub_cohort_zygosity = Concat(*sample_substrings)
                remove_hom = Replace(sub_cohort_zygosity, Value(Zygosity.HOM_ALT), Value(''))
                remove_het = Replace(sub_cohort_zygosity, Value(Zygosity.HET), Value(''))
                remove_ref = Replace(sub_cohort_zygosity, Value(Zygosity.HOM_REF), Value(''))

                hom_count = Length(sub_cohort_zygosity) - Length(remove_hom)
                het_count = Length(sub_cohort_zygosity) - Length(remove_het)
                ref_count = Length(sub_cohort_zygosity) - Length(remove_ref)
                annotation_kwargs[self.hom_count_column] = hom_count
                annotation_kwargs[self.het_count_column] = het_count
                annotation_kwargs[self.ref_count_column] = ref_count

            # Just add all annotations (only those used will be actually executed)
            hom_and_het = F(self.hom_count_column) + F(self.het_count_column)
            annotation_kwargs[self.any_germline_count_column] = hom_and_het
            annotation_kwargs[self.any_zygosity_count_column] = hom_and_het + F(self.ref_count_column)
        return annotation_kwargs

    def _get_node_q(self) -> Optional[Q]:
        cohort, q = self.get_cohort_and_q()
        if cohort:
            q &= self.get_cohort_settings_q(cohort)
            q &= self.get_vcf_locus_filters_q()
        return q

    def get_cohort_settings_q(self, cohort):
        cohort_genotype_collection = cohort.cohort_genotype_collection
        and_q = []
        if self.accordion_panel == self.COUNT:
            # Use minimum filters even for sub-cohorts as the min will always be above sub-cohort min
            # And they are a fast integer test, so may be a quicker path
            and_q.extend(self.get_zygosity_count_q_list())
        elif self.accordion_panel == self.SIMPLE_ZYGOSITY:
            and_q.extend(self.get_cohort_simple_zygosity_q_list())
        elif self.accordion_panel == self.PER_SAMPLE_ZYGOSITY:
            and_q.extend(self.get_cohort_per_sample_zygosity_q_list(cohort_genotype_collection))

        if and_q:
            q = reduce(operator.and_, and_q)
        else:
            q = self.q_all()
        return q

    def get_cohort_simple_zygosity_q_list(self):
        COLUMNS = {
            SimpleZygosity.REF: self.ref_count_column,
            SimpleZygosity.HET: self.het_count_column,
            SimpleZygosity.HOM_ALT: self.hom_count_column,
            SimpleZygosity.ANY_GERMLINE: self.any_germline_count_column,
            SimpleZygosity.ANY_ZYGOSITY: self.any_zygosity_count_column,
        }
        q_and = []
        if self.zygosity:
            column = COLUMNS[self.zygosity]
            min_count = None
            if self.zygosity_op == GroupOperation.ALL:
                min_count = len(self.cohort.get_samples())
            elif self.zygosity_op == GroupOperation.ANY:
                min_count = 1

            if column and min_count:
                q_and = [Q(**{column + "__gte": min_count})]
        return q_and

    def get_cohort_per_sample_zygosity_q_list(self, cohort_genotype_collection):
        sample_zygosities_dict = {}

        # Start with zygosity_sample_matches of all dots (as need to handle all
        # samples in cohort not just those in sub cohort)
        cnzf_collection = CohortNodeZygosityFiltersCollection.objects.get(cohort_node=self, cohort=self.cohort)
        ZYGOSITY_LOOKUP = {'zygosity_ref': [Zygosity.HOM_REF],
                           'zygosity_het': [Zygosity.HET],
                           'zygosity_hom': [Zygosity.HOM_ALT],
                           'zygosity_none': [Zygosity.MISSING, Zygosity.UNKNOWN_ZYGOSITY]}
        for cnzf in cnzf_collection.cohortnodezygosityfilter_set.all():
            zygosities = set()
            for f, zyg in ZYGOSITY_LOOKUP.items():
                if getattr(cnzf, f):
                    zygosities.update(zyg)

            sample_zygosities_dict[cnzf.cohort_sample.sample] = zygosities

        return [cohort_genotype_collection.get_zygosity_q(sample_zygosities_dict)]

    def _get_method_summary(self):
        return f"Cohort ({self.cohort.name})"

    def get_node_name(self):
        name = ''
        if self.cohort:
            num_samples = self.cohort.cohortsample_set.count()
            if self.minimum_count:
                num_message = f"({self.minimum_count} of {num_samples}"
            else:
                num_message = f"({num_samples})"

            if self.zygosity:
                zygosity_message = f" ({self.get_zygosity_display()})"
            else:
                zygosity_message = ''
            name = f"{self.cohort.name}{zygosity_message}\n{num_message}"

            filter_description = self.get_filter_description()
            if filter_description:
                name += f"\n({filter_description})"
        return name

    @staticmethod
    def get_node_class_label():
        return "Cohort"

    def _get_configuration_errors(self) -> List:
        errors = super()._get_configuration_errors()
        if not self.cohort:
            errors.append("No cohort selected.")
        else:
            errors.extend(self._get_genome_build_errors("cohort", self.cohort.genome_build))
            try:
                self.cohort.cohort_genotype_collection
            except:
                cohort_name = self.cohort.name or "cohort"
                url = reverse("view_cohort", kwargs={"cohort_id": self.cohort.pk})
                msg = "Your cohort has not finished processing - please visit cohort page for "
                msg += f"<a href='{url}'>{cohort_name}</a>"
                errors.append(msg)

        return errors

    def _get_node_extra_columns(self):
        extra_columns = super()._get_node_extra_columns()

#        extra_columns.append(self.hom_count_column)
#        extra_columns.append(self.het_count_column)
        return extra_columns

    def _get_node_extra_colmodel_overrides(self):
        extra_colmodel_overrides = super()._get_node_extra_colmodel_overrides()
        if self.cohort and self.cohort.is_sub_cohort():
            labels = ["hom_count", "het_count"]
            for c, l in zip([self.hom_count_column, self.het_count_column], labels):
                override = extra_colmodel_overrides.get(c, {})
                override["label"] = l
                override["name"] = c
                override["model_fields"] = False
                override["queryset_fields"] = True
                extra_colmodel_overrides[c] = override

        return extra_colmodel_overrides

    def save(self, **kwargs):
        is_new = self.version == 0
        vcf = self._get_vcf()

        if is_new and vcf:
            somatic_samples = vcf.sample_set.filter(variants_type__in=VariantsType.SOMATIC_TYPES)
            if somatic_samples.exists():
                # Somatic mode...
                self.accordion_panel = self.SIMPLE_ZYGOSITY

        return super().save(**kwargs)

    def __str__(self):
        if self.cohort:
            name = self.cohort.name
        else:
            name = ''

        return name


class CohortNodeZygosityFiltersCollection(models.Model):
    """ You may want multiple records for a node at the same time, as
        you can start editing before saving, then perhaps switch back etc. """
    cohort_node = models.ForeignKey(CohortNode, on_delete=CASCADE)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    cohort_version = models.IntegerField()

    class Meta:
        unique_together = ('cohort_node', 'cohort')

    @staticmethod
    def get_for_node_and_cohort(node, cohort):
        defaults = {"cohort_version": cohort.version}
        cnzfc, created = CohortNodeZygosityFiltersCollection.objects.get_or_create(cohort_node=node, cohort=cohort,
                                                                                   defaults=defaults)

        if created or cnzfc.cohort_version != cohort.version:
            # Deleted CohortSample will cascade, but we need to create new entries for added samples
            for cohort_sample in cohort.cohortsample_set.all().order_by("sample"):
                defaults = {"show_in_grid": not cohort_sample.sample.no_dna_control}
                CohortNodeZygosityFilter.objects.get_or_create(collection=cnzfc,
                                                               cohort_sample=cohort_sample,
                                                               defaults=defaults)
            if cnzfc.cohort_version != cohort.version:
                cnzfc.cohort_version = cohort.version
                cnzfc.save()
        return cnzfc

    def __str__(self):
        return f"CohortNodeZygosityFiltersCollection: {self.cohort_node}, {self.cohort}"


class CohortNodeZygosityFilter(models.Model):
    collection = models.ForeignKey(CohortNodeZygosityFiltersCollection, on_delete=CASCADE)
    cohort_sample = models.ForeignKey(CohortSample, on_delete=CASCADE)
    show_in_grid = models.BooleanField(default=True)
    zygosity_ref = models.BooleanField(default=True)
    zygosity_het = models.BooleanField(default=True)
    zygosity_hom = models.BooleanField(default=True)
    zygosity_none = models.BooleanField(default=True)

    unique_together = ("collection", "cohort_sample")
