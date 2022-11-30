import operator
import re
from functools import reduce
from typing import List, Dict, Optional, Set, Tuple

from django.db.models import Q

from analysis.models.enums import GroupOperation
from analysis.models.nodes.analysis_node import NodeVCFFilter, NodeAlleleFrequencyFilter
from patients.models_enums import Zygosity
from snpdb.models import VCFFilter, Cohort
from upload.models import UploadedVCF


class CohortMixin:
    """ Since a Cohort is based off a VCF we also  """

    def _get_cohort(self):
        """ Each subclass needs to implement the way to get their Cohort """
        raise NotImplementedError()

    def _get_vcf(self):
        vcf = None
        if cohort := self._get_cohort():
            vcf = cohort.get_vcf()
        return vcf

    def _get_cache_key(self) -> str:
        """ Use cohort genotype in the key, as that can change if a VCF is reloaded """
        cache_key = super()._get_cache_key()
        cgc_id = 0
        if cgc := self.cohort_genotype_collection:
            cgc_id = cgc.pk
        return "_".join((cache_key, str(cgc_id)))

    @property
    def cohort_genotype_collection(self):
        cohort = self._get_cohort()
        if cohort:
            cdc = cohort.cohort_genotype_collection
        else:
            cdc = None
        return cdc

    def _get_annotation_kwargs_for_node(self) -> Dict:
        annotation_kwargs = super()._get_annotation_kwargs_for_node()
        if cgc := self.cohort_genotype_collection:
            annotation_kwargs.update(cgc.get_annotation_kwargs())
        return annotation_kwargs

    def _get_cohorts_and_sample_visibility(self):
        # Overrides AnalysisNode
        cohorts = []
        visibility = {}
        if cohort := self._get_cohort():
            cohorts = [cohort]
            visibility = {s: cohort.has_genotype for s in cohort.get_samples()}
        return cohorts, visibility

    @property
    def count_column_prefix(self):
        cohort = self._get_cohort()
        if cohort and cohort.is_sub_cohort():
            column_prefix = "sub_cohort_"
        else:
            cgc = self.cohort_genotype_collection
            column_prefix = f"{cgc.cohortgenotype_alias}__"
        return column_prefix

    @property
    def count_annotation_arg(self):
        """ key in annotation_kwargs """
        return self.cohort_genotype_collection.cohortgenotype_alias

    @property
    def ref_count_column(self):
        return self.count_column_prefix + "ref_count"

    @property
    def hom_count_column(self):
        return self.count_column_prefix + "hom_count"

    @property
    def het_count_column(self):
        return self.count_column_prefix + "het_count"

    @property
    def any_zygosity_count_column(self):
        return self.count_column_prefix + "any_zygosity"

    @property
    def any_germline_count_column(self):
        return self.count_column_prefix + "any_germline"

    def _get_q_and_list(self) -> List[Q]:
        """ Collects node editor filters. Overridden below """
        return self.get_allele_frequency_q_list()

    def get_cohort_and_arg_q_dict(self) -> Tuple[Cohort, Dict[Optional[str], Dict[str, Q]]]:
        arg_q_dict = {}
        cohort = self._get_cohort()
        if cohort:
            cgc = self.cohort_genotype_collection
            q_and = []
            if cohort.is_sub_cohort():
                missing = [Zygosity.UNKNOWN_ZYGOSITY, Zygosity.MISSING]
                sample_zygosities_dict = {s: missing for s in cohort.get_samples()}
                q_sub = cgc.get_zygosity_q(sample_zygosities_dict, exclude=True)
                q_and.append(q_sub)
            q_and.extend(self._get_q_and_list())
            if q_and:
                print(q_and)
                q = reduce(operator.and_, q_and)
                arg_q_dict[cgc.cohortgenotype_alias] = {str(q): q}
        else:
            q_none = self.q_none()
            arg_q_dict[None] = {str(q_none): q_none}
        return cohort, arg_q_dict

    def get_allele_frequency_q_list(self):
        """ Anything that subclasses this (eg TrioNode/PedigreeNode) must also implement
            self.get_samples() and reduce to what is used there so filter is only applied on those samples """
        try:
            naff = self.nodeallelefrequencyfilter
            if not naff.nodeallelefrequencyrange_set.exists():
                return []
        except NodeAlleleFrequencyFilter.DoesNotExist:
            return []

        filters = []
        cgc = self.cohort_genotype_collection

        for sample in self.get_samples():
            # Indexes are handled by cohortgenotype (sub cohorts etc)
            array_index = cgc.get_array_index_for_sample_id(sample.pk)
            # https://docs.djangoproject.com/en/2.1/ref/contrib/postgres/fields/#index-transforms
            allele_frequency_column = f"{cgc.cohortgenotype_alias}__samples_allele_frequency__{array_index}"
            q = naff.get_q(allele_frequency_column, sample.vcf.allele_frequency_percent)
            if q:
                # logging.info("%s => %s => %s", sample_id, allele_frequency_column, q)
                filters.append(q)

        q_and = []
        if filters:
            q_and.append(GroupOperation.reduce(filters, naff.group_operation))
        return q_and

    def get_vcf_locus_filters_arg_q_dict(self) -> Dict[Optional[str], Dict[str, Q]]:
        arg_q_dict = {}
        if self.has_filters:
            vcf = self._get_vcf()
            alias = self.cohort_genotype_collection.cohortgenotype_alias

            node_vcf_filters_qs = NodeVCFFilter.filter_for_node(self, vcf)
            filter_codes = set(node_vcf_filters_qs.values_list("vcf_filter__filter_code", flat=True).distinct())
            if filter_codes:
                q_or = []
                if None in filter_codes:  # Pass
                    filter_codes.remove(None)
                    q_or.append(Q(**{f"{alias}__filters__isnull": True}))

                if filter_codes:
                    joined_codes = re.escape(''.join(filter_codes))
                    joined_codes = joined_codes.replace("'", "''")
                    pattern = f"[{joined_codes}]"
                    q_or.append(Q(**{f"{alias}__filters__regex": pattern}))

                if q_or:
                    q = reduce(operator.or_, q_or)
                    arg_q_dict[alias] = {str(q): q}
        return arg_q_dict

    def get_filter_code(self):
        """
            Used for cached label counts

            0 - No Filters
            1 - Pass Only
            2 - Other (will calculate as not cached) """

        vcf = self._get_vcf()
        if vcf:
            node_vcf_filters_qs = NodeVCFFilter.filter_for_node(self, vcf)
            filter_codes = list(node_vcf_filters_qs.values_list("vcf_filter__filter_code", flat=True).distinct())
            if filter_codes:
                if filter_codes == [None]:  # PASS only
                    return 1
                return 2
        return 0

    def get_filter_description(self):
        FILTER_DESCRIPTIONS = {1: "Pass Filters",
                               2: "Custom Filters"}
        filter_code = self.get_filter_code()
        return FILTER_DESCRIPTIONS.get(filter_code)

    @property
    def has_filters(self):
        vcf = self._get_vcf()
        return vcf and vcf.vcffilter_set.exists()

    def _get_node_extra_columns(self):
        """ show filters if we have them and they're not filtered away (no point then) """

        extra_columns = []
        if self.has_filters:
            cgc = self.cohort_genotype_collection
            extra_columns.append(f"{cgc.cohortgenotype_alias}__filters")

        return extra_columns

    def _get_node_extra_colmodel_overrides(self):
        extra_colmodel_overrides = super()._get_node_extra_colmodel_overrides()
        if self.has_filters:
            vcf = self._get_vcf()
            server_side_formatter = VCFFilter.get_formatter(vcf)
            cgc = self.cohort_genotype_collection
            filters_column = f"{cgc.cohortgenotype_alias}__filters"
            extra_colmodel_overrides[filters_column] = {
                'name': filters_column,
                'model_field': False,  # It's an alias
                'queryset_field': True,
                'server_side_formatter': server_side_formatter,
            }

        return extra_colmodel_overrides


class SampleMixin(CohortMixin):
    """ Adds sample to query via annotation kwargs, must have a "sample" field """

    def _get_annotation_kwargs_for_node(self) -> Dict:
        annotation_kwargs = super()._get_annotation_kwargs_for_node()
        if self.sample:
            annotation_kwargs.update(self.sample.get_annotation_kwargs())
        return annotation_kwargs

    def _get_cohort(self):
        cohort = None
        if self.sample:
            cohort = self.sample.vcf.cohort
        return cohort

    def _get_cohorts_and_sample_visibility_for_node(self):
        cohorts = []
        visibility = {}

        if self.sample:
            cohorts = [self._get_cohort()]
            visibility[self.sample] = self.sample.has_genotype
        return cohorts, visibility


class AncestorSampleMixin(SampleMixin):
    """ Must have a "sample" field that is set from ancestor """

    def _set_sample(self, sample):
        self.sample = sample

    def handle_ancestor_input_samples_changed(self):
        """ Auto-set to single sample ancestor (or remove if no longer ancestor) """

        parent_sample_set = set()
        parents, _errors = self.get_parent_subclasses_and_errors()
        for parent in parents:  # Use parent samples not own as own inserts self.sample
            parent_sample_set.update(parent.get_samples())

        modified = False
        # Don't do anything if new as the get_samples won't work
        if self.version != 0:  # Being set in analysis template
            # may have been moved/copied into a different DAG without current sample as ancestor
            if self.sample and self.sample not in parent_sample_set:
                self._set_sample(None)
                modified = True

        if self.sample is None:
            if proband_sample := self.get_proband_sample():
                self._set_sample(proband_sample)
                modified = True
            elif len(parent_sample_set) == 1:
                self._set_sample(parent_sample_set.pop())
                modified = True

        if modified:
            self.appearance_dirty = True
