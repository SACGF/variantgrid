import operator
import re
from functools import cached_property, reduce
from typing import Optional

import simplejson
from django.db.models import Q

from analysis.models.enums import GroupOperation
from analysis.models.nodes.analysis_node import (
    NodeAlleleFrequencyFilter,
    NodeVCFFilter,
    annotate_and_filter_queryset,
    queryset_to_pk_in_q,
)
from library.genomics.vcf_writer import percent_decode_info_value
from patients.models import Patient
from patients.models_enums import SampleSourceLevel, Zygosity
from patients.sample_grouping import get_patient_for_source
from snpdb.archive import DataArchivedError
from snpdb.models import Cohort, CohortGenotypeCollection, ImportStatus, Sample, VCFFilter, VCFInfo
from snpdb.views.datatable_view import CellData, NullOrder, RichColumn
from upload.models import UploadedVCF
from upload.tso500.dragen_all_fusions_parser import (
    FUSION_OBSERVATIONS_INFO,
    format_fusion_observations,
)


def _render_fusion_calls(cell: CellData) -> str:
    """ The caller rows this fusion was merged from. INFO values are stored as the VCF wrote them -
        htslib doesn't decode, so we do. @see upload.tso500.dragen_all_fusions_parser """
    if not (encoded := cell.value):
        return ""
    return format_fusion_observations(simplejson.loads(percent_decode_info_value(encoded)))


def get_sample_annotation_kwargs(sample: Sample, **kwargs) -> dict:
    """ The genotype join for one sample's VCF, plus its zygosity alias """
    annotation_kwargs = dict(sample.cohort_genotype_collection.get_annotation_kwargs(**kwargs))
    annotation_kwargs.update(sample.get_annotation_kwargs(**kwargs))
    return annotation_kwargs


def get_sample_pk_in_q(node, sample: Sample, arg_q_dict: dict[Optional[str], dict[str, Q]]) -> Q:
    """ One sample's variants as pk IN (subquery). Annotated with only its own VCF's genotype join,
        so the subquery doesn't drag the other samples' outer joins through with it """
    qs = node._get_model_queryset()  # pylint: disable=protected-access
    a_kwargs = get_sample_annotation_kwargs(sample)
    qs, q_list = annotate_and_filter_queryset(qs, a_kwargs, arg_q_dict)
    if q_list:
        qs = qs.filter(reduce(operator.and_, q_list))
    return queryset_to_pk_in_q(qs)


def get_sample_any_zygosity_arg_q_dict(sample: Sample) -> dict[Optional[str], dict[str, Q]]:
    """ A sample's rows, unfiltered - the zygosity IN is what restricts the outer join to them, so a
        VCF with nothing to filter on passes through rather than being left out
        (@see analysis/models/nodes/sources/sample_node.py:SampleNode._get_sample_arg_q_dict) """
    alias, field = sample.get_cohort_genotype_alias_and_field("zygosity")
    q = Q(**{f"{field}__in": [code for code, _ in Zygosity.CHOICES]})
    return {alias: {str(q): q}}


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
        """ Use cohort genotype in the key, as that can change if a VCF is reloaded.
            Also include the sub-cohort any-sample-called VariantCollection pk so the cache invalidates
            when a VC build completes (or is dropped) under the same CGC. @see issue #1551 """
        cache_key = super()._get_cache_key()
        cgc_id = 0
        vc_id = 0
        if cgc := self.cohort_genotype_collection:
            cgc_id = cgc.pk
        cohort = self._get_cohort()
        if cohort and cohort.is_sub_cohort:
            if vc := cohort.get_any_sample_called_variant_collection():
                vc_id = vc.pk
        return "_".join((cache_key, str(cgc_id), str(vc_id)))

    @property
    def cohort_genotype_collection(self):
        cohort = self._get_cohort()
        if cohort:
            try:
                cdc = cohort.cohort_genotype_collection
            except (DataArchivedError, CohortGenotypeCollection.DoesNotExist):
                # Surface via _get_configuration_errors → ERROR_CONFIGURATION;
                # keep node-internal "no source" code paths working.
                # A missing CGC (mid-reload / version mismatch) is treated the same
                # as archived data so the node shows a config error instead of 500ing.
                cdc = None
        else:
            cdc = None
        return cdc

    def _get_annotation_kwargs_for_node(self, **kwargs) -> dict:
        annotation_kwargs = super()._get_annotation_kwargs_for_node(**kwargs)
        if cgc := self.cohort_genotype_collection:
            annotation_kwargs.update(cgc.get_annotation_kwargs(**kwargs))
        cohort = self._get_cohort()
        if cohort and cohort.is_sub_cohort:
            # Register the pre-computed any-sample-called VC alias so get_cohort_and_arg_q_dict can join
            # to it instead of running the EXCLUDE regex. @see issue #1551
            if vc := cohort.get_any_sample_called_variant_collection():
                annotation_kwargs.update(vc.get_annotation_kwargs(**kwargs))
        return annotation_kwargs

    def _get_cohorts_and_sample_visibility(self):
        # Overrides AnalysisNode
        cohorts = []
        visibility = {}
        if cohort := self._get_cohort():
            cohorts = [cohort]
            visibility = dict.fromkeys(cohort.get_samples(), cohort.has_sample_columns)
        return cohorts, visibility

    @property
    def count_column_prefix(self):
        cohort = self._get_cohort()
        if cohort and cohort.is_sub_cohort:
            return f"sub_cohort_{cohort.pk}_"
        if cgc := self.cohort_genotype_collection:
            return f"{cgc.cohortgenotype_alias}__"
        return None

    @property
    def non_ref_call_count_annotation_arg(self):
        """
            The ..._annotation_args are used to be able to put queries into the right q_arg_dict key
            You need to be able to apply the Q after the right annotate() - in AllVariantsNode and Cohort
            we just need to add to the same alias (that joins to the table) and then use the SQL columns
            In sub cohorts, we need to build fake ones via annotate and then add Qs there

            See AbstractZygosityCountNode.get_zygosity_count_arg_q_dict and
            CohortNode._get_annotation_kwargs_for_node
        """
        return self.non_ref_call_count_column  # This is always an annotation

    @property
    def ref_count_annotation_arg(self):
        """ key in annotation_kwargs """
        cohort = self._get_cohort()
        if cohort and cohort.is_sub_cohort:
            return self.ref_count_column
        return self.cohort_genotype_collection.cohortgenotype_alias

    @property
    def het_count_annotation_arg(self):
        """ key in annotation_kwargs """
        cohort = self._get_cohort()
        if cohort and cohort.is_sub_cohort:
            return self.het_count_column
        return self.cohort_genotype_collection.cohortgenotype_alias

    @property
    def hom_count_annotation_arg(self):
        """ key in annotation_kwargs """
        cohort = self._get_cohort()
        if cohort and cohort.is_sub_cohort:
            return self.hom_count_column
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
    def any_call_count_column(self):
        return self.count_column_prefix + "any_call"

    @property
    def non_ref_call_count_column(self):
        return self.count_column_prefix + "non_ref"

    def _get_q_and_list(self) -> list[Q]:
        """ Collects node editor filters. Overridden below """
        return self.get_allele_frequency_q_list()

    def get_cohort_and_arg_q_dict(self) -> tuple[Cohort, dict[Optional[str], dict[str, Q]]]:
        arg_q_dict = {}
        cohort = self._get_cohort()
        if cohort:
            cgc = self.cohort_genotype_collection
            q_and = []
            if cohort.is_sub_cohort:
                if vc := cohort.get_any_sample_called_variant_collection():
                    # Pre-computed any-sample-called set => hash join instead of regex seq-scan (#1551).
                    # Keyed under the VC alias (registered in _get_annotation_kwargs_for_node) so it's
                    # filtered after its own annotate, not the cohortgenotype one.
                    q_vc = Q(**{f"{vc.variant_collection_alias}__isnull": False})
                    arg_q_dict[vc.variant_collection_alias] = {str(q_vc): q_vc}
                else:
                    missing = [Zygosity.UNKNOWN_ZYGOSITY, Zygosity.MISSING]
                    sample_zygosities_dict = dict.fromkeys(cohort.get_samples(), missing)
                    q_sub = cgc.get_zygosity_q(sample_zygosities_dict, exclude=True)
                    q_and.append(q_sub)
            q_and.extend(self._get_q_and_list())
            if q_and:
                q = reduce(operator.and_, q_and)
                arg_q_dict[cgc.cohortgenotype_alias] = {str(q): q}
        else:
            q_none = self.q_none()
            arg_q_dict[None] = {str(q_none): q_none}
        return cohort, arg_q_dict

    def get_allele_frequency_q_list(self):
        """ Anything that subclasses this (eg TrioNode/PedigreeNode) must also implement
            self.get_samples_with_genotype() and reduce to what is used there so filter is only
            applied on those samples """
        try:
            naff = self.nodeallelefrequencyfilter
            if not naff.nodeallelefrequencyrange_set.exists():
                return []
        except NodeAlleleFrequencyFilter.DoesNotExist:
            return []

        filters = []
        cgc = self.cohort_genotype_collection
        packed_index_by_sample_id = cgc.get_packed_index_by_sample_id

        for sample in self.get_samples_with_genotype():
            # get_samples_with_genotype() includes ancestor samples (eg a compound het Trio/Quad's parent node)
            # that may not be in this cohort's genotype array - skip those.
            if sample.pk not in packed_index_by_sample_id:
                continue
            # Indexes are handled by cohortgenotype (sub cohorts etc)
            array_index = packed_index_by_sample_id[sample.pk]
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

    def _get_vcf_locus_filters_arg_q_dict(self, vcf, alias: str,
                                          pass_only: Optional[bool] = None) -> dict[Optional[str], dict[str, Q]]:
        """ Filter ids are stored on the node - they resolve into each VCF's own codes here.
            pass_only lets a caller decide PASS for itself (see SampleNode's per sample overrides) """
        arg_q_dict = {}
        filter_codes = NodeVCFFilter.get_filter_codes(self, vcf, pass_only=pass_only)
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

    def get_vcf_locus_filter_vcfs(self) -> list:
        """ The VCFs the node's filter id selection is offered from / resolved against """
        if vcf := self._get_vcf():
            return [vcf]
        return []

    def get_vcf_locus_filters_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        if not self.has_filters:
            return {}
        alias = self.cohort_genotype_collection.cohortgenotype_alias
        return self._get_vcf_locus_filters_arg_q_dict(self._get_vcf(), alias)

    def get_filter_code(self):
        """
            Used for cached label counts

            0 - No Filters
            1 - Pass Only
            2 - Other (will calculate as not cached) """

        filter_codes = set()
        for vcf in self.get_vcf_locus_filter_vcfs():
            filter_codes |= NodeVCFFilter.get_filter_codes(self, vcf)
        if not filter_codes:
            return 0
        if filter_codes == {None}:  # PASS only, which means the same thing in every VCF
            return 1
        return 2

    def get_filter_description(self):
        FILTER_DESCRIPTIONS = {1: "Pass Filters",
                               2: "Custom Filters"}
        filter_code = self.get_filter_code()
        return FILTER_DESCRIPTIONS.get(filter_code)

    @cached_property
    def has_filters(self) -> bool:
        """ Cached: the query build asks this once per sample, and a VCF's filters can't change
            mid request """
        vcfs = self.get_vcf_locus_filter_vcfs()
        return bool(vcfs) and VCFFilter.objects.filter(vcf__in=vcfs).exists()

    def _get_filters_cohort_genotype_collections(self) -> list:
        """ The genotype collections whose record level FILTER to show. Nodes spanning VCFs override """
        if self.has_filters:
            if cgc := self.cohort_genotype_collection:
                return [cgc]
        return []

    def _get_node_extra_columns(self) -> list[RichColumn]:
        """ show filters if we have them and they're not filtered away (no point then) """
        extra_columns = super()._get_node_extra_columns()
        cgcs = self._get_filters_cohort_genotype_collections()
        for cgc in cgcs:
            vcf = cgc.cohort.get_vcf()
            filters_column = f"{cgc.cohortgenotype_alias}__filters"
            # Which VCF's FILTER this is only needs saying when there are several
            label = f"{vcf} Filters" if len(cgcs) > 1 else "Filters"
            extra_columns.append(RichColumn(
                key=filters_column, label=label, width=80, orderable=True, search=False,
                include_in_csv=True,
                # Expanded to the VCF's own filter descriptions server side, so the CSV matches
                renderer=VCFFilter.get_formatter(vcf), csv_rendered=True,
                # Nearly every record passed - the cell fades a pass right down so only a
                # failure reads. @see VariantGridFormat.vcfFilters
                client_renderer='VariantGridFormat.vcfFilters',
                null_order=NullOrder.FIRST_ON_ASC))

        # One gene pair can be several caller rows, all merged onto the one Variant - blank on
        # everything that isn't a fusion, which is why the column only appears for a fusion VCF
        fusion_cgcs = self._get_fusion_calls_cohort_genotype_collections()
        for cgc in fusion_cgcs:
            label = f"{cgc.cohort.get_vcf()} Fusion calls" if len(fusion_cgcs) > 1 else "Fusion calls"
            extra_columns.append(RichColumn(
                key=f"{cgc.cohortgenotype_alias}__info__{FUSION_OBSERVATIONS_INFO}",
                label=label, width=90,
                orderable=False, search=False, include_in_csv=True,
                renderer=_render_fusion_calls, csv_rendered=True,
                client_renderer='VariantGridFormat.fusionCalls'))
        return extra_columns

    def _get_fusion_calls_cohort_genotype_collections(self) -> list:
        """ The genotype collections whose VCF carries fusion observations. Nodes spanning VCFs override """
        if cgc := self.cohort_genotype_collection:
            if VCFInfo.objects.filter(vcf=cgc.cohort.get_vcf(), identifier=FUSION_OBSERVATIONS_INFO).exists():
                return [cgc]
        return []

    def _get_configuration_check_cohorts(self) -> list:
        """ Cohorts to check for missing/archived genotype data. Nodes spanning VCFs override """
        if cohort := self._get_cohort():
            return [cohort]
        return []

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        for cohort in self._get_configuration_check_cohorts():
            if cohort.import_status != ImportStatus.SUCCESS:
                errors.append(f"'{cohort}' has import status: {cohort.get_import_status_display()}")
                continue

            try:
                _ = cohort.cohort_genotype_collection
            except CohortGenotypeCollection.DoesNotExist:
                errors.append("Source data missing: underlying genotype data is no longer available")
            except DataArchivedError as e:
                errors.append(str(e))

            if vcf := cohort.get_vcf():
                try:
                    uv: UploadedVCF = vcf.uploadedvcf
                    if uv.max_variant_id:  # Very old VCFs may not have this set
                        variant_annotation_version = self.analysis.annotation_version.variant_annotation_version
                        if not uv.is_fully_annotated(variant_annotation_version):
                            errors.append(f"VCF '{vcf}' contains variants that have not finished annotation"
                                          f" (in variant annotation version={variant_annotation_version})")
                except UploadedVCF.DoesNotExist:
                    pass
        return errors


class SampleMixin(CohortMixin):
    """ Adds sample to query via annotation kwargs, must have a "sample" field """

    def _get_sample(self) -> Optional[Sample]:
        """ The sample this node's genotype joins hang off - overridden by nodes that group samples """
        return self.sample

    def _get_annotation_kwargs_for_node(self, **kwargs) -> dict:
        kwargs["override"] = False
        annotation_kwargs = super()._get_annotation_kwargs_for_node(**kwargs)
        if sample := self._get_sample():
            annotation_kwargs.update(sample.get_annotation_kwargs(**kwargs))
        return annotation_kwargs

    def _get_cohort(self):
        cohort = None
        if sample := self._get_sample():
            cohort = sample.vcf.cohort
        return cohort

    def _get_cohorts_and_sample_visibility_for_node(self):
        cohorts = []
        visibility = {}

        if sample := self._get_sample():
            cohorts = [self._get_cohort()]
            visibility[sample] = sample.has_sample_columns
        return cohorts, visibility


class AncestorSampleMixin(SampleMixin):
    """ A filter node that applies to either one sample or one patient, set from its ancestors.

        The model needs a "sample" and a "patient" field, at most one of which is set - both null
        means unset. In patient mode the filter applies to every ancestor sample of that patient and
        the node's query is the OR of the per-sample filters, the shape a group level SampleNode
        produces (@see analysis/models/nodes/sources/sample_node.py:SampleNode). "Every ancestor
        sample" is the scope rather than every sample of the patient, so the source node decides the
        reach and the filter follows it. """

    def _set_sample(self, sample):
        self.sample = sample
        self.patient = None

    def _set_patient(self, patient):
        self.patient = patient
        self.sample = None

    def get_filter_patient(self) -> Optional[Patient]:
        """ Who the node is about - the patient it was set to, or the one its sample belongs to """
        if self.patient:
            return self.patient
        return get_patient_for_source(SampleSourceLevel.SAMPLE, self.sample)

    def get_filter_samples(self) -> list[Sample]:
        """ The samples this node's genotype filters apply to - every path that used self.sample
            for a genotype join goes through here """
        if self.sample:
            return [self.sample]
        if self.patient:
            samples = [s for s in self.get_ancestor_samples()
                       if get_patient_for_source(SampleSourceLevel.SAMPLE, s) == self.patient]
            return sorted(samples, key=lambda s: s.pk)
        return []

    def _get_filter_samples_arg_q_dict(self, per_sample) -> dict[Optional[str], dict[str, Q]]:
        """ The one place the per-sample OR is built. per_sample(sample) returns that sample's
            arg_q_dict, keyed on its own alias.

            Sample mode returns it as is. Patient mode wraps each sample's filters in a pk__in
            subquery annotated with only that sample's genotype join - a Q keyed on an alias runs as
            soon as that alias is annotated (@see analysis_node.annotate_and_filter_queryset), so an
            OR across two aliases has nowhere to hang until both are """
        samples = self.get_filter_samples()
        if not samples:
            return {}
        if self.sample:
            return per_sample(samples[0])
        q = reduce(operator.or_, [get_sample_pk_in_q(self, sample, per_sample(sample)) for sample in samples])
        return {None: {self._get_node_q_hash(): q}}

    def _get_cohorts_and_sample_visibility_for_node(self):
        if not self.patient:
            return super()._get_cohorts_and_sample_visibility_for_node()

        cohorts = []
        visibility = {}
        for sample in self.get_filter_samples():
            cohort = sample.vcf.cohort
            if cohort not in cohorts:
                cohorts.append(cohort)
            visibility[sample] = sample.has_sample_columns
        return cohorts, visibility

    def _get_annotation_kwargs_for_node(self, **kwargs) -> dict:
        annotation_kwargs = super()._get_annotation_kwargs_for_node(**kwargs)
        if self.patient:
            kwargs["override"] = False
            for sample in self.get_filter_samples():
                annotation_kwargs.update(get_sample_annotation_kwargs(sample, **kwargs))
        return annotation_kwargs

    def _get_configuration_errors(self) -> list:
        errors = super()._get_configuration_errors()
        if self.sample:
            if self.sample not in self.get_ancestor_samples():
                errors.append(f"Sample: {self.sample} is not set as a sample in any ancestors of this node")
        elif self.patient:
            if not self.get_filter_samples():
                errors.append(f"Patient: {self.patient} has no samples in any ancestors of this node")
        return errors

    def get_ancestor_samples(self) -> set[Sample]:
        """ Get all samples from ancestor nodes, including those from VCFs without genotypes,
            so that variant-only VCFs (has_sample_columns=False) are still valid ancestors """
        parent_sample_set = set()
        parents, _errors = self.get_parent_subclasses_and_errors()
        for parent in parents:  # Use parent samples not own as own inserts self.sample
            parent_sample_set.update(parent.get_samples())
        return parent_sample_set

    def handle_ancestor_input_samples_changed(self):
        """ Auto-set to the ancestors' proband (or remove if no longer reachable from them) """

        parent_sample_set = self.get_ancestor_samples()

        modified = False
        # Don't do anything if new as the get_samples won't work
        if self.version != 0:  # Being set in analysis template
            # may have been moved/copied into a different DAG without current sample as ancestor
            if self.sample and self.sample not in parent_sample_set:
                self._set_sample(None)
                modified = True
            if self.patient and not self.get_filter_samples():
                self._set_patient(None)
                modified = True

        if self.sample is None and self.patient is None:
            if proband_sample := self.get_proband_sample():
                self._set_sample(proband_sample)
                modified = True
            elif proband_patient := self.get_proband_patient():
                # Several callers on the one extraction have no single sample, but are one person
                self._set_patient(proband_patient)
                modified = True
            elif len(parent_sample_set) == 1:
                self._set_sample(parent_sample_set.pop())
                modified = True

        if modified:
            self.appearance_dirty = True
