""" AnalysisNode is the base class that all analysis nodes inherit from. """
import abc
import logging
import operator
from collections import defaultdict
from collections.abc import Sequence
from datetime import timedelta
from functools import cached_property, reduce
from random import random
from time import time
from typing import Optional

from auditlog.context import disable_auditlog
from auditlog.models import AuditlogHistoryField, LogEntry
from auditlog.registry import auditlog
from cache_memoize import cache_memoize
from celery.canvas import Signature
from django.conf import settings
from django.contrib.postgres.fields import ArrayField
from django.core.cache import cache
from django.core.exceptions import FieldError
from django.db import connection, models, transaction
from django.db.models import IntegerField, QuerySet, Value
from django.db.models.aggregates import Count
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.expressions import RawSQL
from django.db.models.query_utils import Q
from django.db.models.signals import post_delete
from django.dispatch import receiver
from django.utils import timezone
from django_dag.models import edge_factory, node_factory
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager

from analysis.exceptions import (
    NodeConfigurationException,
    NodeNotFoundException,
    NodeOutOfDateException,
    NodeParentErrorsException,
    NodeParentNotReadyException,
    NonFatalNodeError,
)
from analysis.models.enums import (
    AnalysisTemplateType,
    GroupOperation,
    NodeColors,
    NodeErrorSource,
    NodeStatus,
)
from analysis.models.models_analysis import Analysis
from analysis.models.nodes.node_counts import get_node_counts_and_labels_dict, get_node_extra_filters_q
from analysis.models.nodes.node_display import NodeChip, NodeIcon, NodeMenuEntry
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from classification.models import Classification
from library.constants import DAY_SECS, MINUTE_SECS
from library.django_utils import thread_safe_unique_together_get_or_create
from library.django_utils.django_postgres import get_backend_pid
from library.log_utils import log_traceback
from library.utils import add_exception_note, format_percent
from library.utils.database_utils import queryset_to_sql
from library.utils.django_utils import get_model_content_type_dict
from snpdb.models import (
    AlleleSource,
    BuiltInFilters,
    Cohort,
    Contig,
    GenomeBuild,
    ProcessingStatus,
    Sample,
    SampleFilePath,
    Variant,
    VariantCollection,
    VCFFilter,
    Wiki,
)
from snpdb.variant_collection import write_sql_to_variant_collection
from snpdb.views.datatable_view import RichColumn

# How long a node's lease is good for. The window is (re)started when a worker claims the node for
# loading, so it measures actual load time rather than how long the task sat in the queue.
LEASE_SECONDS = MINUTE_SECS * 10


def queryset_to_pk_in_q(qs: QuerySet) -> Q:
    """ Embed a queryset as pk IN (subquery), NOT list(qs.values_list("pk")).
        Callers reach this with querysets that are large by construction - small ones are substituted
        to a literal PK list upstream by get_small_parent_arg_q_dict - so list() here could pull an
        unbounded number of PKs into Python (e.g. a 7.4M-row cohort).
        We render to RawSQL, capturing the compiled SQL + params (incl. the partition table rewrite the
        TransformerQuerySet applies in as_sql), so it runs as a single DB-side semi-join. RawSQL also keeps
        arg_q_dict picklable for the q-cache cache.set: a live TransformerQuerySet embedded in the Q is
        unpicklable (closure-local compiler classes) which previously forced the materialise-to-list
        workaround (issue #546, #240, ad35a7fb1). """
    pk_qs = qs.values_list("pk", flat=True)
    sql, params = pk_qs.query.sql_with_params()
    return Q(pk__in=RawSQL(sql, params))


def annotate_and_filter_queryset(qs: QuerySet, a_kwargs: dict, arg_q_dict: dict) -> tuple[QuerySet, list[Q]]:
    """ Apply each annotation then the filters that use it, so that it forces an inner query - applying
        them all at once can join to the same table twice. Returns the queryset plus the filters that
        don't rely on an annotation (arg=None), for the caller to combine with its own.

        arg_q_dict is consumed - anything left over means a filter had no annotation to hang off. """
    for k, v in a_kwargs.items():
        qs = qs.annotate(**{k: v})
        if q_and_list := list(arg_q_dict.pop(k, {}).values()):
            q = reduce(operator.and_, q_and_list)
            try:
                qs = qs.filter(q)
            except FieldError as fe:
                add_exception_note(fe, f"annotation kwarg: {k}: {q=}.")
                raise

    q_list = []
    # Anything stored under None means filters that don't rely on annotation - do afterwards
    if q_dict := arg_q_dict.pop(None, {}):
        q_list.extend(q_dict.values())

    if arg_q_dict:
        raise ValueError(f"arg_q_dict filters {arg_q_dict.keys()} not applied (missing in {a_kwargs=})")
    return qs, q_list


def _default_position():
    return 10 + random() * 50


class NodeInheritanceManager(InheritanceManager):
    def get_queryset(self):
        queryset = super()._queryset_class(self.model)
        return queryset.select_related("analysis",
                                       "analysis__user",
                                       "analysis__genome_build",
                                       "analysis__annotation_version",
                                       "analysis__annotation_version__variant_annotation_version",
                                       "analysis__annotation_version__variant_annotation_version__gene_annotation_release")


class NodeAuditLogMixin:
    @abc.abstractmethod
    def _get_node(self):
        pass

    def get_additional_data(self):
        """ For django-audit-log """
        node = self._get_node()
        return {
            "analysis_id": node.analysis_id,
            "node_id": node.pk,
        }


class AnalysisNode(NodeAuditLogMixin, node_factory('AnalysisEdge', base_model=TimeStampedModel)):
    model = Variant
    objects = NodeInheritanceManager()
    history = AuditlogHistoryField()
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE)
    name = models.TextField(blank=True)
    x = models.IntegerField(default=_default_position)
    y = models.IntegerField(default=_default_position)
    version = models.IntegerField(default=0)  # Queryset version
    appearance_version = models.IntegerField(default=0)
    auto_node_name = models.BooleanField(default=True)
    output_node = models.BooleanField(default=False)
    hide_node_and_descendants_upon_template_configuration_error = models.BooleanField(default=False)
    # Waives the errors on get_ignorable_error_fields() - they're reported as warnings instead
    ignore_field_errors = models.BooleanField(default=False)
    ready = models.BooleanField(default=True)
    valid = models.BooleanField(default=False)
    visible = models.BooleanField(default=True)
    count = models.IntegerField(null=True, default=None)
    errors = models.TextField(null=True)
    shadow_color = models.TextField(null=True)
    load_seconds = models.FloatField(null=True)
    # This is set to node/version you cloned - cleared upon modification
    cloned_from = models.ForeignKey('NodeVersion', null=True, on_delete=SET_NULL)
    status = models.CharField(max_length=1, choices=NodeStatus.choices, default=NodeStatus.DIRTY)

    PARENT_CAP_NOT_SET = -1
    min_inputs = 1
    max_inputs = 1
    uses_parent_queryset = True
    disabled = False

    UPDATE_TASK = "analysis.tasks.node_update_tasks.update_node_task"
    NODE_CACHE_TASK = "analysis.tasks.node_update_tasks.node_cache_task"
    WAIT_FOR_CACHE_TASK = "analysis.tasks.node_update_tasks.wait_for_cache_task"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.appearance_dirty = False
        self.ancestor_input_samples_changed = False
        self.parents_changed = False
        self.queryset_dirty = False
        self.update_children = True
        self._cached_parents = None
        self._cached_analysis_errors = None
        self._cache_node_q = settings.ANALYSIS_NODE_CACHE_Q  # Disable for unit tests

    def __lt__(self, other):
        return self.pk < other.pk

    @cached_property
    def num_samples_for_build(self) -> int:
        return Sample.objects.filter(vcf__genome_build=self.analysis.genome_build).count()

    def get_subclass(self):
        """ Returns the node loaded as a subclass """
        return AnalysisNode.objects.get_subclass(pk=self.pk)

    def log_entry_qs(self) -> QuerySet[LogEntry]:
        # This is put on there via AnalysisNode.get_additional_data
        return LogEntry.objects.filter(additional_data__analysis_id=self.analysis_id,
                                       additional_data__node_id=self.pk)

    def has_audit_log(self) -> bool:
        return self.log_entry_qs().exists()

    def _get_node(self):
        return self

    def check_still_valid(self):
        """ Checks that the node is still there and has the version we expect - or throw exception """
        version_qs = AnalysisNode.objects.filter(pk=self.pk).values_list("version", flat=True)
        if version_qs:
            db_version = version_qs[0]
            if db_version > self.version:
                raise NodeOutOfDateException()
        else:
            raise NodeNotFoundException(self.pk)

    def _get_cohorts_and_sample_visibility_for_node(self) -> tuple[Sequence[Cohort], dict]:
        """ Visibility = can see on grid """
        return [], {}

    @staticmethod
    def _get_visible_samples_from_cohort(cohorts, visibility):
        samples = set()
        for c in cohorts:
            for s in c.get_samples():
                if visibility.get(s):
                    samples.add(s)
        return sorted(samples)

    def _get_model_queryset(self):
        self.analysis.check_valid()
        return get_variant_queryset_for_annotation_version(self.analysis.annotation_version)

    def get_cohorts_and_sample_visibility(self, sort=True) -> tuple[Sequence[Cohort], dict]:
        """ Returns all node + ancestor cohorts (and visibilities of their samples)

            The underlying data for all samples/cohorts/sub-cohorts/trios/pedigrees is Cohorts, so need to know which
            to retrieve from DB (and what sample info to extract from packed columns) to filter + show on grid """

        cohorts, visibility = self._get_cohorts_and_sample_visibility_for_node()
        cohorts = set(cohorts)
        if self.has_input():
            parents, _ = self.get_parent_subclasses_and_errors()
            for parent in parents:
                c, v = parent.get_cohorts_and_sample_visibility(sort=False)
                cohorts.update(c)
                visibility.update(v)

        # May have sub-cohorts, so get unique base cohorts
        cohorts = {c.get_base_cohort() for c in cohorts}
        if sort:
            cohorts = sorted(cohorts)
        return cohorts, visibility

    @cache_memoize(DAY_SECS, args_rewrite=lambda s: (s.pk, s.version))
    def get_sample_ids_with_genotype(self) -> list[int]:
        return [s.pk for s in self.get_samples_with_genotype()]

    def get_samples_from_node_only_not_ancestors(self):
        _, visibility = self._get_cohorts_and_sample_visibility_for_node()
        return sorted(visibility)

    def _get_proband_sample_for_node(self) -> Optional[Sample]:
        """ Sample of the object of a study, if known """
        return None

    def get_proband_sample(self) -> Optional[Sample]:
        """ Sample of the object of a study if known """
        proband_samples = set()
        if proband_sample := self._get_proband_sample_for_node():
            proband_samples.add(proband_sample)

        if self.has_input():
            parents, _ = self.get_parent_subclasses_and_errors()
            for parent in parents:
                if parent_proband_sample := parent.get_proband_sample():
                    proband_samples.add(parent_proband_sample)

        proband_sample = None
        if len(proband_samples) == 1:  # If ambiguous, then just give up
            proband_sample = proband_samples.pop()
        return proband_sample

    def get_samples_with_genotype(self) -> list[Sample]:
        """ Node + ancestor samples whose genotype we can show/filter on - ie variant-only VCFs
            (has_genotype=False) are left out. Use get_samples() for sample level data """
        cohorts, visibility = self.get_cohorts_and_sample_visibility(sort=False)
        return self._get_visible_samples_from_cohort(cohorts, visibility)

    def get_samples(self) -> list[Sample]:
        """ Every node + ancestor sample, including those from variant-only VCFs. Sample level data
            that doesn't come from the genotype - gene lists, coverage, BAMs, patients - and
            restricting sample fields to ancestors @see AncestorSampleMixin """
        _, visibility = self.get_cohorts_and_sample_visibility(sort=False)
        return sorted(visibility)  # Every sample the node knows about is a key, genotype or not

    def get_sample_ids(self) -> list[int]:
        return [s.pk for s in self.get_samples()]

    def get_bams_dict(self):
        bams_dict = defaultdict(set)
        sfp_qs = SampleFilePath.objects.filter(sample__in=self.get_samples())
        for sample_id, file_path in sfp_qs.values_list("sample_id", "file_path"):
            bams_dict[sample_id].add(file_path)
        return {k: list(v) for k, v in bams_dict.items()}

    def get_connection_data(self, parent):
        """ Return dict of source_id/target_id for sending as JSON """
        return {"source_id": parent.get_css_id(),
                "target_id": self.get_css_id()}

    def get_rendering_args(self):
        return {}

    def get_css_id(self):
        if self.pk:
            css_id = f"analysis-node-{self.pk}"
        else:
            css_id = None

        return css_id

    def get_update_task(self):
        return Signature(self.UPDATE_TASK, args=(self.pk, self.version), immutable=True)

    def get_cache_task_args_set(self, force_cache=False):
        """ returns Celery tasks which are called in node_utils.get_analysis_update_task before children are loaded
            Uses tasks not signatures so they are hashable in a set to be able to remove dupes """

        task_args_set = set()
        if self.is_valid and (force_cache or self.use_cache):
            if parent := self.get_unmodified_single_parent_node():
                return parent.get_cache_task_args_set(force_cache=force_cache)

            node_cache, created = NodeCache.get_or_create_for_node(self)
            if created:
                task_args_set.add((self.NODE_CACHE_TASK, (self.pk, self.version)))
            # else: cache is already being built by another node - do NOT add WAIT_FOR_CACHE_TASK.
            #       The scheduler's dependency gate (_node_ready_to_lease) holds this node back while
            #       that cache is actively building (PROCESSING); node_cache_task re-triggers the
            #       dispatcher on completion (issue #346).
        return task_args_set

    def get_parent_subclasses_and_errors(self):
        if self._cached_parents is None:
            qs = AnalysisNode.objects.filter(children=self.id, children__isnull=False)
            self._cached_parents = list(qs.select_subclasses())

        parents = self._cached_parents

        num_parents = len(parents)
        errors = []
        if self.min_inputs != AnalysisNode.PARENT_CAP_NOT_SET and num_parents < self.min_inputs:
            errors.append((NodeErrorSource.CONFIGURATION, f"{num_parents} parents < minimum of {self.min_inputs}"))
        elif self.max_inputs != AnalysisNode.PARENT_CAP_NOT_SET and num_parents > self.max_inputs:
            errors.append((NodeErrorSource.CONFIGURATION, f"{num_parents} parents > maximum of {self.max_inputs}"))

        for parent in parents:
            if NodeStatus.is_error(parent.status):
                errors.append((NodeErrorSource.PARENT, "Parent has errors"))
                break

        return parents, errors

    def get_parent_subclasses(self):
        """ Gets parents, throws an Exception if any errors """
        parents, errors = self.get_parent_subclasses_and_errors()
        if errors:
            AnalysisNode.throw_errors_exception(errors)
        return parents

    def get_non_empty_parents(self, require_parents_ready=True) -> list['AnalysisNode']:
        """ Returns non-empty (count > 0) parents.
            If require_parents_ready=True, die if parents not ready
            Otherwise, return them as we don't know if they're empty or not """
        non_empty_parents = []
        for p in self.get_parent_subclasses():
            if p.is_ready:
                if p.count == 0:
                    continue
            elif require_parents_ready:
                raise NodeParentNotReadyException(f"Parent {p} is not ready!")
            non_empty_parents.append(p)
        return non_empty_parents

    def _get_live_data_sources(self) -> dict[str, int]:
        """ Mutable tables this node's own query reads, keyed by a stable source key, with the source's
            version at the time of the call. Empty means the node's query is a pure function of versioned
            inputs (annotation version, node settings, parents) """
        return {}

    def get_live_data_sources(self) -> dict[str, int]:
        sources = dict(self._get_live_data_sources())
        for parent in self.get_non_empty_parents():
            sources.update(parent.get_live_data_sources())
        return sources

    @property
    def live_data_sources(self) -> dict[str, int]:
        """ The live sources recorded when this node version was loaded """
        return self.node_version.live_data_sources

    @property
    def count_is_deterministic(self) -> bool:
        """ False means the count is advisory - the node reads tables that change under it """
        return not self.live_data_sources

    def get_live_data_notes(self) -> list[str]:
        """ Editor notes explaining a node's live data sources - unlike get_warnings() these aren't faults """
        return []

    def _raise_or_warn_count_mismatch(self, detail: str):
        """ A deterministic node re-running the same query must give the same answer, so a mismatch is a
            query bug. A live-source node's data moves under it, so it's expected drift """
        msg = f"Node {self}(pk={self.pk}) count mismatch: {detail}"
        if self.count_is_deterministic:
            raise ValueError(msg)
        logging.warning("%s (live sources: %s)", msg, self.live_data_sources)

    def get_single_parent(self):
        if self.min_inputs != 1:
            msg = "get_single_parent() should only be called for single parent nodes"
            raise ValueError(msg)

        parents, errors = self.get_parent_subclasses_and_errors()
        if errors:
            errors = AnalysisNode.flatten_errors(errors)
            msg = "Parent had errors: " + ', '.join(errors)
            raise NonFatalNodeError(msg)
        num_parents = len(parents)
        if num_parents != 1:
            msg = f"get_single_parent() called for node with {num_parents} parents"
            raise ValueError(msg)
        return parents[0]

    def get_single_parent_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = {}
        parent = self.get_single_parent()
        if parent.is_ready:
            if parent.count == 0:
                q_none = self.q_none()
                arg_q_dict[None] = {str(q_none): q_none}
            elif (small_arg_q_dict := AnalysisNode.get_small_parent_arg_q_dict(parent)) is not None:
                arg_q_dict = small_arg_q_dict
            else:
                arg_q_dict = parent.get_arg_q_dict()
        else:
            # This should never happen...
            raise ValueError("get_single_parent_q called when single parent not ready!!!")
        return arg_q_dict

    def get_single_parent_contigs(self):
        parent = self.get_single_parent()
        if parent.is_ready:
            if parent.count == 0:
                contigs = set()
            else:
                contigs = parent.get_contigs()
        else:
            # This should never happen...
            raise ValueError("get_single_parent_contigs called when single parent not ready!!!")
        return contigs

    def _get_kwargs_for_parent_annotation_kwargs(self, **kwargs) -> dict:
        """ Use this to pass messages up through to parents """
        return {}

    def _get_annotation_kwargs_for_node(self, **kwargs) -> dict:
        """ Override this method per-node.
            Any key/values in here MUST be consistent - as annotation_kwargs from multiple
            nodes may be combined in the MergeNode
        """
        annotation_kwargs = {}
        if self.node_cache:
            annotation_kwargs.update(self.node_cache.variant_collection.get_annotation_kwargs(**kwargs))
        return annotation_kwargs

    def _has_common_variants(self) -> bool:
        if self.has_input():
            return any(parent._has_common_variants() for parent in self.get_non_empty_parents())
        return True

    def get_annotation_kwargs(self, **kwargs) -> dict:
        """ Passed to Variant QuerySet annotate()
            Can be used w/FilteredRelation to force a join to a partition, in which case you need to use
            the alias given in annotate. @see https://github.com/SACGF/variantgrid/wiki/Data-Partitioning """

        kwargs.update(self._get_kwargs_for_parent_annotation_kwargs(**kwargs))
        # Only apply parent annotation kwargs if you actually use their queryset
        a_kwargs = {}
        if self.has_input() and self.uses_parent_queryset:
            if self._has_common_variants():
                kwargs["common_variants"] = True

            # Pass annotation gnomAD version so common filter version can be validated
            # @see https://github.com/SACGF/variantgrid/issues/1119
            if "annotation_gnomad_version" not in kwargs:
                try:
                    kwargs["annotation_gnomad_version"] = self.analysis.annotation_version.variant_annotation_version.gnomad
                except AttributeError:
                    pass

            for parent in self.get_non_empty_parents():
                a_kwargs.update(parent.get_annotation_kwargs(**kwargs))

        kwargs["existing_annotation_kwargs"] = set(a_kwargs.keys())
        a_kwargs.update(self._get_annotation_kwargs_for_node(**kwargs))
        return a_kwargs

    @property
    def queryset_requires_distinct(self):
        if self._queryset_requires_distinct():
            return True

        if self.has_input() and self.uses_parent_queryset:
            for parent in self.get_non_empty_parents():
                if parent.queryset_requires_distinct:
                    return True
        return False

    def _queryset_requires_distinct(self):
        """ Override if you need this - don't do by default as it's slow """
        return False

    @staticmethod
    def q_all():
        return Q(pk__isnull=False)

    @staticmethod
    def q_none():
        return ~AnalysisNode.q_all()

    def _get_cache_key(self) -> str:
        return str(self.node_version.pk)

    def _get_arg_q_dict_from_parents_and_node(self):
        arg_q_dict = self.get_parent_arg_q_dict()
        if self.modifies_parents():
            if node_arg_q_dict := self._get_node_arg_q_dict():
                self.merge_arg_q_dicts(arg_q_dict, node_arg_q_dict)
        return arg_q_dict

    def get_arg_q_dict(self, disable_cache=False) -> dict[Optional[str], dict[str, Q]]:
        """ A Django Q object representing the Variant filters for this node.
            This is the method to override in subclasses - not get_queryset()

            @see https://github.com/SACGF/variantgrid/wiki/Analysis-Nodes#node-q-objects
        """
        # We need this for node counts, and doing a grid query (each page) - and it can take a few secs to generate
        # for some nodes (Comp HET / pheno) so cache it
        cache_key = self._get_cache_key() + f"q_cache={disable_cache}"
        arg_q_dict: dict[Optional[str], dict[str, Q]] = {}
        if self._cache_node_q:
            arg_q_dict = cache.get(cache_key)

        if not arg_q_dict:
            if disable_cache is False:
                if cache_arg_q_dict := self._get_node_cache_arg_q_dict():
                    return cache_arg_q_dict

            if self.has_input():
                arg_q_dict = self._get_arg_q_dict_from_parents_and_node()
            else:
                if node_arg_q_dict := self._get_node_arg_q_dict():
                    arg_q_dict = node_arg_q_dict
                else:
                    arg_q_dict = {None: {}}
            if self._cache_node_q:
                try:
                    cache.set(cache_key, arg_q_dict)
                except Exception:
                    log_traceback()
        return arg_q_dict

    def get_contigs(self) -> set[Contig]:
        """ A set of contigs that contain variants for the node """

        cache_key = self._get_cache_key() + "_contigs"
        contigs: set[Contig] = cache.get(cache_key)

        if contigs is None:
            if self.has_input():
                contigs = self.get_parent_contigs()
                if self.modifies_parents():
                    node_contigs = self._get_narrowing_contigs()
                    if node_contigs is not None:
                        contigs &= node_contigs
            else:
                node_contigs = self._get_narrowing_contigs()
                if node_contigs is not None:
                    contigs = node_contigs
                else:
                    contigs = set(self.analysis.genome_build.contigs)
            cache.set(cache_key, contigs)
        return contigs

    def get_parent_arg_q_dict(self):
        if self.min_inputs == 1:
            return self.get_single_parent_arg_q_dict()
        raise NotImplementedError("Implement a non-default 'get_parent_arg_q_dict' if you have more than 1 parent")

    def get_parent_contigs(self) -> set[Contig]:
        contigs = set()
        for parent in self.get_non_empty_parents():
            contigs.update(parent.get_contigs())
        return contigs

    @property
    def use_cache(self):
        return False

    def write_cache(self, variant_collection: VariantCollection):
        qs = self.get_queryset(disable_cache=True)
        qs = qs.annotate(variant_collection_id=Value(variant_collection.pk, output_field=IntegerField()))
        sql = queryset_to_sql(qs.values_list('pk', 'variant_collection_id'))
        write_sql_to_variant_collection(variant_collection, sql)

    @cached_property
    def columns_version(self):
        return self.analysis.annotation_version.variant_annotation_version.columns_version

    @cached_property
    def node_version(self):
        return NodeVersion.get(self)

    @cached_property
    def node_cache(self) -> Optional['NodeCache']:
        if parent := self.get_unmodified_single_parent_node():
            return parent.node_cache
        return NodeCache.objects.filter(node_version=self.node_version,
                                        variant_collection__status=ProcessingStatus.SUCCESS).first()

    def _get_node_cache_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = {}
        if self.node_cache:
            arg_q_dict = self.node_cache.variant_collection.get_arg_q_dict()
        return arg_q_dict

    @staticmethod
    def get_cached_node_pks(node) -> Optional[list[int]]:
        """ The exact PK set stored at load for nodes <= ANALYSIS_NODE_STORE_ID_SIZE_MAX (@see node_counts),
            or None for a large node - or one loaded before the PKs were stored with the count """
        try:
            node_count = NodeCount.load_for_node(node, BuiltInFilters.TOTAL)
        except (NodeCount.DoesNotExist, NodeVersion.DoesNotExist):
            return None
        return node_count.variant_ids

    @staticmethod
    def get_small_parent_arg_q_dict(parent) -> Optional[dict[Optional[str], dict[str, Q]]]:
        """ Issue #546 explicit-PK substitution. When the parent holds only a small number of variants,
            substitute its contribution with a literal Q(pk__in=[...]) of the PKs it stored at load, so
            Postgres plans a tight bitmap-or over the variant PK index instead of re-running the parent's
            full filter chain wrapped in pk IN (subquery).

            Returns None when the parent has no stored PKs, in which case callers fall back to
            parent.get_arg_q_dict() - a subquery, which is always self-consistent. """
        if (variant_ids := AnalysisNode.get_cached_node_pks(parent)) is not None:
            q = Q(pk__in=variant_ids)
            return {None: {q: q}}
        return None

    def _get_node_q(self) -> Optional[Q]:
        return None

    def get_extra_filters_tag_q(self, tag_ids: list[str]) -> Optional[Q]:
        """ Tag scope for a tag 'extra_filters' selection, or None to use the analysis-scoped
            default @see get_node_extra_filters_q """
        return None

    @staticmethod
    def merge_arg_q_dicts(arg_q_dict, other_arg_q_dict):
        for k, other_q_dict in other_arg_q_dict.items():
            existing_dict = arg_q_dict.get(k, {})
            existing_dict.update(other_q_dict)
            arg_q_dict[k] = existing_dict

    def _get_node_q_hash(self) -> str:
        """" A Hash such that the same value equals the same Q filter being applied
             This is so merge node can remove duplicate filters - Q objects that use querysets don't hash the same
             Default implementation is to use something unique so will never merge them
        """
        return self.get_identifier()

    def _get_node_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        """ By default - we assume node implements _get_node_q and none of the filters apply to annotations """
        node_arg_q_dict = {}
        if node_q := self._get_node_q():
            node_arg_q_dict[None] = {self._get_node_q_hash(): node_q}
        return node_arg_q_dict

    def _get_node_contigs(self) -> Optional[set[Contig]]:
        """ Return the contigs we filter for in this node. None means we don't know how to describe that """
        return None

    def _get_narrowing_contigs(self) -> Optional[set[Contig]]:
        """ _get_node_contigs + always letting the gene-level contig through.

            Nodes can restrict queries to contigs (though we don't actually use this yet)
            Gene-level variants are on a contig of their own so should always be included
            @see snpdb.gene_level_variants """
        node_contigs = self._get_node_contigs()
        if node_contigs is not None:
            node_contigs = node_contigs | {Contig.get_gene_level()}
        return node_contigs

    def get_queryset(self, extra_filters_q=None, extra_annotation_kwargs=None, arg_q_dict=None,
                     inner_query_distinct=False, disable_cache=False):
        if extra_annotation_kwargs is None:
            extra_annotation_kwargs = {}

        qs = self._get_model_queryset()
        a_kwargs = self.get_annotation_kwargs()
        a_kwargs.update(extra_annotation_kwargs)
        if arg_q_dict is None:
            arg_q_dict = self.get_arg_q_dict(disable_cache=disable_cache)
            # print(arg_q_dict)

        qs, q_list = annotate_and_filter_queryset(qs, a_kwargs, arg_q_dict)

        if self.analysis.node_queryset_filter_contigs:
            q_list.append(Q(locus__contig__in=self.get_contigs()))

        if extra_filters_q:
            q_list.append(extra_filters_q)
        if q_list:
            q = reduce(operator.and_, q_list)
            filtered_qs = qs.filter(q)

            if self.queryset_requires_distinct:
                if inner_query_distinct:
                    qs = qs.filter(pk__in=filtered_qs.values_list("pk", flat=True))
                else:
                    qs = filtered_qs.distinct()
            else:
                qs = filtered_qs

        # Clear ordering, @see
        # https://docs.djangoproject.com/en/3.0/topics/db/aggregation/#interaction-with-default-ordering-or-order-by
        return qs.order_by()

    def get_grid_post_data(self) -> dict:
        """ Per-request state the node grid page sends back as its ajax params """
        return {}

    def get_class_name(self):
        return self.__class__.__name__

    def get_identifier(self):
        return f"{self.get_class_name()}-{self.pk}"

    def get_css_classes(self):
        """ returns list of css classes - set on "node > .node-overlay" on node appearance update """
        css_classes = []
        if self.output_node:
            css_classes.append("output-node")
        if self.analysis.template_type == AnalysisTemplateType.TEMPLATE and self.analysisvariable_set.exists():
            css_classes.append("variable-node")
        return css_classes

    def get_input_count(self):
        parents = self.get_non_empty_parents()
        return sum(p.get_output_count() for p in parents)

    def get_output_count(self):
        # TODO: Move the if not modify parents code in here.
        if self.count is not None:
            return self.count
        count = self.get_queryset().count()
        self.count = count
        self.save()
        return count

    def _get_method_summary(self):
        raise NotImplementedError()

    def get_method_summary(self):
        errors = self.get_errors(flat=True)
        if not errors:
            html_summary = self._get_method_summary()
        else:
            html_summary = "<b>incorrectly configured</b><ul>"
            html_summary += "".join([f"<li>{error}</li>" for error in errors])
            html_summary += "</ul>"
        return html_summary

    def get_node_name(self):
        """ Automatic node name """
        raise NotImplementedError(f"Node Class: {self.get_class_name()}")

    def get_name_or_identifier(self):
        return self.get_node_name() or self.get_identifier()

    @staticmethod
    def get_help_text() -> str:
        raise NotImplementedError()

    @staticmethod
    def get_node_class_label():
        """ Used in create node dropdown """
        raise NotImplementedError("get_node_class_label not implemented - this is probably due to a new class, or a reverse migration wiping out the subclass leaving just the AnalysisNode")

    @classmethod
    def get_node_class_label_short(cls) -> str:
        """ Class strip on the node card - override where the long label doesn't fit 128px """
        return cls.get_node_class_label()

    def get_node_strip_label(self) -> str:
        """ Class strip text for this node - override where a node's configuration changes what it
            is rather than just what it reads (eg SampleNode at extraction level) """
        return self.get_node_class_label_short()

    @classmethod
    def get_node_class_icon(cls) -> NodeIcon:
        """ Badge icon for a node of this class - also what the create node dropdown shows """
        return NodeIcon(fa="fa-solid fa-circle-nodes")

    @classmethod
    def get_menu_entries(cls) -> list[NodeMenuEntry]:
        """ The add node dropdown rows this class provides. One per class unless a class is
            configurable enough that users look for its configurations by name """
        return [NodeMenuEntry(key=cls.__name__, label=cls.get_node_class_label(),
                              icon=cls.get_node_class_icon())]

    def get_node_icon(self) -> NodeIcon:
        """ Badge icon for this node - the class default unless config changes it (eg SampleNode
            draws the patient's pedigree shape) """
        return self.get_node_class_icon()

    def get_node_chips(self) -> list[NodeChip]:
        """ Pills under the node name - what this node is reading, from its saved config """
        return []

    def _get_genome_build_errors(self, field_name, field_genome_build: GenomeBuild) -> list:
        """ Used to quickly add errors about genome build mismatches
            This only happens in templates (ran template on sample with different build than hardcoded data)
            In normal analyses, autocomplete restrictions should not allow you to configure data from other builds """
        errors = []
        if field_genome_build != self.analysis.genome_build:
            msg = f"{field_name} genome build: {field_genome_build} different from analysis build: {self.analysis.genome_build}"
            errors.append(msg)
        return errors

    def _get_configuration_errors(self) -> list:
        return []

    def get_parents_and_errors(self):
        """ Returns error array, includes any min/max parent error and node config error """
        if self.has_input():
            return self.get_parent_subclasses_and_errors()
        return [], []

    def _get_analysis_errors(self) -> list[str]:
        if self._cached_analysis_errors is None:
            self._cached_analysis_errors = self.analysis.get_errors()
        return self._cached_analysis_errors

    def get_analysis_errors(self, flat=False):
        """ Analysis level errors (eg no annotation version) - these can't be fixed by editing the node """
        errors = [(NodeErrorSource.ANALYSIS, e) for e in self._get_analysis_errors()]
        if flat:
            errors = AnalysisNode.flatten_errors(errors)
        return errors

    def get_ignorable_error_fields(self) -> list[str]:
        """ Fields whose _get_field_errors() the user may waive with ignore_field_errors """
        return []

    def _get_field_errors(self) -> dict[str, list[str]]:
        """ Configuration errors keyed by the node field at fault, so the editor shows them against
            that field. Templates are exempt - they're configured however the author likes """
        return {}

    def _get_ignored_error_fields(self) -> set[str]:
        if self.ignore_field_errors:
            return set(self.get_ignorable_error_fields())
        return set()

    def get_field_errors(self) -> dict[str, list[str]]:
        """ The field errors that stop the node running """
        if self.analysis.template_type == AnalysisTemplateType.TEMPLATE:
            return {}
        ignored = self._get_ignored_error_fields()
        return {field: errors for field, errors in self._get_field_errors().items()
                if errors and field not in ignored}

    def get_ignored_field_errors(self) -> list[str]:
        if self.analysis.template_type == AnalysisTemplateType.TEMPLATE:
            return []
        ignored = self._get_ignored_error_fields()
        return [e for field, errors in self._get_field_errors().items() if field in ignored for e in errors]

    def get_warnings(self) -> list[str]:
        return list(self.get_ignored_field_errors())

    def _get_non_field_errors(self, include_parent_errors=True) -> list[tuple[NodeErrorSource, str]]:
        errors = self.get_analysis_errors()
        _, parent_errors = self.get_parents_and_errors()
        if include_parent_errors:
            errors.extend(parent_errors)
        if self.errors:
            errors.append((NodeErrorSource.INTERNAL_ERROR, self.errors))
        errors.extend((NodeErrorSource.CONFIGURATION, ce) for ce in self._get_configuration_errors())
        return errors

    def get_errors(self, include_parent_errors=True, flat=False):
        """ returns a tuple of (NodeError, str) unless flat=True where it's only string """
        errors = self._get_non_field_errors(include_parent_errors=include_parent_errors)
        errors.extend((NodeErrorSource.CONFIGURATION, e) for field_errors in self.get_field_errors().values()
                      for e in field_errors)
        if flat:
            errors = AnalysisNode.flatten_errors(errors)
        return errors

    def can_ignore_errors(self) -> bool:
        """ Everything wrong with the node is a field error ignore_field_errors would waive - what the
            Template tab offers when it reveals a branch the template run hid """
        field_errors = self.get_field_errors()
        if not field_errors or self._get_non_field_errors():
            return False
        return set(field_errors) <= set(self.get_ignorable_error_fields())

    @staticmethod
    def flatten_errors(errors):
        return [f"{NodeErrorSource(nes).label}: {error}" for nes, error in errors]

    @staticmethod
    def get_status_from_errors(errors):
        ERROR_STATUS = {
            NodeErrorSource.INTERNAL_ERROR: NodeStatus.ERROR,
            NodeErrorSource.ANALYSIS: NodeStatus.ERROR_WITH_PARENT,
            NodeErrorSource.PARENT: NodeStatus.ERROR_WITH_PARENT,
            NodeErrorSource.CONFIGURATION: NodeStatus.ERROR_CONFIGURATION,
        }
        if not errors:
            raise ValueError("Passed in empty errors!")

        error_sources = {s for s, _ in errors}
        for source, status in ERROR_STATUS.items():
            if source in error_sources:
                return status
        raise ValueError("No error source found")

    @staticmethod
    def throw_errors_exception(errors):
        ERROR_EXCEPTIONS = {
            NodeErrorSource.INTERNAL_ERROR: ValueError,
            NodeErrorSource.ANALYSIS: NonFatalNodeError,
            NodeErrorSource.PARENT: NodeParentErrorsException,
            NodeErrorSource.CONFIGURATION: NodeConfigurationException,
        }
        if not errors:
            raise ValueError("Passed in empty errors!")
        error_sources = {s for s, _ in errors}
        for source, exception_klass in ERROR_EXCEPTIONS.items():
            if source in error_sources:
                raise exception_klass()
        raise ValueError("No error source found")

    def inherits_parent_columns(self):
        return self.min_inputs == 1 and self.max_inputs == 1

    def _get_node_extra_columns(self) -> list[RichColumn]:
        """ Subclasses override to add their own columns to the node grid """
        return []

    def _get_inherited_columns(self) -> list[RichColumn]:
        extra_columns = []
        if self.inherits_parent_columns():
            parent = self.get_single_parent()
            extra_columns.extend(parent.get_extra_columns())
        return extra_columns

    def get_extra_columns(self) -> list[RichColumn]:
        extra_columns = []
        if self.is_valid:
            extra_columns.extend(self._get_inherited_columns())
        # Only add columns that are unique, as otherwise filters get added twice.
        for col in self._get_node_extra_columns():
            if col not in extra_columns:
                extra_columns.append(col)
        return extra_columns

    def get_node_classification(self):
        if self.is_source:
            classification = "source"
        else:
            classification = "filter"
        return classification

    def has_input(self):
        return self.max_inputs != 0

    @property
    def is_source(self):
        return self.has_input() is False

    @property
    def is_valid(self):
        return not self.get_errors()

    @property
    def is_ready(self):
        return NodeStatus.is_ready(self.status)

    def bump_version(self):
        self.version += 1
        self.status = NodeStatus.DIRTY
        self.count = None
        self.errors = None
        self.cloned_from = None

    def modifies_parents(self):
        """ Can overwrite and set to False to use parent counts """
        return True

    def parent_count(self) -> int:
        return AnalysisEdge.objects.filter(child=self).count()

    def get_unmodified_single_parent_node(self) -> Optional['AnalysisNode']:
        """ If a node doesn't modify single parent - can use that in some places to re-use cache """

        if self.has_input() and self.parent_count() == 1 and self.is_valid and not self.modifies_parents():
            try:
                return self.get_single_parent()
            except ValueError:
                pass
        return None

    def _get_cached_label_count(self, label) -> Optional[int]:
        """ Override for optimisation.
            Returning None means we need to run the SQL to get the count """

        try:
            if self.cloned_from:
                # If cloned (and we or original haven't changed) - use those counts
                try:
                    node_count = NodeCount.load_for_node_version(self.cloned_from, label)
                    return node_count.count
                except NodeCount.DoesNotExist:
                    # Should only ever happen if original bumped version since we were loaded
                    # otherwise should have cascade set cloned_from to NULL
                    pass

            if self.has_input():
                parent_non_zero_label_counts = []
                for parent in self.get_non_empty_parents():
                    if parent.count != 0:  # count=0 has 0 for all labels
                        parent_node_count = NodeCount.load_for_node(parent, label)
                        if parent_node_count.count != 0:
                            parent_non_zero_label_counts.append(parent_node_count.count)

                if not parent_non_zero_label_counts:
                    # logging.info("all parents had 0 %s counts", label)
                    return 0

                if not self.modifies_parents():
                    if len(parent_non_zero_label_counts) == 1:
                        # logging.info("Single parent, no modification, using that")
                        return parent_non_zero_label_counts[0]
        except NodeCount.DoesNotExist:
            pass
        except Exception as e:
            logging.warning("Trouble getting cached %s count: %s", label, e)

        return None

    def get_grid_node_id_and_version(self):
        """ Uses parent node_id/version if possible to re-use cache """
        node_id = self.pk
        version = self.version
        if self.cloned_from:
            node_id = self.cloned_from.node_id
            version = self.cloned_from.version

        if parent := self.get_unmodified_single_parent_node():
            node_id, version = parent.get_grid_node_id_and_version()
        return node_id, version

    def node_counts(self):
        """ This is inside Celery task """

        self.count = None
        # Record provenance first, so the checks below know whether this node's data can move under it
        live_data_sources = self.get_live_data_sources()
        NodeVersion.objects.filter(pk=self.node_version.pk).update(live_data_sources=live_data_sources)
        self.node_version.live_data_sources = live_data_sources

        counts_to_get = {BuiltInFilters.TOTAL}
        counts_to_get.update([i[0] for i in self.analysis.get_node_count_types()])
        label_counts = {}

        for label in counts_to_get:
            label_count = self._get_cached_label_count(label)
            if label_count is not None:
                label_counts[label] = label_count

        counts_to_get -= set(label_counts)

        logging.debug("Node %d.%d cached counts: %s", self.pk, self.version, label_counts)
        if counts_to_get:
            logging.debug("Node %d.%d needs DB request for %s", self.pk, self.version, counts_to_get)
            retrieved_label_counts = get_node_counts_and_labels_dict(self, counts_to_get)
            label_counts.update(retrieved_label_counts)

        total_count = label_counts[BuiltInFilters.TOTAL]
        variant_ids = self._get_variant_ids_to_store(total_count)
        if variant_ids is not None:
            # For a small node the PK list is the truth - it and the count came from the same load
            total_count = len(variant_ids)

        node_counts = []
        for label, count in label_counts.items():
            node_counts.append(NodeCount(node_version=self.node_version, label=label, count=count))
        if node_counts:
            total_node_count = next(nc for nc in node_counts if nc.label == BuiltInFilters.TOTAL)
            total_node_count.count = total_count
            total_node_count.variant_ids = variant_ids
            # Counts are a cache of a query against an immutable node_version, so a re-load (eg after a
            # backoff retry that failed once these were already written) can safely overwrite them
            NodeCount.objects.bulk_create(node_counts, update_conflicts=True,
                                          update_fields=["count", "variant_ids", "modified"],
                                          unique_fields=["node_version", "label"])

        # Every label count is a subset of the total - a bigger one means the query fanned out over a
        # multi-valued join, or a cached count is out of sync with the live query
        if bigger_than_total := {l: c for l, c in label_counts.items() if c > total_count}:
            self._raise_or_warn_count_mismatch(f"label counts {bigger_than_total} > total count={total_count}")

        # Single parent nodes should always reduce the number of variants - run a check to make sure the
        # query wasn't bad and returned more results than it should have
        parents = list(self.get_non_empty_parents())
        if len(parents) == 1:
            parent = parents[0]
            if parent.count < total_count:
                self._raise_or_warn_count_mismatch(f"count={total_count} > {parent=}(pk={parent.pk}) count={parent.count}")

        return NodeStatus.READY, total_count

    def _get_variant_ids_to_store(self, total_count: int) -> Optional[list[int]]:
        """ The node's exact PK set, for nodes small enough to hold it - taken after the count so a node
            whose data moved between the two is caught here rather than storing a silently short list """
        max_size = settings.ANALYSIS_NODE_STORE_ID_SIZE_MAX
        if not (max_size and total_count <= max_size):
            return None

        variant_ids = list(self.get_queryset().values_list("pk", flat=True)[:max_size + 1])
        if len(variant_ids) > max_size:
            self._raise_or_warn_count_mismatch(f"{len(variant_ids)} pks > max_size={max_size}")
            return None
        return variant_ids

    def _load(self):
        """ Override to do anything interesting.
            Return a dict of node fields to write as part of the load - setting them on self isn't
            enough, as load() persists via update() (see #431 - no save() in celery tasks) """

    def load(self):
        """ load is called after parents are run """
        # logging.debug("node %d (%d) load()", self.id, self.version)
        start = time()
        load_update_kwargs = self._load() or {}  # Do before counts in case it affects anything
        status, count = self.node_counts()
        load_seconds = time() - start
        self.update(status=status, count=count, load_seconds=load_seconds, **load_update_kwargs)

    def add_parent(self, parent, *args, **kwargs):
        if not parent.visible:
            raise NonFatalNodeError("Not connecting children to invisible nodes!")

        existing_connect = parent.children.through.objects.filter(parent=parent, child=self)
        if not existing_connect.exists():
            super().add_parent(parent)
            self.parents_changed = True
        else:
            logging.error("Node(pk=%d).add_parent(pk=%d) already exists!", self.pk, parent.pk)

    def remove_parent(self, parent):
        """ disconnects parent by deleting edge """
        # Ok to have multiple, just delete first
        edge = parent.children.through.objects.filter(parent=parent, child=self).first()
        if edge:  # could be some kind of race condition?
            edge.delete()
        self.parents_changed = True

    def handle_ancestor_input_samples_changed(self):
        pass

    def update(self, **kwargs):
        """ Updates Node if self.version matches DB - otherwise throws NodeOutOfDateException """
        # Subclass queryset so node subclass fields (eg has_gene_coverage) can be written as well as base ones
        self_qs = type(self).objects.filter(pk=self.pk, version=self.version)
        updated = self_qs.update(**kwargs)
        if not updated:
            raise NodeOutOfDateException()

    def save(self, *args, **kwargs):
        """ To avoid race conditions, don't use save() in a celery task (unless running in scheduling_single_worker)
            instead use update() method above """
        with transaction.atomic():
            # Lock the analysis row first so concurrent node-tree saves serialize in a consistent order.
            # Saving a node cascades version bumps to its children (see below), and each NodeVersion insert
            # takes a FOR KEY SHARE lock on the referenced analysisnode row. Without this, two concurrent
            # saves over overlapping subtrees grab those row locks in opposite orders and deadlock.
            Analysis.objects.select_related(None).select_for_update().get(pk=self.analysis_id)
            return self._save(*args, **kwargs)

    def _save(self, *args, **kwargs):
        # logging.debug("save: pk=%s kwargs=%s", self.pk, str(kwargs))
        super_save = super().save

        if self.parents_changed or self.ancestor_input_samples_changed:
            self.handle_ancestor_input_samples_changed()

        if self.auto_node_name:
            self.name = self.get_node_name()

        # TODO: This causes lots of DB queries... should we change this?
        self.valid = self.is_valid
        if not self.valid:
            self.shadow_color = NodeColors.ERROR
            self.appearance_dirty = True
        elif self.shadow_color == NodeColors.ERROR:  # Need to allow nodes to set to warning
            self.shadow_color = NodeColors.VALID
            self.appearance_dirty = True

        if self.appearance_dirty:
            self.appearance_version += 1

        if self.parents_changed or self.queryset_dirty:
            self.bump_version()

            super_save(*args, **kwargs)
            if self.update_children:
                # We also need to bump if node has it's own sample - as in templates, we set fields in toposort order
                # So we could go from having multiple proband samples to only one later (thus can set descendants)
                for kid in self.children.select_subclasses():
                    kid.ancestor_input_samples_changed = self.is_source or self.ancestor_input_samples_changed or \
                                                         self.get_samples_from_node_only_not_ancestors()
                    kid.appearance_dirty = False
                    kid.queryset_dirty = True
                    kid.save()  # Will bump versions
        else:
            super_save(*args, **kwargs)

        # Make sure this always exists
        NodeVersion.objects.get_or_create(node=self, version=self.version)
        # Modify our analyses last updated time
        Analysis.objects.filter(pk=self.analysis.pk).update(modified=timezone.now())

    def set_node_task_and_status(self, celery_task, status):
        with connection.cursor() as cursor:
            db_pid = get_backend_pid(cursor)
        self.update(status=status)

        NodeTask.objects.filter(node_version__node=self, node_version__version=self.version) \
            .update(celery_task=celery_task, db_pid=db_pid)

    def claim_for_load(self, celery_task) -> bool:
        """ Take ownership of loading this node/version, returning False if someone else has it.

            A lease that expires while its update task is still sitting in the analysis_workers
            queue is reclaimed and re-dispatched, so two update tasks can exist for the same
            node/version. The conditional UPDATE is the arbiter - only one of them moves the node
            out of CLAIMABLE_STATUSES, and the loser leaves the winner's status/lease alone. """

        claimed = AnalysisNode.objects.filter(pk=self.pk, version=self.version,
                                              status__in=NodeStatus.CLAIMABLE_STATUSES) \
            .update(status=NodeStatus.LOADING)
        if not claimed:
            return False
        self.status = NodeStatus.LOADING
        with connection.cursor() as cursor:
            db_pid = get_backend_pid(cursor)
        # Restart the lease window now the load is actually starting - a node dispatched into a
        # backlog would otherwise spend most of its lease queued, and be perma-failed as a lost
        # worker part way through a perfectly healthy load
        NodeTask.objects.filter(node_version__node=self, node_version__version=self.version) \
            .update(celery_task=celery_task, db_pid=db_pid,
                    lease_expires=timezone.now() + timedelta(seconds=LEASE_SECONDS))
        return True

    def adjust_cloned_parents(self, old_new_map):
        """ If you need to do something with old/new parents """

    def save_clone(self):
        with disable_auditlog():
            node_id = self.pk
            try:
                # Have sometimes had race condition where we try to clone a node that has been updated
                # In that case we'll just miss out on the cache
                original_node_version = self.node_version
            except NodeVersion.DoesNotExist:
                original_node_version = None

            copy = self
            # Have to set both id/pk to None when using model inheritance
            copy.id = None
            copy.pk = None
            copy.version = 1  # 0 is for those being constructed in analysis templates
            # Store cloned_from so we can use original's NodeCounts
            copy.cloned_from = original_node_version
            copy.save()

            for npf in NodeVCFFilter.objects.filter(node_id=node_id):
                npf.pk = None
                npf.node = copy
                npf.save()

            naff = NodeAlleleFrequencyFilter.objects.filter(node_id=node_id).first()  # 1-to-1
            if naff:
                af_frequency_ranges = list(naff.nodeallelefrequencyrange_set.all().values_list("min", "max"))
                # Use existing if already created for node (eg AlleleFrequencyNode always makes one)
                copy_naff, created = NodeAlleleFrequencyFilter.objects.get_or_create(node=copy)
                if not created:
                    # Wipe out defaults to clear way for clone
                    copy_naff.nodeallelefrequencyrange_set.all().delete()
                copy_naff.group_operation = naff.group_operation
                copy_naff.save()

                for min_value, max_value in af_frequency_ranges:
                    copy_naff.nodeallelefrequencyrange_set.create(min=min_value, max=max_value)

            return copy

    def __str__(self):
        return self.name

    @classmethod
    def depth_first(cls, node):
        parents = node.get_parent_subclasses()
        l = []
        for p in parents:
            l.extend(cls.depth_first(p))
        l.append(node)
        return l


class AnalysisEdge(NodeAuditLogMixin, edge_factory(AnalysisNode, concrete=False)):
    def _get_node(self):
        return self.child

    def get_additional_data(self):
        """ For django-audit-log """

        additional_data = super().get_additional_data()
        additional_data.update({
            "parent": {
                "id": self.parent.pk,
                "content_type": get_model_content_type_dict(self.parent),
                "object_repr": str(self.parent),
            },
            "child": {
                "id": self.child.pk,
                "content_type": get_model_content_type_dict(self.child),
                "object_repr": str(self.child),
            }
        })
        return additional_data


class NodeTask(TimeStampedModel):
    """ Tracks/locks the celery update task for a single node-version, and acts as the lease
        record for the state-driven scheduler (issue #346, mocha-style): the row is claimed by a
        dispatcher (lease_ready_nodes), stamped with a lease_expires, and reclaimed if the owning
        worker dies.

        OneToOne against NodeVersion so it shares that row's lifecycle - the (node, version)
        identity comes from the NodeVersion, and stale tasks are cleaned up automatically when old
        versions are deleted (delete_analysis_old_node_versions cascades). """

    node_version = models.OneToOneField("NodeVersion", on_delete=CASCADE)
    analysis_update_uuid = models.UUIDField()
    celery_task = models.CharField(max_length=36, null=True)
    db_pid = models.IntegerField(null=True)
    # Lease + backoff fields (mocha lease record):
    leased_by = models.CharField(max_length=64, null=True)
    lease_expires = models.DateTimeField(null=True)
    attempt_count = models.IntegerField(default=0)
    last_attempt = models.DateTimeField(null=True)
    run_after = models.DateTimeField(null=True)  # backoff gate: not leasable before this

    def __str__(self):
        return f"NodeTask: {self.analysis_update_uuid} - {self.node_version}"


class NodeWiki(Wiki):
    node = models.OneToOneField(AnalysisNode, on_delete=CASCADE)

    def _get_restricted_object(self):
        return self.node.analysis


class AnalysisNodeAlleleSource(AlleleSource):
    """ Used to link a nodes variants to alleleles and liftover to other builds """
    node = models.ForeignKey(AnalysisNode, null=True, on_delete=SET_NULL)

    def get_genome_build(self):
        if self.node:
            genome_build = self.node.analysis.genome_build
        else:
            genome_build = None
        return genome_build

    def get_variant_qs(self):
        if self.node:
            qs = self.node.get_subclass().get_queryset()
        else:
            qs = Variant.objects.none()
        return qs


class NodeVersion(TimeStampedModel):
    """ This will be deleted once a node updates, so make all version specific caches cascade delete from this """
    node = models.ForeignKey(AnalysisNode, on_delete=CASCADE)
    version = models.IntegerField(null=False)
    # {source_key: data_version} of the mutable tables this node read at load. Empty = deterministic
    live_data_sources = models.JSONField(default=dict)

    class Meta:
        unique_together = ("node", "version")

    @staticmethod
    def get(node: AnalysisNode):
        try:
            return NodeVersion.objects.get(node=node, version=node.version)
        except NodeVersion.DoesNotExist:
            node.check_still_valid()
            raise

    def __str__(self):
        return f"{self.node.pk} (v{self.version})"


class NodeCache(models.Model):
    """ This is only used by Intersection node now """
    node_version = models.OneToOneField(NodeVersion, on_delete=CASCADE)
    variant_collection = models.OneToOneField(VariantCollection, on_delete=CASCADE)

    @staticmethod
    def get_or_create_for_node(node: AnalysisNode) -> tuple['NodeCache', bool]:
        variant_collection = VariantCollection.objects.create(name=f"NodeCache {node.node_version}")
        defaults = {"variant_collection": variant_collection}
        node_cache, created = thread_safe_unique_together_get_or_create(NodeCache, node_version=node.node_version,
                                                                        defaults=defaults)
        if not created:
            variant_collection.delete()
        return node_cache, created

    def __str__(self):
        return f"NodeCache {self.node_version}: {self.variant_collection.get_status_display()}"


@receiver(post_delete, sender=NodeCache)
def post_delete_node_cache(sender, instance, **kwargs):  # pylint: disable=unused-argument
    """ This can sometimes be called multiple times - if node updated again before previous delete is done """
    try:
        if instance.variant_collection:
            instance.variant_collection.delete_related_objects()
            instance.variant_collection.delete()
    except VariantCollection.DoesNotExist:
        # Deleted already
        pass


class NodeCount(TimeStampedModel):
    node_version = models.ForeignKey(NodeVersion, on_delete=CASCADE)
    label = models.CharField(max_length=100)
    count = models.IntegerField(null=False)
    # Only on the TOTAL row, and only for nodes <= ANALYSIS_NODE_STORE_ID_SIZE_MAX. When present,
    # count == len(variant_ids) and this is the exact set the node held at load
    variant_ids = ArrayField(models.IntegerField(), null=True)

    class Meta:
        unique_together = ("node_version", "label")

    @staticmethod
    def load_for_node_version(node_version: NodeVersion, label: str) -> 'NodeCount':
        return NodeCount.objects.get(node_version=node_version, label=label)

    @staticmethod
    def load_for_node(node: AnalysisNode, label: str) -> 'NodeCount':
        return NodeCount.load_for_node_version(node.node_version, label=label)

    def __str__(self):
        return f"NodeCount({self.node_version}, {self.label}) = {self.count}"


class NodeColumnSummaryCacheCollection(models.Model):
    node_version = models.ForeignKey(NodeVersion, on_delete=CASCADE)
    variant_column = models.TextField(null=False)
    extra_filters = models.TextField(null=False)

    @staticmethod
    def get_counts_for_node(node, variant_column, extra_filters):
        ncscc, created = NodeColumnSummaryCacheCollection.objects.get_or_create(node_version=node.node_version,
                                                                                variant_column=variant_column,
                                                                                extra_filters=extra_filters)
        if created:
            extra_filters_q = get_node_extra_filters_q(node, extra_filters)
            queryset = node.get_queryset(extra_filters_q)
            count_qs = queryset.values_list(variant_column).distinct().annotate(Count('id'))
            data_list = []
            for value, count in count_qs:
                data = NodeColumnSummaryData(collection=ncscc,
                                             value=value,
                                             count=count)
                data_list.append(data)

            if data_list:
                NodeColumnSummaryData.objects.bulk_create(data_list)
        else:
            data_list = ncscc.nodecolumnsummarydata_set.all()

        counts = {}
        for ncsd in data_list:
            counts[ncsd.value] = ncsd.count

        return counts


class NodeColumnSummaryData(models.Model):
    collection = models.ForeignKey(NodeColumnSummaryCacheCollection, on_delete=CASCADE)
    value = models.TextField(null=True)
    count = models.IntegerField(null=False)


class NodeVCFFilter(NodeAuditLogMixin, models.Model):
    """ If these exist, they mean use that filter """
    node = models.ForeignKey(AnalysisNode, on_delete=CASCADE)
    vcf_filter = models.ForeignKey(VCFFilter, on_delete=CASCADE, null=True)  # null = 'PASS'

    def _get_node(self):
        return self.node

    @staticmethod
    def get_filter_ids(node):
        """ Returns PASS as None """
        all_nvf_qs = NodeVCFFilter.objects.filter(node_id=node.pk)
        return set(all_nvf_qs.values_list("vcf_filter__filter_id", flat=True))

    @staticmethod
    def has_pass(node) -> bool:
        """ PASS is stored as the node level row with no vcf_filter - it's the one FILTER value
            that means the same thing in every VCF """
        return NodeVCFFilter.objects.filter(node_id=node.pk, vcf_filter__isnull=True).exists()

    @staticmethod
    def get_vcf_filter_ids(node, vcf=None) -> list[tuple]:
        """ (vcf_id, filter_id) for the node's non-PASS rows, optionally for one VCF """
        nvf_qs = NodeVCFFilter.objects.filter(node_id=node.pk, vcf_filter__isnull=False)
        if vcf is not None:
            nvf_qs = nvf_qs.filter(vcf_filter__vcf=vcf)
        return list(nvf_qs.values_list("vcf_filter__vcf_id", "vcf_filter__filter_id"))

    @staticmethod
    def get_filter_codes(node, vcf, pass_only: Optional[bool] = None):
        """ What this VCF lets through: its own ticked codes, plus PASS if the node's PASS row is set.

            PASS is the one FILTER value that means the same thing in every VCF, so a code is never
            translated across them - 'LowDepth' in a DRAGEN small variant VCF is not 'LowDepth' in
            its CNV VCF. pass_only stands in for the node's PASS row, which is how SampleNode lets
            one of its samples decide for itself. """
        nvf_qs = NodeVCFFilter.objects.filter(node_id=node.pk, vcf_filter__vcf=vcf)
        filter_codes = set(nvf_qs.values_list("vcf_filter__filter_code", flat=True))
        if pass_only is None:
            pass_only = NodeVCFFilter.has_pass(node)
        if pass_only:
            filter_codes.add(None)
        return filter_codes


class NodeAlleleFrequencyFilter(NodeAuditLogMixin, models.Model):
    """ Used for various nodes """
    node = models.OneToOneField(AnalysisNode, on_delete=CASCADE)
    group_operation = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ANY)

    def _get_node(self):
        return self.node

    def get_q(self, allele_frequency_path: str, allele_frequency_percent: bool) -> Optional[Q]:
        af_q = None
        try:
            filters = []
            for af_range in self.nodeallelefrequencyrange_set.all():
                # Only apply filter if restricted range.
                # Missing value (historical data) == -1 so those will come through
                and_filters = []
                if af_range.min > 0:
                    min_value = af_range.min
                    if allele_frequency_percent:
                        min_value *= 100.0
                    and_filters.append(Q(**{allele_frequency_path + "__gte": min_value}))
                if af_range.max < 1:
                    max_value = af_range.max
                    if allele_frequency_percent:
                        max_value *= 100.0
                    and_filters.append(Q(**{allele_frequency_path + "__lte": max_value}))

                if and_filters:
                    and_q = reduce(operator.and_, and_filters)
                    filters.append(and_q)
            if filters:
                group_op = GroupOperation.get_operation(self.group_operation)
                af_q = reduce(group_op, filters)
        except NodeAlleleFrequencyFilter.DoesNotExist:
            pass

        return af_q

    @staticmethod
    def get_sample_arg_q_dict(node: AnalysisNode, sample: Sample) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = {}
        if sample:
            try:
                alias, allele_frequency_path = sample.get_cohort_genotype_alias_and_field("allele_frequency")
                allele_frequency_percent = sample.vcf.allele_frequency_percent
                if af_q := node.nodeallelefrequencyfilter.get_q(allele_frequency_path, allele_frequency_percent):
                    arg_q_dict[alias] = {str(af_q): af_q}
            except NodeAlleleFrequencyFilter.DoesNotExist:
                pass
        return arg_q_dict

    def get_description(self):
        # TODO: do this properly with group operators etc
        af_ranges = list(self.nodeallelefrequencyrange_set.all())
        if len(af_ranges) == 1:
            description = str(af_ranges[0])
        else:
            description = f"{self.get_group_operation_display()} of {len(af_ranges)} filters"
        return description

    def __str__(self):
        return self.get_description()


class NodeAlleleFrequencyRange(NodeAuditLogMixin, models.Model):
    MIN_VALUE = 0
    MAX_VALUE = 1

    filter = models.ForeignKey(NodeAlleleFrequencyFilter, on_delete=CASCADE)
    min = models.FloatField(null=False)
    max = models.FloatField(null=False)

    def _get_node(self):
        return self.filter.node

    def __str__(self):
        has_min = self.min is not None and self.min > self.MIN_VALUE
        has_max = self.max is not None and self.max < self.MAX_VALUE

        min_perc = format_percent(self.min, is_unit=True)
        max_perc = format_percent(self.max, is_unit=True)

        if has_min and has_max:
            return f"{min_perc} - {max_perc}"
        if has_min:
            return f">={min_perc}"
        if has_max:
            return f"<={max_perc}"
        return ""


class AnalysisClassification(models.Model):
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE)
    classification = models.ForeignKey(Classification, on_delete=CASCADE)


auditlog.register(AnalysisEdge)
auditlog.register(NodeVCFFilter)
auditlog.register(NodeAlleleFrequencyFilter)
auditlog.register(NodeAlleleFrequencyRange)
