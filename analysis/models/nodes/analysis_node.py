""" AnalysisNode is the base class that all analysis nodes inherit from. """
import logging
import operator
from collections import defaultdict
from functools import cached_property, reduce
from random import random
from time import time
from typing import Tuple, Sequence, List, Dict, Optional, Set

from cache_memoize import cache_memoize
from celery.canvas import Signature
from django.conf import settings
from django.core.cache import cache
from django.db import connection, models
from django.db.models import Value, IntegerField
from django.db.models.aggregates import Count
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.query_utils import Q
from django.dispatch import receiver
from django.utils import timezone
from django_dag.models import node_factory, edge_factory
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager

from analysis.exceptions import NonFatalNodeError, NodeParentErrorsException, NodeConfigurationException, \
    NodeParentNotReadyException, NodeNotFoundException, NodeOutOfDateException
from analysis.models.enums import GroupOperation, NodeStatus, NodeColors, NodeErrorSource, AnalysisTemplateType
from analysis.models.models_analysis import Analysis
from analysis.models.nodes.node_counts import get_extra_filters_q, get_node_counts_and_labels_dict
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from classification.models import Classification, post_delete
from library.constants import DAY_SECS
from library.django_utils import thread_safe_unique_together_get_or_create
from library.log_utils import report_event, log_traceback
from library.utils import format_percent
from library.utils.database_utils import queryset_to_sql
from snpdb.models import BuiltInFilters, Sample, Variant, VCFFilter, Wiki, Cohort, VariantCollection, \
    ProcessingStatus, GenomeBuild, AlleleSource, Contig, SampleFilePath
from snpdb.variant_collection import write_sql_to_variant_collection


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


class AnalysisNode(node_factory('AnalysisEdge', base_model=TimeStampedModel)):
    model = Variant
    objects = NodeInheritanceManager()
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE)
    name = models.TextField(blank=True)
    x = models.IntegerField(default=_default_position)
    y = models.IntegerField(default=_default_position)
    version = models.IntegerField(default=0)  # Queryset version
    appearance_version = models.IntegerField(default=0)
    auto_node_name = models.BooleanField(default=True)
    output_node = models.BooleanField(default=False)
    hide_node_and_descendants_upon_template_configuration_error = models.BooleanField(default=False)
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

    def get_subclass(self):
        """ Returns the node loaded as a subclass """
        return AnalysisNode.objects.get_subclass(pk=self.pk)

    def check_still_valid(self):
        """ Checks that the node is still there and has the version we expect - or throw exception """
        version_qs = AnalysisNode.objects.filter(pk=self.pk).values_list("version", flat=True)
        if version_qs:
            db_version = version_qs[0]
            if db_version > self.version:
                raise NodeOutOfDateException()
        else:
            raise NodeNotFoundException(self.pk)

    def _get_cohorts_and_sample_visibility_for_node(self) -> Tuple[Sequence[Cohort], Dict]:
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

    def get_cohorts_and_sample_visibility(self, sort=True) -> Tuple[Sequence[Cohort], Dict]:
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
    def get_sample_ids(self) -> List[Sample]:
        return [s.pk for s in self.get_samples()]

    def get_samples_from_node_only_not_ancestors(self):
        cohorts, visibility = self._get_cohorts_and_sample_visibility_for_node()
        return self._get_visible_samples_from_cohort(cohorts, visibility)

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

    def get_samples(self) -> List[Sample]:
        """ Return all ancestor samples for a node"""
        cohorts, visibility = self.get_cohorts_and_sample_visibility(sort=False)
        return self._get_visible_samples_from_cohort(cohorts, visibility)

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

    def get_cache_task_args_objs_set(self, force_cache=False):
        """ returns Celery tasks which are called in node_utils.get_analysis_update_task before children are loaded
            Uses tasks not signatures so they are hashable in a set to be able to remove dupes """

        task_args_objs_set = set()
        if self.is_valid() and (force_cache or self.use_cache):
            if parent := self.get_unmodified_single_parent_node():
                return parent.get_cache_task_args_objs_set(force_cache=force_cache)

            node_cache, created = NodeCache.get_or_create_for_node(self)
            if created:
                task_args_objs_set.add((self.NODE_CACHE_TASK, (self.pk, self.version), node_cache))
            else:
                # Cache has been launched already, we just need to make sure it's ready, so launch a task
                # waiting on it, to be used as a dependency
                task_args_objs_set.add((self.WAIT_FOR_CACHE_TASK, (node_cache.pk, ), node_cache))
        return task_args_objs_set

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

    def get_non_empty_parents(self, require_parents_ready=True) -> List['AnalysisNode']:
        """ Returns non-empty (count > 0) parents.
            If require_parents_ready=True, die if parents not ready
            Otherwise, return them as we don't know if they're empty or not """
        non_empty_parents = []
        for p in self.get_parent_subclasses():
            if p.is_ready():
                if p.count == 0:
                    continue
            elif require_parents_ready:
                raise NodeParentNotReadyException(f"Parent {p} is not ready!")
            non_empty_parents.append(p)
        return non_empty_parents

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

    def get_single_parent_arg_q_dict(self) -> Dict[Optional[str], Dict[str, Q]]:
        arg_q_dict = {}
        parent = self.get_single_parent()
        if parent.is_ready():
            if parent.count == 0:
                q_none = self.q_none()
                arg_q_dict[None] = {str(q_none): q_none}
            else:
                arg_q_dict = parent.get_arg_q_dict()
        else:
            # This should never happen...
            raise ValueError("get_single_parent_q called when single parent not ready!!!")
        return arg_q_dict

    def get_single_parent_contigs(self):
        parent = self.get_single_parent()
        if parent.is_ready():
            if parent.count == 0:
                contigs = set()
            else:
                contigs = parent.get_contigs()
        else:
            # This should never happen...
            raise ValueError("get_single_parent_contigs called when single parent not ready!!!")
        return contigs

    def _get_kwargs_for_parent_annotation_kwargs(self, **kwargs) -> Dict:
        """ Use this to pass messages up through to parents """
        return {}

    def _get_annotation_kwargs_for_node(self, **kwargs) -> Dict:
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

    def get_annotation_kwargs(self, **kwargs) -> Dict:
        """ Passed to Variant QuerySet annotate()
            Can be used w/FilteredRelation to force a join to a partition, in which case you need to use
            the alias given in annotate. @see https://github.com/SACGF/variantgrid/wiki/Data-Partitioning """

        kwargs.update(self._get_kwargs_for_parent_annotation_kwargs(**kwargs))
        # Only apply parent annotation kwargs if you actually use their queryset
        a_kwargs = {}
        if self.has_input() and self.uses_parent_queryset:
            if self._has_common_variants():
                kwargs["common_variants"] = True

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

    def get_arg_q_dict(self, disable_cache=False) -> Dict[Optional[str], Dict[str, Q]]:
        """ A Django Q object representing the Variant filters for this node.
            This is the method to override in subclasses - not get_queryset()

            @see https://github.com/SACGF/variantgrid/wiki/Analysis-Nodes#node-q-objects
        """
        # We need this for node counts, and doing a grid query (each page) - and it can take a few secs to generate
        # for some nodes (Comp HET / pheno) so cache it
        cache_key = self._get_cache_key() + f"q_cache={disable_cache}"
        arg_q_dict: Dict[Optional[str], Dict[str, Q]] = {}
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
                except:
                    log_traceback()
        return arg_q_dict

    def get_contigs(self) -> Set[Contig]:
        """ A set of contigs that contain variants for the node """

        cache_key = self._get_cache_key() + "_contigs"
        contigs: Set[Contig] = cache.get(cache_key)

        if contigs is None:
            if self.has_input():
                contigs = self.get_parent_contigs()
                if self.modifies_parents():
                    node_contigs = self._get_node_contigs()
                    if node_contigs is not None:
                        contigs &= node_contigs
            else:
                node_contigs = self._get_node_contigs()
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

    def get_parent_contigs(self) -> Set[Contig]:
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
    def node_version(self):
        return NodeVersion.get(self)

    @cached_property
    def node_cache(self) -> Optional['NodeCache']:
        if parent := self.get_unmodified_single_parent_node():
            return parent.node_cache
        return NodeCache.objects.filter(node_version=self.node_version,
                                        variant_collection__status=ProcessingStatus.SUCCESS).first()

    def _get_node_cache_arg_q_dict(self) -> Dict[Optional[str], Dict[str, Q]]:
        arg_q_dict = {}
        if self.node_cache:
            arg_q_dict = self.node_cache.variant_collection.get_arg_q_dict()
        return arg_q_dict

    def _get_node_q(self) -> Optional[Q]:
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

    def _get_node_arg_q_dict(self) -> Dict[Optional[str], Dict[str, Q]]:
        """ By default - we assume node implements _get_node_q and none of the filters apply to annotations """
        node_arg_q_dict = {}
        if node_q := self._get_node_q():
            node_arg_q_dict[None] = {self._get_node_q_hash(): node_q}
        return node_arg_q_dict

    def _get_node_contigs(self) -> Optional[Set[Contig]]:
        """ Return the contigs we filter for in this node. None means we don't know how to describe that """
        return None

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

        if a_kwargs:
            # If we apply the kwargs at the same time, it can join to the same table twice.
            # We want to go through and apply each annotation then the filters that use it, so that it forces
            # an inner query. Then do next annotation etc

            for k, v in a_kwargs.items():
                qs = qs.annotate(**{k: v})
                for q in arg_q_dict.pop(k, {}).values():
                    qs = qs.filter(q)

        q_list = []
        # Anything stored under None means filters that don't rely on annotation - do afterwards
        if q_dict := arg_q_dict.pop(None, {}):
            # print(f"q_dict(None): {q_dict}")
            q_list.extend(q_dict.values())

        if arg_q_dict:
            raise Exception(f"arg_q_dict filters {arg_q_dict.keys()} not applied (missing in {a_kwargs=})")

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

    def get_extra_grid_config(self):
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

    @staticmethod
    def get_help_text() -> str:
        raise NotImplementedError()

    @staticmethod
    def get_node_class_label():
        """ Used in create node dropdown """
        raise NotImplementedError()

    def _get_genome_build_errors(self, field_name, field_genome_build: GenomeBuild) -> List:
        """ Used to quickly add errors about genome build mismatches
            This only happens in templates (ran template on sample with different build than hardcoded data)
            In normal analyses, autocomplete restrictions should not allow you to configure data from other builds """
        errors = []
        if field_genome_build != self.analysis.genome_build:
            msg = f"{field_name} genome build: {field_genome_build} different from analysis build: {self.analysis.genome_build}"
            errors.append(msg)
        return errors

    def _get_configuration_errors(self) -> List:
        return []

    def get_parents_and_errors(self):
        """ Returns error array, includes any min/max parent error and node config error """
        if self.has_input():
            return self.get_parent_subclasses_and_errors()
        return [], []

    def _get_analysis_errors(self) -> List[str]:
        if self._cached_analysis_errors is None:
            self._cached_analysis_errors = self.analysis.get_errors()
        return self._cached_analysis_errors

    def get_errors(self, include_parent_errors=True, flat=False):
        """ returns a tuple of (NodeError, str) unless flat=True where it's only string """
        errors = []
        for analysis_error in self._get_analysis_errors():
            errors.append((NodeErrorSource.ANALYSIS, analysis_error))

        _, parent_errors = self.get_parents_and_errors()
        if include_parent_errors:
            errors.extend(parent_errors)
        if self.errors:
            errors.append((NodeErrorSource.INTERNAL_ERROR, self.errors))
        errors.extend((NodeErrorSource.CONFIGURATION, ce) for ce in self._get_configuration_errors())
        if flat:
            errors = AnalysisNode.flatten_errors(errors)
        return errors

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

    def _get_node_extra_columns(self):
        return []

    def _get_inherited_columns(self):
        extra_columns = []
        if self.inherits_parent_columns():
            parent = self.get_single_parent()
            extra_columns.extend(parent.get_extra_columns())
        return extra_columns

    def get_extra_columns(self):
        cache_key = self._get_cache_key() + "_extra_columns"
        extra_columns = cache.get(cache_key)
        if extra_columns is None:
            extra_columns = []
            if self.is_valid():
                extra_columns.extend(self._get_inherited_columns())
            # Only add columns that are unique, as otherwise filters get added twice.
            node_extra_columns = self._get_node_extra_columns()
            for col in node_extra_columns:
                if col not in extra_columns:
                    extra_columns.append(col)
            cache.set(cache_key, extra_columns)
        return extra_columns

    def _get_node_extra_colmodel_overrides(self):
        """ Subclasses should override to add colmodel overrides for JQGrid """
        return {}

    def _get_inherited_colmodel_overrides(self):
        extra_overrides = {}
        if self.inherits_parent_columns():
            parent = self.get_single_parent()
            extra_overrides.update(parent.get_extra_colmodel_overrides())
        return extra_overrides

    def get_extra_colmodel_overrides(self):
        """ For JQGrid - subclasses should override _get_node_extra_colmodel_overrides """

        extra_overrides = {}
        if self.is_valid() and self.uses_parent_queryset:
            extra_overrides.update(self._get_inherited_colmodel_overrides())
        extra_overrides.update(self._get_node_extra_colmodel_overrides())
        return extra_overrides

    def get_node_classification(self):
        if self.is_source():
            classification = "source"
        else:
            classification = "filter"
        return classification

    def has_input(self):
        return self.max_inputs != 0

    def is_source(self):
        return self.has_input() is False

    def is_valid(self):
        return not self.get_errors()

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

    def get_unmodified_single_parent_node(self) -> Optional['AnalysisNode']:
        """ If a node doesn't modify single parent - can use that in some places to re-use cache """
        if self.is_valid() and self.has_input() and not self.modifies_parents():
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

        node_counts = []
        for label, count in label_counts.items():
            node_counts.append(NodeCount(node_version=self.node_version, label=label, count=count))
        if node_counts:
            NodeCount.objects.bulk_create(node_counts)

        total_count = label_counts[BuiltInFilters.TOTAL]

        # Single parent nodes should always reduce the number of variants - run a check to make sure the
        # query wasn't bad and returned more results than it should have
        parents = list(self.get_non_empty_parents())
        if len(parents) == 1:
            parent = parents[0]
            if parent.count < total_count:
                raise ValueError(f"Single parent node {self}(pk={self.pk}) had count={total_count} > {parent=}(pk={parent.pk}) count={parent.count=}")

        return NodeStatus.READY, total_count

    def _load(self):
        """ Override to do anything interesting """

    def load(self):
        """ load is called after parents are run """
        # logging.debug("node %d (%d) load()", self.id, self.version)
        start = time()
        self._load()  # Do before counts in case it affects anything
        status, count = self.node_counts()
        load_seconds = time() - start
        self.update(status=status, count=count, load_seconds=load_seconds)

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
        self_qs = AnalysisNode.objects.filter(pk=self.pk, version=self.version)
        updated = self_qs.update(**kwargs)
        if not updated:
            raise NodeOutOfDateException()

    def save(self, **kwargs):
        """ To avoid race conditions, don't use save() in a celery task (unless running in scheduling_single_worker)
            instead use update() method above """
        # logging.debug("save: pk=%s kwargs=%s", self.pk, str(kwargs))
        super_save = super().save

        if self.parents_changed or self.ancestor_input_samples_changed:
            self.handle_ancestor_input_samples_changed()

        if self.auto_node_name:
            self.name = self.get_node_name()

        # TODO: This causes lots of DB queries... should we change this?
        self.valid = self.is_valid()
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

            super_save(**kwargs)
            if self.update_children:
                # We also need to bump if node has it's own sample - as in templates, we set fields in toposort order
                # So we could go from having multiple proband samples to only one later (thus can set descendants)
                for kid in self.children.select_subclasses():
                    kid.ancestor_input_samples_changed = self.is_source() or self.ancestor_input_samples_changed or \
                                                         self.get_samples_from_node_only_not_ancestors()
                    kid.appearance_dirty = False
                    kid.queryset_dirty = True
                    kid.save()  # Will bump versions
        else:
            super_save(**kwargs)

        # Make sure this always exists
        NodeVersion.objects.get_or_create(node=self, version=self.version)
        # Modify our analyses last updated time
        Analysis.objects.filter(pk=self.analysis.pk).update(modified=timezone.now())

    def set_node_task_and_status(self, celery_task, status):
        cursor = connection.cursor()
        db_pid = cursor.db.connection.get_backend_pid()
        self.update(status=status)

        NodeTask.objects.filter(node=self, version=self.version).update(celery_task=celery_task, db_pid=db_pid)

    def adjust_cloned_parents(self, old_new_map):
        """ If you need to do something with old/new parents """

    def save_clone(self):
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


class AnalysisEdge(edge_factory(AnalysisNode, concrete=False)):
    pass


class NodeTask(TimeStampedModel):
    """ Used to track/lock celery update tasks for nodes (uses DB constraints to ensure 1 per node/version) """

    node = models.ForeignKey(AnalysisNode, on_delete=CASCADE)
    version = models.IntegerField(null=False)
    analysis_update_uuid = models.UUIDField()
    celery_task = models.CharField(max_length=36, null=True)
    db_pid = models.IntegerField(null=True)

    class Meta:
        unique_together = ("node", "version")

    def __str__(self):
        return f"NodeTask: {self.analysis_update_uuid} - {self.node.pk}/{self.version}"


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

    def liftover_complete(self, genome_build: GenomeBuild):
        report_event('Completed AnalysisNode liftover',
                     extra_data={'node_id': self.node_id, 'allele_count': self.get_allele_qs().count()})


class NodeVersion(models.Model):
    """ This will be deleted once a node updates, so make all version specific caches cascade delete from this """
    node = models.ForeignKey(AnalysisNode, on_delete=CASCADE)
    version = models.IntegerField(null=False)

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
    node_version = models.OneToOneField(NodeVersion, on_delete=CASCADE)
    variant_collection = models.OneToOneField(VariantCollection, on_delete=CASCADE)

    @staticmethod
    def get_or_create_for_node(node: AnalysisNode) -> Tuple['NodeCache', bool]:
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


class NodeCount(models.Model):
    node_version = models.ForeignKey(NodeVersion, on_delete=CASCADE)
    label = models.CharField(max_length=100)
    count = models.IntegerField(null=False)

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
            extra_filters_q = get_extra_filters_q(node.analysis.user, node.analysis.genome_build, extra_filters)
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


class NodeVCFFilter(models.Model):
    """ If these exist, they mean use that filter """
    node = models.ForeignKey(AnalysisNode, on_delete=CASCADE)
    vcf_filter = models.ForeignKey(VCFFilter, on_delete=CASCADE, null=True)  # null = 'PASS'

    @staticmethod
    def filter_for_node(node, vcf):
        """ returns vfc but also where vcf_filter is NULL (for pass) """
        q_vcf_filter = Q(vcf_filter__isnull=True) | Q(vcf_filter__vcf=vcf)
        return NodeVCFFilter.objects.filter(q_vcf_filter, node=node)


class NodeAlleleFrequencyFilter(models.Model):
    """ Used for various nodes """
    node = models.OneToOneField(AnalysisNode, on_delete=CASCADE)
    group_operation = models.CharField(max_length=1, choices=GroupOperation.choices, default=GroupOperation.ANY)

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
    def get_sample_arg_q_dict(node: AnalysisNode, sample: Sample) -> Dict[Optional[str], Dict[str, Q]]:
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


class NodeAlleleFrequencyRange(models.Model):
    MIN_VALUE = 0
    MAX_VALUE = 1

    filter = models.ForeignKey(NodeAlleleFrequencyFilter, on_delete=CASCADE)
    min = models.FloatField(null=False)
    max = models.FloatField(null=False)

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
