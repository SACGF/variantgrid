""" AnalysisNode is the base class that all analysis nodes inherit from. """
import logging
import operator
from functools import reduce
from random import random
from time import time
from typing import Tuple, Sequence, List, Dict, Optional

from celery.canvas import Signature
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
from lazy import lazy
from model_utils.managers import InheritanceManager

from analysis.exceptions import NonFatalNodeError, NodeParentErrorsException, NodeConfigurationException, \
    NodeParentNotReadyException
from analysis.models.enums import GroupOperation, NodeStatus, NodeColors, NodeErrorSource, AnalysisTemplateType
from analysis.models.models_analysis import Analysis
from analysis.models.nodes.node_counts import get_extra_filters_q, get_node_counts_and_labels_dict
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from library.database_utils import queryset_to_sql
from library.django_utils import thread_safe_unique_together_get_or_create
from library.log_utils import report_event
from snpdb.models import BuiltInFilters, Sample, Variant, VCFFilter, Wiki, Cohort, VariantCollection, \
    ProcessingStatus, GenomeBuild, AlleleSource
from snpdb.variant_collection import write_sql_to_variant_collection
from classification.models import Classification, post_delete
from variantgrid.celery import app


def _default_position():
    return 10 + random() * 50


class AnalysisNode(node_factory('AnalysisEdge', base_model=TimeStampedModel)):
    model = Variant
    objects = InheritanceManager()
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
    parents_should_cache = models.BooleanField(default=False)  # Node suggests parents use a cache

    PARENT_CAP_NOT_SET = -1
    min_inputs = 1
    max_inputs = 1
    uses_parent_queryset = True
    disabled = False

    # Task Update fields
    analysis_update_uuid = models.UUIDField(null=True, default=None)
    status = models.CharField(max_length=1, choices=NodeStatus.choices, default=NodeStatus.DIRTY)
    celery_task = models.CharField(max_length=36, null=True)
    db_pid = models.IntegerField(null=True)

    UPDATE_TASK = "analysis.tasks.node_update_tasks.update_node_task"
    NODE_CACHE_TASK = "analysis.tasks.node_update_tasks.node_cache_task"
    WAIT_FOR_CACHE_TASK = "analysis.tasks.node_update_tasks.wait_for_cache_task"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.appearance_dirty = False
        self.ancestor_input_samples_changed = False
        self.parents_changed = False
        self.queryset_dirty = False

    def get_subclass(self):
        """ Returns the node loaded as a subclass """
        return AnalysisNode.objects.get_subclass(pk=self.pk)

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
        bams_dict = {}
        for sample in self.get_samples():
            if sample.bam_file_path:
                bams_dict[sample.pk] = sample.bam_file_path
        return bams_dict

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
        qs = AnalysisNode.objects.filter(children=self.id, children__isnull=False)
        parents = list(qs.select_subclasses())
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

    def get_non_empty_parents(self, require_parents_ready=True):
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

    def get_single_parent_q(self):
        parent = self.get_single_parent()
        if parent.is_ready():
            if parent.count == 0:
                q = self.q_none()
            else:
                q = parent.get_q()
        else:
            # This should never happen...
            raise ValueError("get_single_parent_q called when single parent not ready!!!")
        return q

    def _get_annotation_kwargs_for_node(self) -> Dict:
        """ Override this method per-node.
            Any key/values in here MUST be consistent - as annotation_kwargs from multiple
            nodes may be combined in the MergeNode
        """
        annotation_kwargs = {}
        if self.node_cache:
            annotation_kwargs.update(self.node_cache.variant_collection.get_annotation_kwargs())
        return annotation_kwargs

    def get_annotation_kwargs(self) -> Dict:
        """ Passed to Variant QuerySet annotate()
            Can be used w/FilteredRelation to force a join to a partition, in which case you need to use
            the alias given in annotate. @see https://github.com/SACGF/variantgrid/wiki/Data-Partitioning """
        a_kwargs = {}
        # Only apply parent annotation kwargs if you actually use their queryset
        if self.has_input() and self.uses_parent_queryset:
            for parent in self.get_non_empty_parents():
                a_kwargs.update(parent.get_annotation_kwargs())

        a_kwargs.update(self._get_annotation_kwargs_for_node())
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

    def get_q(self, disable_cache=False):
        """ A Django Q object representing the Variant filters for this node.
            This is the method to override in subclasses - not get_queryset() as:

            Chains of filters to a reverse foreign key relationship causes
            Multiple joins, so use Q objects which are combined at the end

            qs = qs.filter(table_1__val=1)
            qs = qs.filter(table_2__val=2)

            This is not necessarily equal to:

            qs.filter(table_1__val=1, table_2__val=2)

            @see https://docs.djangoproject.com/en/2/topics/db/queries/#spanning-multi-valued-relationships
        """

        if disable_cache is False:
            if cache_q := self._get_node_cache_q():
                return cache_q

        if self.has_input():
            q = self.get_parent_q()
            if self.modifies_parents():
                if node_q := self._get_node_q():
                    q &= node_q
        else:
            q = self.q_all()
            if node_q := self._get_node_q():
                q = node_q
        return q

    def get_parent_q(self):
        if self.min_inputs == 1:
            return self.get_single_parent_q()
        raise NotImplementedError("You need to implement a non-default 'get_parent_q' if you have more than 1 parent")

    @property
    def use_cache(self):
        """ At the moment we only cache when a child requests it """
        return AnalysisEdge.objects.filter(parent=self, child__parents_should_cache=True).exists()

    def write_cache(self, variant_collection: VariantCollection):
        qs = self.get_queryset(disable_cache=True)
        qs = qs.annotate(variant_collection_id=Value(variant_collection.pk, output_field=IntegerField()))
        sql = queryset_to_sql(qs.values_list('pk', 'variant_collection_id'))
        write_sql_to_variant_collection(variant_collection, sql)

    @lazy
    def node_version(self):
        return thread_safe_unique_together_get_or_create(NodeVersion, node=self, version=self.version)[0]

    @lazy
    def node_cache(self) -> Optional['NodeCache']:
        if parent := self.get_unmodified_single_parent_node():
            return parent.node_cache
        return NodeCache.objects.filter(node_version=self.node_version,
                                        variant_collection__status=ProcessingStatus.SUCCESS).first()

    def _get_node_cache_q(self) -> Optional[Q]:
        q = None
        if self.node_cache:
            q = self.node_cache.variant_collection.get_q()
        return q

    def _get_node_q(self) -> Optional[Q]:
        raise NotImplementedError()

    def _get_unfiltered_queryset(self, **extra_annotation_kwargs):
        """ Unfiltered means before the get_q() is applied
            extra_annotation_kwargs is applied AFTER node's annotation kwargs
        """
        qs = self._get_model_queryset()
        a_kwargs = self.get_annotation_kwargs()
        a_kwargs.update(extra_annotation_kwargs)
        if a_kwargs:
            # Clear ordering, @see
            # https://docs.djangoproject.com/en/3.0/topics/db/aggregation/#interaction-with-default-ordering-or-order-by
            qs = qs.annotate(**a_kwargs).order_by()
        return qs

    def get_queryset(self, extra_filters_q=None, extra_annotation_kwargs=None,
                     inner_query_distinct=False, disable_cache=False):
        if extra_annotation_kwargs is None:
            extra_annotation_kwargs = {}
        qs = self._get_unfiltered_queryset(**extra_annotation_kwargs)
        q = self.get_q(disable_cache=disable_cache)
        if extra_filters_q:
            q &= extra_filters_q
        filtered_qs = qs.filter(q)

        if self.queryset_requires_distinct:
            if inner_query_distinct:
                qs = qs.filter(pk__in=filtered_qs.values_list("pk", flat=True))
            else:
                qs = filtered_qs.distinct()
        else:
            qs = filtered_qs

        return qs

    def get_extra_grid_config(self):
        return {}

    def get_class_name(self):
        return self.__class__.__name__

    def get_identifier(self):
        return f"{self.get_class_name()}-{self.pk}"

    def get_css_classes(self):
        """ returns list of css classes """
        css_classes = ["window"]
        if self.output_node:
            css_classes.append("output-node")
        if self.analysis.template_type == AnalysisTemplateType.TEMPLATE and self.analysisvariable_set.exists():
            css_classes.append("variable-node")
        return css_classes

    def get_input_count(self):
        parents = self.get_non_empty_parents()
        return sum([p.get_output_count() for p in parents])

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
            for error in errors:
                html_summary += f"<li>{error}</li>"
            html_summary += "</ul>"
        return html_summary

    def get_node_name(self):
        """ Automatic node name """
        raise NotImplementedError(f"Node Class: {self.get_class_name()}")

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

    def get_errors(self, include_parent_errors=True, flat=False):
        """ returns a tuple of (NodeError, str) unless flat=True where it's only string """
        errors = []
        try:
            self.analysis.check_valid()
        except ValueError as ve:
            errors.append((NodeErrorSource.ANALYSIS, str(ve)))
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
        cache_key = f"{self.pk}_{self.version}_extra_columns"
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
        if self.version:
            DELETE_CACHE_TASK = "analysis.tasks.node_update_tasks.delete_old_node_versions"
            app.send_task(DELETE_CACHE_TASK, args=(self.pk, self.version))

        self.version += 1
        self.status = NodeStatus.DIRTY
        self.count = None
        self.errors = None

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

    def get_cached_label_count(self, label):
        """ Override for optimisation """

        try:
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
            label_count = self.get_cached_label_count(label)
            if label_count is not None:
                label_counts[label] = label_count

        counts_to_get -= set(label_counts)

        logging.debug("%s cached counts: %s", self, label_counts)
        if counts_to_get:
            logging.debug("%s needs DB request for %s", self, counts_to_get)
            retrieved_label_counts = get_node_counts_and_labels_dict(self)
            label_counts.update(retrieved_label_counts)

        node_version = NodeVersion.get(self)
        for label, count in label_counts.items():
            NodeCount.objects.create(node_version=node_version, label=label, count=count)

        return NodeStatus.READY, label_counts[BuiltInFilters.TOTAL]

    def load(self):
        """ load is called after parents are run """
        # logging.debug("node %d (%d) load()", self.id, self.version)
        start = time()
        status, count = self.node_counts()
        logging.debug("node_counts returned (%s, %d)", status, count)

        this = AnalysisNode.objects.filter(pk=self.pk, version=self.version)
        load_seconds = time() - start
        this.update(status=status,
                    count=count,
                    celery_task=None,
                    db_pid=None,
                    load_seconds=load_seconds)

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

    def save(self, **kwargs):
        # print("save: pk=%s args=%s, kwargs=%s" % (self.pk, str(args), str(kwargs)))
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
            for kid in self.children.select_subclasses():
                # We also need to bump if node has it's own sample - as in templates, we set fields in toposort order
                # So we could go from having multiple proband samples to only one later (thus can set descendants)
                kid.ancestor_input_samples_changed = self.is_source() or self.ancestor_input_samples_changed or \
                                                     self.get_samples_from_node_only_not_ancestors()
                kid.appearance_dirty = False
                kid.queryset_dirty = True
                kid.save()  # Will bump versions
        else:
            super_save(**kwargs)

        # Modfiy our analyses last updated time
        Analysis.objects.filter(pk=self.analysis.pk).update(modified=timezone.now())

    def set_node_task_and_status(self, celery_task, status):
        qs = AnalysisNode.objects.filter(pk=self.pk, version=self.version)
        cursor = connection.cursor()
        db_pid = cursor.db.connection.get_backend_pid()
        qs.update(celery_task=celery_task, status=status, db_pid=db_pid)

    def adjust_cloned_parents(self, old_new_map):
        """ If you need to do something with old/new parents """
        pass

    def save_clone(self):
        node_id = self.pk

        copy = self
        # Have to set both id/pk to None when using model inheritance
        copy.id = None
        copy.pk = None
        copy.version = 1  # 0 is for those being constructed in analysis templates
        copy.save()
        copy_node_version = NodeVersion.get(copy)

        # Copy node counts
        for nc in NodeCount.objects.filter(node_version__node_id=node_id, node_version__version=self.version):
            nc.pk = None
            nc.node_version = copy_node_version
            nc.save()

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
        return NodeVersion.objects.get_or_create(node=node, version=node.version)[0]

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
def post_delete_intersection_cache(sender, instance, *args, **kwargs):
    if instance.variant_collection:
        instance.variant_collection.delete_related_objects()
        instance.variant_collection.delete()


class NodeCount(models.Model):
    node_version = models.ForeignKey(NodeVersion, on_delete=CASCADE)
    label = models.CharField(max_length=100)
    count = models.IntegerField(null=False)

    class Meta:
        unique_together = ("node_version", "label")

    @staticmethod
    def load_for_node(node, label):
        return NodeCount.objects.get(node_version=NodeVersion.get(node), label=label)

    def __str__(self):
        return f"NodeCount({self.node_version}, {self.label}) = {self.count}"


class NodeColumnSummaryCacheCollection(models.Model):
    node_version = models.ForeignKey(NodeVersion, on_delete=CASCADE)
    variant_column = models.TextField(null=False)
    extra_filters = models.TextField(null=False)

    @staticmethod
    def get_counts_for_node(node, variant_column, extra_filters):
        node_version, _ = thread_safe_unique_together_get_or_create(NodeVersion,
                                                                    node=node,
                                                                    version=node.version)

        ncscc, created = NodeColumnSummaryCacheCollection.objects.get_or_create(node_version=node_version,
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

    def get_q(self, allele_frequency_path):
        af_q = None
        try:
            filters = []
            for af_range in self.nodeallelefrequencyrange_set.all():
                # Only apply filter if restricted range.
                # Missing value (historical data) == -1 so those will come through
                and_filters = []
                if af_range.min > 0:
                    and_filters.append(Q(**{allele_frequency_path + "__gte": af_range.min}))
                if af_range.max < 100:
                    and_filters.append(Q(**{allele_frequency_path + "__lte": af_range.max}))

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
    def get_sample_q(node: AnalysisNode, sample: Sample) -> Optional[Q]:
        af_q = None
        if sample:
            try:
                af_q = node.nodeallelefrequencyfilter.get_q(sample.get_cohort_genotype_field("allele_frequency"))
            except NodeAlleleFrequencyFilter.DoesNotExist:
                pass
        return af_q

    def get_description(self):
        # TODO: do this properly with group operators etc
        af_ranges = list(self.nodeallelefrequencyrange_set.all())
        if len(af_ranges) == 1:
            description = str(af_ranges[0])
        else:
            description = f"{self.get_group_operation_display()} of {len(af_ranges)} filters"
        return description


class NodeAlleleFrequencyRange(models.Model):
    filter = models.ForeignKey(NodeAlleleFrequencyFilter, on_delete=CASCADE)
    min = models.IntegerField(null=False)
    max = models.IntegerField(null=False)

    def __str__(self):
        has_min = self.min is not None and self.min > 0
        has_max = self.max is not None and self.max < 100

        if has_min and has_max:
            return f"{self.min} - {self.max}%"
        if has_min:
            return f">={self.min}%"
        if has_max:
            return f"<={self.max}%"
        return ""


class AnalysisClassification(models.Model):
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE)
    classification = models.ForeignKey(Classification, on_delete=CASCADE)
