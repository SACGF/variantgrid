""" AnalysisNode is the base class that all analysis nodes inherit from. """
import logging
from collections import defaultdict
from datetime import timedelta
from random import random
from time import perf_counter, time
from typing import Optional

from auditlog.context import disable_auditlog
from auditlog.models import AuditlogHistoryField, LogEntry
from auditlog.registry import auditlog
from celery.canvas import Signature
from django.conf import settings
from django.db import connection, models, transaction
from django.db.models import QuerySet, Subquery
from django.db.models.deletion import CASCADE, SET_NULL
from django.utils import timezone
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager

from analysis.exceptions import (
    NodeNotFoundException,
    NodeOutOfDateException,
    NodeParentNotReadyException,
    NonFatalNodeError,
)
from analysis.models.enums import NodeErrorSource, NodeStatus
from analysis.models.models_analysis import Analysis
from analysis.models.nodes.node_allele_frequency import NodeAlleleFrequencyFilter
from analysis.models.nodes.node_audit_log_mixin import NodeAuditLogMixin
from analysis.models.nodes.node_counts import get_node_counts_and_labels_dict
from analysis.models.nodes.node_display_mixin import NodeDisplayMixin
from analysis.models.nodes.node_errors_mixin import NodeErrorsMixin
from analysis.models.nodes.node_query_mixin import NodeQueryMixin
from analysis.models.nodes.node_samples_mixin import NodeSamplesMixin
from analysis.models.nodes.node_vcf_filter import NodeVCFFilter
from analysis.models.nodes.node_version import NodeCache, NodeTask, NodeVersion
from library.constants import MINUTE_SECS
from library.django_utils.django_postgres import get_backend_pid
from library.django_utils.model_utils import get_model_content_type_dict
from snpdb.models import AlleleSource, BuiltInFilters, Variant, Wiki
from snpdb.views.datatable_view import RichColumn

# How long a node's lease is good for. The window is (re)started when a worker claims the node for
# loading, so it measures actual load time rather than how long the task sat in the queue.
LEASE_SECONDS = MINUTE_SECS * 10


def _phase_seconds(phase_start: float) -> float:
    """ How long a load phase took @see NodeVersion.load_data["timings"] """
    return round(perf_counter() - phase_start, 3)


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


class AnalysisNode(NodeSamplesMixin, NodeQueryMixin, NodeDisplayMixin, NodeErrorsMixin, NodeAuditLogMixin,  # pylint: disable=abstract-method
                   TimeStampedModel):
    model = Variant
    objects = NodeInheritanceManager()
    history = AuditlogHistoryField()
    analysis = models.ForeignKey(Analysis, on_delete=CASCADE)
    children = models.ManyToManyField('self', blank=True, symmetrical=False, through='AnalysisEdge',
                                      related_name='_parents')
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
        self._saved_by_parent = False  # Cascade save - the parent's save holds the analysis lock
        self._cached_parents = None
        self._cached_analysis_errors = None
        self._cache_node_q = settings.ANALYSIS_NODE_CACHE_Q  # Disable for unit tests

    def __lt__(self, other):
        return self.pk < other.pk

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
                # A missing count should only ever happen if original bumped version since we were
                # loaded, otherwise the cascade should have set cloned_from to NULL
                if (cloned_count := self.cloned_from.counts.get(label)) is not None:
                    return cloned_count

            if self.has_input():
                parent_non_zero_label_counts = []
                for parent in self.get_non_empty_parents():
                    if parent.count != 0:  # count=0 has 0 for all labels
                        parent_count = parent.node_version.counts.get(label)
                        if parent_count is None:
                            return None  # Parent loaded before this label was configured - run the SQL
                        if parent_count != 0:
                            parent_non_zero_label_counts.append(parent_count)

                if not parent_non_zero_label_counts:
                    # logging.info("all parents had 0 %s counts", label)
                    return 0

                if not self.modifies_parents():
                    if len(parent_non_zero_label_counts) == 1:
                        # logging.info("Single parent, no modification, using that")
                        return parent_non_zero_label_counts[0]
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

    def node_counts(self, timings: Optional[dict] = None):
        """ This is inside Celery task.
            timings - {phase: seconds} the load fills in as it goes, stored on the NodeVersion """

        if timings is None:
            timings = {}

        self.count = None
        # Record provenance first, so the checks below know whether this node's data can move under it
        phase_start = perf_counter()
        live_data_sources = self.get_live_data_sources()
        NodeVersion.objects.filter(pk=self.node_version.pk).update(live_data_sources=live_data_sources)
        self.node_version.live_data_sources = live_data_sources
        timings["live_data_sources"] = _phase_seconds(phase_start)

        phase_start = perf_counter()
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

        timings["counts"] = _phase_seconds(phase_start)

        total_count = label_counts[BuiltInFilters.TOTAL]
        phase_start = perf_counter()
        variant_ids = self._get_variant_ids_to_store(total_count)
        timings["variant_ids"] = _phase_seconds(phase_start)
        if variant_ids is not None:
            # For a small node the PK list is the truth - it and the count came from the same load
            total_count = len(variant_ids)

        # Every label count is a subset of the total - a bigger one means the query fanned out over a
        # multi-valued join, or a cached count is out of sync with the live query
        bigger_than_total = {l: c for l, c in label_counts.items() if c > total_count}

        label_counts[BuiltInFilters.TOTAL] = total_count
        load_data = {"counts": label_counts}
        phase_start = perf_counter()
        load_data.update(self._get_load_data())
        timings["load_data"] = _phase_seconds(phase_start)
        load_data["timings"] = timings
        # Counts are a cache of a query against an immutable node_version, so a re-load (eg after a
        # backoff retry that failed once these were already written) can safely overwrite them.
        # "modified" is the client's signal that counts landed @see nodes_status
        NodeVersion.objects.filter(pk=self.node_version.pk).update(load_data=load_data, variant_ids=variant_ids,
                                                                   modified=timezone.now())
        self.node_version.load_data = load_data
        self.node_version.variant_ids = variant_ids

        if bigger_than_total:
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

    def _get_load_data(self) -> dict:
        """ Extend (via super) to snapshot anything else the node worked out at load - merged into
            NodeVersion.load_data alongside "counts" @see node_counts.
            "warnings" are drawn on the node card, so they're visible without opening the editor (#347) """
        return {"warnings": self.get_warnings()}

    def _load(self):
        """ Override to do anything interesting.
            Return a dict of node fields to write as part of the load - setting them on self isn't
            enough, as load() persists via update() (see #431 - no save() in celery tasks) """

    def load(self):
        """ load is called after parents are run """
        # logging.debug("node %d (%d) load()", self.id, self.version)
        start = time()
        timings = {}
        phase_start = perf_counter()
        load_update_kwargs = self._load() or {}  # Do before counts in case it affects anything
        timings["load"] = _phase_seconds(phase_start)
        status, count = self.node_counts(timings=timings)
        load_seconds = time() - start
        # load_seconds is one number - the phase breakdown is what says where a slow load went
        if slow_seconds := settings.ANALYSIS_NODE_SLOW_LOAD_SECONDS:
            if load_seconds > slow_seconds:
                logging.warning("Node %d.%d slow load %.1fs: %s", self.pk, self.version, load_seconds, timings)
        self.update(status=status, count=count, load_seconds=load_seconds, **load_update_kwargs)

    def add_parent(self, parent, *args, **kwargs):
        if not parent.visible:
            raise NonFatalNodeError("Not connecting children to invisible nodes!")

        existing_connect = AnalysisEdge.objects.filter(parent=parent, child=self)
        if not existing_connect.exists():
            AnalysisEdge.objects.create(parent=parent, child=self)
            self.parents_changed = True
            self._cached_parents = None
        else:
            logging.error("Node(pk=%d).add_parent(pk=%d) already exists!", self.pk, parent.pk)

    def add_child(self, child):
        """ Only creates the edge - unlike add_parent, leaves the child's parents_changed alone """
        AnalysisEdge.objects.create(parent=self, child=child)

    def remove_parent(self, parent):
        """ disconnects parent by deleting edge """
        # Ok to have multiple, just delete first
        edge = AnalysisEdge.objects.filter(parent=parent, child=self).first()
        if edge:  # could be some kind of race condition?
            edge.delete()
        self.parents_changed = True
        self._cached_parents = None

    def parents(self):
        """ Base AnalysisNode rows, so parents of every node type are found """
        return AnalysisNode.objects.filter(children=self)

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
            if not self._saved_by_parent:
                # Lock the analysis row first so concurrent node-tree saves serialize in a consistent order.
                # Saving a node cascades version bumps to its children (see below), and each NodeVersion insert
                # takes a FOR KEY SHARE lock on the referenced analysisnode row. Without this, two concurrent
                # saves over overlapping subtrees grab those row locks in opposite orders and deadlock.
                Analysis.objects.select_related(None).select_for_update().get(pk=self.analysis_id)
            result = self._save(*args, **kwargs)
            if not self._saved_by_parent:
                Analysis.objects.filter(pk=self.analysis_id).update(modified=timezone.now())
            return result

    def _validity_inputs_changed(self) -> bool:
        """ valid is stored at save time - recomputing it (get_errors) costs queries for the node, its
            parents and the analysis, so only do so when something it is computed from has changed.
            Edges made behind the node's back (add_child, cloning) leave these unset: set valid yourself """
        return self.pk is None or self._state.adding or self.parents_changed or self.queryset_dirty or \
            self.ancestor_input_samples_changed

    def _save(self, *args, **kwargs):
        super_save = super().save

        if self.parents_changed or self.ancestor_input_samples_changed:
            self.handle_ancestor_input_samples_changed()

        if self.auto_node_name:
            self.name = self.get_node_name()

        if self._validity_inputs_changed():
            self.valid = self.is_valid
        elif self.errors:
            self.valid = False

        if self.appearance_dirty:
            self.appearance_version += 1

        if self.parents_changed or self.queryset_dirty:
            self.bump_version()

            super_save(*args, **kwargs)
            if self.update_children:
                self._save_children()
        else:
            super_save(*args, **kwargs)

        # Make sure this always exists
        NodeVersion.objects.bulk_create([NodeVersion(node=self, version=self.version)], ignore_conflicts=True)

    def _save_children(self):
        """ Bump every child's version, which cascades on down. Each child is handed what this save already
            knows - the analysis and its errors, and this node where it is the child's only parent - so it
            doesn't query them again, nor re-walk every ancestor above it to find its samples """
        kids = list(self.children.select_subclasses())
        if not kids:
            return

        # We also need to bump if node has it's own sample - as in templates, we set fields in toposort order
        # So we could go from having multiple proband samples to only one later (thus can set descendants)
        ancestor_input_samples_changed = self.is_source or self.ancestor_input_samples_changed or \
            bool(self.get_samples_from_node_only_not_ancestors())
        analysis_errors = self._get_analysis_errors()
        parent_ids_by_kid_id = defaultdict(set)
        for kid_id, parent_id in AnalysisEdge.objects.filter(child__in=kids).values_list("child_id", "parent_id"):
            parent_ids_by_kid_id[kid_id].add(parent_id)
        # A parent is read through its subclass (samples, status) - a base row would answer for no node type
        is_subclass = type(self) is not AnalysisNode

        for kid in kids:
            kid.analysis = self.analysis
            kid._cached_analysis_errors = analysis_errors
            if is_subclass and parent_ids_by_kid_id[kid.pk] == {self.pk}:
                kid._cached_parents = [self]
            kid._saved_by_parent = True
            kid.ancestor_input_samples_changed = ancestor_input_samples_changed
            kid.appearance_dirty = False
            kid.queryset_dirty = True
            kid.save()  # Will bump versions

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
            # Store cloned_from so we can use original's counts
            copy.cloned_from = original_node_version
            copy.save()
            copy._cached_parents = None  # The caller connects the copy's edges

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


class AnalysisEdge(NodeAuditLogMixin, models.Model):
    parent = models.ForeignKey(AnalysisNode, related_name="analysisnode_child", on_delete=CASCADE)
    child = models.ForeignKey(AnalysisNode, related_name="analysisnode_parent", on_delete=CASCADE)

    def __str__(self):
        return f"{self.parent_id} -> {self.child_id}"

    def save(self, *args, **kwargs):
        dag = AnalysisDag.for_node(self.parent_id)
        if self.parent_id == self.child_id or self.child_id in dag.ancestor_ids(self.parent_id):
            raise NonFatalNodeError(f"Not connecting node {self.parent_id} to {self.child_id} - it would create a cycle")
        super().save(*args, **kwargs)

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


class AnalysisDag:
    """ An analysis's edges loaded in one query, so walks are by pk rather than a parents() query per node """

    def __init__(self, analysis_id):
        """ analysis_id can be a Subquery - see for_node """
        self.parent_ids = defaultdict(set)
        self.child_ids = defaultdict(set)
        edges = AnalysisEdge.objects.filter(child__analysis_id=analysis_id).values_list("parent_id", "child_id")
        for parent_id, child_id in edges:
            self.parent_ids[child_id].add(parent_id)
            self.child_ids[parent_id].add(child_id)

    @classmethod
    def for_node(cls, node_id) -> 'AnalysisDag':
        return cls(Subquery(AnalysisNode.objects.filter(pk=node_id).values("analysis_id")[:1]))

    @staticmethod
    def _walk(node_id, links: dict[int, set[int]]) -> set[int]:
        found = set()
        stack = [node_id]
        while stack:
            for linked_id in links.get(stack.pop(), ()):
                if linked_id not in found:
                    found.add(linked_id)
                    stack.append(linked_id)
        return found

    def ancestor_ids(self, node_id) -> set[int]:
        return self._walk(node_id, self.parent_ids)

    def descendant_ids(self, node_id) -> set[int]:
        return self._walk(node_id, self.child_ids)

    def root_ids(self, node_id) -> set[int]:
        """ Ancestors without parents - empty if node_id has none itself """
        return {pk for pk in self.ancestor_ids(node_id) if not self.parent_ids.get(pk)}


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


auditlog.register(AnalysisEdge)
