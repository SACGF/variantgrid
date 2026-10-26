""" NodeQueryMixin - how an AnalysisNode becomes a Variant queryset: annotation kwargs and the arg Q dict
    (filters keyed by the annotation alias they need) composed down the DAG, the contigs a node can
    hold, the node cache, and get_queryset which applies them. Module helpers are for nodes building Qs
    (queryset_to_pk_in_q, querysets_to_pk_any_array_q) and for code evaluating a node's queryset
    (node_query_planner_settings). """
import operator
from contextlib import contextmanager
from functools import cached_property, reduce
from typing import Optional

from django.conf import settings
from django.core.cache import cache
from django.core.exceptions import FieldError
from django.db.models import BooleanField, F, Func, IntegerField, QuerySet, Value
from django.db.models.expressions import RawSQL
from django.db.models.query_utils import Q

from analysis.models.nodes.node_version import NodeCache, NodeVersion
from annotation.annotation_version_querysets import get_variant_queryset_for_annotation_version
from library.django_utils.database_utils import queryset_to_sql
from library.django_utils.major_operation import planner_join_collapse_limit
from library.log_utils import log_traceback
from snpdb.enums import ProcessingStatus
from snpdb.models import Contig, VariantCollection
from snpdb.variant_collection import write_sql_to_variant_collection


@contextmanager
def node_query_planner_settings():
    """ Wrap anything that evaluates a node's Variant queryset - grid pages, exports, loads, tag recounts.
        With the grid's columns selected a node query joins 50+ relations, past the planner's default
        join_collapse_limit, and the SQL's join order starts from the whole variant table with the node's
        selective filter applied last (a 421 variant gene search grid took over two minutes).
        @see settings.ANALYSIS_NODE_QUERY_JOIN_COLLAPSE_LIMIT """
    with planner_join_collapse_limit(settings.ANALYSIS_NODE_QUERY_JOIN_COLLAPSE_LIMIT):
        yield


def queryset_to_pk_in_q(qs: QuerySet) -> Q:
    """ Embed a queryset as pk IN (subquery), NOT list(qs.values_list("pk")).
        Callers reach this with querysets that are large by construction - a small node answers
        get_arg_q_dict with the literal PK list it stored at load - so list() here could pull an
        unbounded number of PKs into Python (e.g. a 7.4M-row cohort).
        We render to RawSQL, capturing the compiled SQL + params (incl. the partition table rewrite the
        TransformerQuerySet applies in as_sql), so it runs as a single DB-side semi-join. RawSQL also keeps
        arg_q_dict picklable for the q-cache cache.set: a live TransformerQuerySet embedded in the Q is
        unpicklable (closure-local compiler classes) which previously forced the materialise-to-list
        workaround (issue #546, #240, ad35a7fb1). """
    pk_qs = qs.values_list("pk", flat=True)
    sql, params = pk_qs.query.sql_with_params()
    return Q(pk__in=RawSQL(sql, params))


def querysets_to_pk_any_array_q(querysets: list[QuerySet]) -> Q:
    """ The union of several querysets as pk = ANY(ARRAY(<qs1> UNION <qs2> ...)), rendered to RawSQL as
        queryset_to_pk_in_q does. UNION also dedupes a pk more than one queryset returns.

        An OR of pk IN (subquery) has no index to use, so Postgres walks the whole pkey index checking
        every row against each hashed subplan (#1894). pk IN (<union>) fixes count() but the planner
        misestimates the union's size and hash joins every variant for a sorted page. ARRAY() runs the
        union once as an InitPlan, so the plan drives from the pkey. The pk is a Django column rather
        than SQL text, so it's relabelled when this Q ends up inside a subquery (U0) """
    sql_list = []
    params = []
    for qs in querysets:
        sql, qs_params = qs.values_list("pk", flat=True).query.sql_with_params()
        sql_list.append(f"({sql})")
        params.extend(qs_params)
    union = RawSQL(" UNION ".join(sql_list), params)
    # Not Q(pk=...) - the exact lookup wraps its rhs in parentheses and "= (ANY(...))" isn't valid SQL
    pk_any_array = Func(F("pk"), union, arg_joiner=" = ANY(ARRAY", template="%(expressions)s)",
                        output_field=BooleanField())
    return Q(pk_any_array)


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
                fe.add_note(f"annotation kwarg: {k}: {q=}.")
                raise

    q_list = []
    # Anything stored under None means filters that don't rely on annotation - do afterwards
    if q_dict := arg_q_dict.pop(None, {}):
        q_list.extend(q_dict.values())

    if arg_q_dict:
        raise ValueError(f"arg_q_dict filters {arg_q_dict.keys()} not applied (missing in {a_kwargs=})")
    return qs, q_list


class NodeQueryMixin:
    """ Part of AnalysisNode. A node is a filter - override these, not get_queryset:
        - _get_node_q: the node's own filter, or _get_node_arg_q_dict when it needs an annotation alias
          (with _get_annotation_kwargs_for_node to add that annotation)
        - _get_node_q_hash: lets MergeNode drop a filter two parents share
        - _get_node_contigs: the contigs the node restricts to, if it can say
        - _get_kwargs_for_parent_annotation_kwargs / _has_common_variants: messages passed up to parents
        - _get_arg_q_dict_from_parents_and_node / get_parent_arg_q_dict: a node that combines several parents
        - _get_cache_key: when the Q depends on more than the NodeVersion
        - use_cache / write_cache / _get_node_cache_arg_q_dict: materialise the node into a VariantCollection
        - get_extra_filters_tag_q: a node that scopes tag 'extra_filters' itself """

    def _get_model_queryset(self):
        self.analysis.check_valid()
        return get_variant_queryset_for_annotation_version(self.analysis.annotation_version)

    def get_single_parent_arg_q_dict(self) -> dict[Optional[str], dict[str, Q]]:
        arg_q_dict = {}
        parent = self.get_single_parent()
        if parent.is_ready:
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

    @staticmethod
    def q_all():
        return Q(pk__isnull=False)

    @staticmethod
    def q_none():
        return ~NodeQueryMixin.q_all()

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
        # Issue #546 explicit-PK substitution, in the one place every consumer passes through: a node
        # whose load stored its exact pk list answers with a literal Q(pk__in=[...]) - for its own grid,
        # export and recounts, and for every child composing it - so Postgres plans a bitmap-or over the
        # pk index instead of re-running the filter chain. Checked ahead of the Redis Q cache, which
        # holds the real filter the load itself ran (the pks only exist once that load has finished)
        if (variant_ids := NodeQueryMixin.get_cached_node_pks(self)) is not None:
            q = Q(pk__in=variant_ids)
            return {None: {q: q}}

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
            or None for a large node - or one loaded before the PKs were stored with the count.
            A list over the current setting is treated as absent, so lowering it (or the profiler's
            --pk-substitution off) takes effect without a reload """
        try:
            variant_ids = node.node_version.variant_ids
        except NodeVersion.DoesNotExist:
            return None
        if variant_ids is not None and len(variant_ids) > settings.ANALYSIS_NODE_STORE_ID_SIZE_MAX:
            return None
        return variant_ids

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
        if (node_contigs := self._get_node_contigs()) is None:
            return None
        return node_contigs | {Contig.get_gene_level()}

    def get_queryset(self, extra_filters_q=None, extra_annotation_kwargs=None, arg_q_dict=None,
                     disable_cache=False):
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
            qs = qs.filter(q)

        # Clear ordering, @see
        # https://docs.djangoproject.com/en/3.0/topics/db/aggregation/#interaction-with-default-ordering-or-order-by
        return qs.order_by()
