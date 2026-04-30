""" Invalidate analysis-node caches when their underlying input data is deleted.

When a Sample / Cohort / Trio / Pedigree / Quad row is deleted, any analysis source
node that pinned it would otherwise keep returning a stale cached arg_q_dict that
references annotations (e.g. cohortgenotype_<id>) that no longer exist - producing a
FieldError at queryset build time.

For each affected node we:
  1. Bump the node's version (cascading queryset_dirty + version bumps to children
     via AnalysisNode.save()), so the q-dict cache misses on next access.
  2. Schedule create_and_launch_analysis_tasks for each affected analysis once the
     deletion transaction commits, so DIRTY nodes actually re-evaluate and land in
     ERROR_CONFIGURATION instead of sitting in DIRTY forever (which leaves the grid
     spinner hanging).
"""
import logging

from django.db import transaction
from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode


def _bump_nodes(q):
    """ Bump version + cascade to children for AnalysisNodes matching `q`.
        Returns the set of analysis_ids touched. """
    analysis_ids = set()
    for node in AnalysisNode.objects.filter(q).select_subclasses():
        try:
            node.queryset_dirty = True
            node.update_children = True
            node.save()
            analysis_ids.add(node.analysis_id)
        except Exception:
            logging.exception("Failed to bump version for node %s during source-data deletion", node.pk)
    return analysis_ids


def _schedule_analysis_updates(analysis_ids):
    if not analysis_ids:
        return
    # pylint: disable=import-outside-toplevel
    from analysis.tasks.analysis_update_tasks import create_and_launch_analysis_tasks

    def _launch():
        for analysis_id in analysis_ids:
            try:
                create_and_launch_analysis_tasks.si(analysis_id).apply_async()
            except Exception:
                logging.exception("Failed to launch analysis update for analysis %s", analysis_id)

    transaction.on_commit(_launch)


def handle_sample_pre_delete(sender, instance, **kwargs):
    # Source SampleNode + AncestorSampleMixin filter nodes (Zygosity/GeneList/AlleleFrequency/MOI)
    # all carry a `sample` FK. Walk them via the multi-table inheritance reverse relations
    # on AnalysisNode in a single query.
    q = (
        Q(samplenode__sample=instance)
        | Q(zygositynode__sample=instance)
        | Q(genelistnode__sample=instance)
        | Q(allelefrequencynode__sample=instance)
        | Q(moinode__sample=instance)
    )
    _schedule_analysis_updates(_bump_nodes(q))


def handle_cohort_pre_delete(sender, instance, **kwargs):
    _schedule_analysis_updates(_bump_nodes(Q(cohortnode__cohort=instance)))


def handle_trio_pre_delete(sender, instance, **kwargs):
    _schedule_analysis_updates(_bump_nodes(Q(trionode__trio=instance)))


def handle_pedigree_pre_delete(sender, instance, **kwargs):
    _schedule_analysis_updates(_bump_nodes(Q(pedigreenode__pedigree=instance)))


def handle_quad_pre_delete(sender, instance, **kwargs):
    _schedule_analysis_updates(_bump_nodes(Q(quadnode__quad=instance)))
