from django.contrib.auth.models import User
from django.http import Http404
from django.shortcuts import get_object_or_404

from analysis.exceptions import NodeNotFoundException, NodeOutOfDateException
from analysis.models import Analysis, AnalysisNode


def get_analysis_or_404(user: User, analysis_id, write=False) -> Analysis:
    """ Ensures user has view permission on analysis """
    analysis = get_object_or_404(Analysis, pk=analysis_id)
    if write:
        analysis.check_can_write(user)
    else:
        analysis.check_can_view(user)
    return analysis


def get_node_subclass_or_404(user: User, node_id, write=False, version=None) -> AnalysisNode:
    """ Ensures user has view permission on analysis
        If version is passed in, node is loaded, and if version in DB is after version
        throws NodeOutOfDateException """
    try:
        node = AnalysisNode.objects.get_subclass(pk=node_id)
        if version is not None:
            if node.version > version:
                raise NodeOutOfDateException(f"{node.pk}: requested: {version}, DB: {node.version}")
        if write:
            node.analysis.check_can_write(user)
        else:
            node.analysis.check_can_view(user)
        return node
    except AnalysisNode.DoesNotExist:
        msg = f"No AnalysisNode of pk={node_id}"
        raise Http404(msg)


def get_node_subclass_or_non_fatal_exception(user: User, node_id, write=False, **kwargs):
    try:
        node = AnalysisNode.objects.get_subclass(pk=node_id, **kwargs)
        if write:
            node.analysis.check_can_write(user)
        else:
            node.analysis.check_can_view(user)
        return node
    except AnalysisNode.DoesNotExist:  # @UndefinedVariable
        raise NodeNotFoundException(node_id)
