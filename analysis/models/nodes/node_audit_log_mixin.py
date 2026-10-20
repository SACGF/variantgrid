""" NodeAuditLogMixin puts the owning node's analysis/node ids on django-auditlog entries, so a change to
    a node-owned row (an edge, a filter, a node's own fields) shows up in that analysis' history. """
import abc


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
