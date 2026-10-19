""" NodeDisplayMixin - how an AnalysisNode is drawn: its name, label, icon and chips on the DAG card, CSS
    classes and ids, the method summary and the args its editor and grid pages render with.
    The NodeIcon / NodeChip value types are in node_display.py. """
from analysis.models.enums import AnalysisTemplateType
from analysis.models.nodes.node_display import NodeChip, NodeIcon


class NodeDisplayMixin:
    """ Part of AnalysisNode. Each concrete node class implements get_node_name, _get_method_summary,
        get_help_text, get_node_class_label and get_node_class_icon. Optional hooks:
        - get_node_class_label_short / get_node_strip_label: the class strip on the node card
        - get_node_icon / get_node_chips: config-dependent badge icon and the pills under the name
        - get_css_classes, get_rendering_args, get_connection_data, get_grid_post_data """

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

    def get_node_icon(self) -> NodeIcon:
        """ Badge icon for this node - the class default unless config changes it (eg SampleNode
            draws the patient's pedigree shape) """
        return self.get_node_class_icon()

    def get_node_chips(self) -> list[NodeChip]:
        """ Pills under the node name - what this node is reading, from its saved config """
        return []
