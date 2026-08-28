""" How a node draws itself on the analysis canvas.

    Everything the card shows comes from here rather than per-class CSS selectors - see
    node_utils.get_rendering_dict() for what gets serialized and analysis_nodes.js for the DOM.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class NodeIcon:
    """ Set exactly one of fa (FontAwesome classes) or symbol (an id in node_icon_sprite.html) """
    fa: Optional[str] = None
    symbol: Optional[str] = None


@dataclass(frozen=True)
class NodeChip:
    """ Pill under the node name saying what the node is actually reading - the details that
        change between two nodes of the same class """
    text: str
    icon: Optional[str] = None  # FontAwesome classes
    title: Optional[str] = None  # hover text
    css_class: Optional[str] = None


def significance_chips(selected: list, field_count: int, short_labels: dict, long_labels: dict,
                       css_class_func) -> list[NodeChip]:
    """ Chips for a row of significance pills - an all-on row isn't filtering, so it says nothing """
    if len(selected) == field_count:
        return []
    return [NodeChip(text=short_labels[value], title=long_labels[value], css_class=css_class_func(value))
            for value in selected]
