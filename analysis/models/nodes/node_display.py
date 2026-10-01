""" How a node draws itself on the analysis canvas.

    Everything the card shows comes from here rather than per-class CSS selectors - see
    node_utils.get_rendering_dict() for what gets serialized and analysis_nodes.js for the DOM.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class NodeIcon:
    """ Set exactly one of fa (FontAwesome classes) or symbol (an id in svg_icon_sprite.html) """
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
    count: Optional[int] = None  # Small bubble after the text, eg the x3 on a "VCF" chip
    row_break: bool = False  # draw this chip and the ones after it on a new line
    # Chips nested the way the relations are - a specimen wrapping its extractions wrapping their VCFs
    children: tuple["NodeChip", ...] = ()


def significance_chips(selected: list, field_count: int, short_labels: dict, long_labels: dict,
                       css_class_func) -> list[NodeChip]:
    """ Chips for a row of significance pills - an all-on row isn't filtering, so it says nothing """
    if len(selected) == field_count:
        return []
    return [NodeChip(text=short_labels[value], title=long_labels[value], css_class=css_class_func(value))
            for value in selected]


def grouped_chips(selected: list, groups: dict[str, list], labels: dict, title_prefix: str,
                  exclude: bool = False, icon: Optional[str] = None) -> list[NodeChip]:
    """ One chip per group with anything selected - a whole group is named for itself, part of one lists
        the members picked. Nothing selected isn't filtering, so says nothing """
    negation = "not " if exclude else ""
    chips = []
    for group_name, members in groups.items():
        chosen = [labels[member] for member in members if member in selected]
        if not chosen:
            continue
        chosen_text = ", ".join(chosen)
        if len(chosen) == len(members):
            text = group_name
            title = f"{title_prefix}: {negation}{group_name}"
            if len(members) > 1:
                title += f" ({chosen_text})"
        else:
            text = chosen_text
            title = f"{title_prefix}: {negation}{chosen_text} ({len(chosen)} of {len(members)} {group_name})"
        chips.append(NodeChip(text=f"{negation}{text}", icon=icon, title=title))
    return chips
