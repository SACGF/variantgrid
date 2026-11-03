""" How a node draws itself on the analysis canvas.

    Everything the card shows comes from here rather than per-class CSS selectors - see
    node_utils.get_rendering_dict() for what gets serialized and analysis_nodes.js for the DOM.
"""

from dataclasses import dataclass
from typing import Optional

# One line of a node name. The card is 64px wide and clamps each line to about 2 rows of 9 characters
# (analysis_nodes.js:renderNodeName); this cap only keeps the hover title and the editor heading readable
NAME_LINE_MAX_LENGTH = 30
CHIP_TEXT_MAX_LENGTH = 20


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


def truncate_name_line(text: str, max_length: int = NAME_LINE_MAX_LENGTH) -> str:
    """ Cut to fit the card - whatever is cut off goes in a chip's title """
    if len(text) > max_length:
        text = text[:max_length - 3] + "…"
    return text


def first_and_more_name_line(names: list[str]) -> str:
    """ Line 1 when several things share it - the first, with "+N more" for the rest, which the chips list """
    if len(names) <= 1:
        return truncate_name_line(names[0]) if names else ""
    more = f" +{len(names) - 1} more"
    return truncate_name_line(names[0], NAME_LINE_MAX_LENGTH - len(more)) + more


def named_chip(name: str, title: Optional[str] = None, icon: Optional[str] = None) -> NodeChip:
    """ A chip for one of several things a node reads (a gene list, a term) - the full name is the hover """
    return NodeChip(text=truncate_name_line(name, CHIP_TEXT_MAX_LENGTH), icon=icon, title=title or name)


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
