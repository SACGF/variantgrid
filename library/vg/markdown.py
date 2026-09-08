"""
Shared output shapes for `vg`: a map is a list of MapTable, rendered as Markdown by default or JSON
with --json. Rendering is deterministic so committed maps diff cleanly and `--check` can compare bytes.
"""
import json
from dataclasses import asdict, dataclass, field


@dataclass
class MapTable:
    title: str
    columns: list[str]
    rows: list[list[str]] = field(default_factory=list)
    anchor: str = ""
    note: str = ""

    def __post_init__(self):
        if not self.anchor:
            self.anchor = self.title.lower().replace(" ", "-")


def _cell(value) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def render_markdown(heading: str, intro: str, tables: list[MapTable]) -> str:
    lines = [f"# {heading}", "", intro, ""]
    for table in tables:
        lines.append(f"## {table.title}")
        if table.anchor != table.title.lower().replace(" ", "-"):
            lines[-1] += f" <a id=\"{table.anchor}\"></a>"
        lines.append("")
        if table.note:
            lines += [table.note, ""]
        if not table.rows:
            lines += ["(none)", ""]
            continue
        lines.append("| " + " | ".join(table.columns) + " |")
        lines.append("|" + "|".join("---" for _ in table.columns) + "|")
        for row in table.rows:
            lines.append("| " + " | ".join(_cell(c) for c in row) + " |")
        lines.append("")
    return "\n".join(lines)


def render_json(tables: list[MapTable]) -> str:
    return json.dumps([asdict(t) for t in tables], indent=1, sort_keys=True) + "\n"
