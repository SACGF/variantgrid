"""
`vg inspect <kind> <key>`: one object's whole graph by domain kind (claude/domain.md nouns), not by table.

Owns the KINDS registry (kind -> inspector module under library/vg/inspect/), the lookup rules (a pk or a
natural key per kind), the depth / list-cap contract every inspector follows and the shared Markdown
renderer. Entry points: `inspect(kind, key, depth)` returns nested dicts (the --json shape in
claude/plans/agent_system.md Appendix D), `render_inspection` prints them.

Every inspector runs inside a read-only transaction that is rolled back, so an inspector can never write.
Lists are capped at CAP items and carry `truncated: n` for the rest; relations are followed to `depth`
and beyond it appear as a one-line summary. Inspector modules are imported on demand so `vg` itself stays
cheap to import.
"""
from importlib import import_module
from typing import Any

from django.db import transaction

CAP = 10
KINDS = ("variant", "allele", "sample", "vcf", "classification", "analysis", "gene", "transcript", "user", "lab")


def inspect(kind: str, key: str, depth: int = 2) -> dict[str, Any]:
    if kind not in KINDS:
        raise LookupError(f"Unknown kind '{kind}'; one of {', '.join(KINDS)}")
    inspector = import_module(f"library.vg.inspect.{kind}")
    with transaction.atomic():
        data = inspector.inspect(key, depth)
        transaction.set_rollback(True)
    return {"kind": kind, **data}


def capped(items, render, cap: int = CAP) -> dict[str, Any]:
    """ {"count": n, "items": [render(x) ...][:cap], "truncated": n - cap} from a queryset or list """
    if hasattr(items, "count") and not isinstance(items, (list, tuple)):
        count = items.count()
        rows = list(items[:cap])
    else:
        rows = list(items)
        count = len(rows)
        rows = rows[:cap]
    result: dict[str, Any] = {"count": count, "items": [render(x) for x in rows]}
    if count > len(rows):
        result["truncated"] = count - len(rows)
    return result


def ref(kind: str, obj, label: str | None = None) -> dict[str, Any]:
    """ A one-line pointer to another inspectable thing: what to run to follow it """
    return {"kind": kind, "id": obj.pk, "label": label or str(obj)}


def render_inspection(data: dict[str, Any]) -> str:
    lines: list[str] = []
    _render_value(lines, data, indent=0)
    return "\n".join(lines)


def _scalar(value) -> bool:
    return value is None or isinstance(value, (str, int, float, bool))


def _inline(value) -> str:
    if isinstance(value, dict) and set(value) >= {"kind", "id"}:
        return f"{value['kind']} {value['id']}" + (f"  {value['label']}" if value.get("label") else "")
    if isinstance(value, dict) and all(_scalar(v) for v in value.values()):
        return ", ".join(f"{k}={_fmt(v)}" for k, v in value.items())
    return _fmt(value)


def _fmt(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.4g}"
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


def _item_line(item: dict[str, Any]) -> str:
    """ A ref-shaped item reads `sample 4  name, zygosity=E`; any other flat dict as k=v pairs """
    rest = {k: v for k, v in item.items() if k not in ("kind", "id", "label")}
    pairs = ", ".join(f"{k}={_inline(v)}" for k, v in rest.items())
    if set(item) >= {"kind", "id"}:
        return _inline({k: item[k] for k in ("kind", "id", "label") if k in item}) + (f", {pairs}" if pairs else "")
    return pairs


def _render_value(lines: list[str], data: dict[str, Any], indent: int):
    pad = "  " * indent
    for key, value in data.items():
        if key in ("kind", "id") and indent == 0:
            continue
        if _scalar(value):
            lines.append(f"{pad}{key}: {_fmt(value)}")
        elif isinstance(value, dict) and "items" in value and "count" in value:
            head = f"{pad}{key} ({value['count']:,})"
            if value.get("truncated"):
                head += f"  showing {len(value['items'])}, {value['truncated']:,} more"
            lines.append(head)
            for item in value["items"]:
                if _scalar(item):
                    lines.append(f"{pad}  - {_fmt(item)}")
                elif isinstance(item, dict) and all(_scalar(v) or (isinstance(v, dict) and "kind" in v) for v in item.values()):
                    lines.append(f"{pad}  - " + _item_line(item))
                else:
                    lines.append(f"{pad}  -")
                    _render_value(lines, item, indent + 2)
        elif isinstance(value, dict) and set(value) >= {"kind", "id"} and len(value) <= 3:
            lines.append(f"{pad}{key}: {_inline(value)}")
        elif isinstance(value, dict):
            if all(_scalar(v) for v in value.values()) and len(value) <= 6:
                lines.append(f"{pad}{key}: {_inline(value)}")
            else:
                lines.append(f"{pad}{key}:")
                _render_value(lines, value, indent + 1)
        elif isinstance(value, (list, tuple)):
            if all(_scalar(v) for v in value):
                lines.append(f"{pad}{key}: " + (", ".join(_fmt(v) for v in value) or "-"))
            else:
                lines.append(f"{pad}{key}:")
                for item in value:
                    if isinstance(item, dict):
                        lines.append(f"{pad}  -")
                        _render_value(lines, item, indent + 2)
                    else:
                        lines.append(f"{pad}  - {_fmt(item)}")
        else:
            lines.append(f"{pad}{key}: {value}")
