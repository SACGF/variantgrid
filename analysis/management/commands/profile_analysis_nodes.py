"""Profile AnalysisNode querysets and write a CSV + EXPLAIN plans bundle.

See https://github.com/SACGF/variantgrid/issues/1494

Usage:
    python3 manage.py profile_analysis_nodes \\
        --analysis 12345 67890 \\
        --sample 555 \\
        --rerun --explain \\
        --out /tmp/prof_$(date +%Y%m%d_%H%M%S)

Output bundle:
    nodes.csv                           single index row per profiled node/pattern
    explain_plans/<source>_<id>.json    full EXPLAIN ANALYZE plans (when --explain)
    meta.txt                            run metadata (host, args, PG version, VAVs)
"""
import csv
import json
import os
import platform
import socket
import sys
import time
import traceback
from datetime import datetime

from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from django.db.models import Q

from analysis.models.nodes.analysis_node import AnalysisNode


CSV_FIELDS = [
    "source",
    "analysis_id",
    "node_id",
    "node_type",
    "node_name",
    "config_summary",
    "parent_input_count",
    "count",
    "cached_load_seconds",
    "rerun_count_seconds",
    "explain_planning_ms",
    "explain_execution_ms",
    "explain_plan_file",
    "sql_truncated",
    "error",
]

SQL_TRUNCATE = 1000


class Command(BaseCommand):
    help = "Profile AnalysisNode querysets per-node and emit a CSV + EXPLAIN plans"

    def add_arguments(self, parser):
        parser.add_argument("--analysis", type=int, nargs="*", default=[],
                            help="Analysis ID(s) — every node in each analysis is profiled")
        parser.add_argument("--sample", type=int, nargs="*", default=[],
                            help="Sample ID(s) — synthetic canonical single-node patterns are run per sample")
        parser.add_argument("--rerun", action="store_true",
                            help="Re-execute qs.count() per node and time it (in addition to cached load_seconds)")
        parser.add_argument("--explain", action="store_true",
                            help="Run EXPLAIN (ANALYZE, BUFFERS, FORMAT JSON) per node and dump plans to JSON files")
        parser.add_argument("--node-type", nargs="*", default=None,
                            help="Restrict to specific subclass names (e.g. SampleNode GeneListNode PopulationNode)")
        parser.add_argument("--limit-per-analysis", type=int, default=None,
                            help="Cap nodes per analysis")
        parser.add_argument("--out", type=str, default=None,
                            help="Output directory (default ./profile_<timestamp>)")

    def handle(self, *args, **options):
        if not options["analysis"] and not options["sample"]:
            raise CommandError("Provide at least one of --analysis or --sample")

        out_dir = options["out"] or os.path.join(
            os.getcwd(), f"profile_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
        os.makedirs(out_dir, exist_ok=True)
        plans_dir = os.path.join(out_dir, "explain_plans")
        if options["explain"]:
            os.makedirs(plans_dir, exist_ok=True)

        csv_path = os.path.join(out_dir, "nodes.csv")
        meta_path = os.path.join(out_dir, "meta.txt")

        rows = []
        node_type_filter = set(options["node_type"]) if options["node_type"] else None

        for analysis_id in options["analysis"]:
            self.stdout.write(f"== Analysis {analysis_id} ==")
            rows.extend(self._profile_analysis(
                analysis_id,
                rerun=options["rerun"],
                explain=options["explain"],
                node_type_filter=node_type_filter,
                limit=options["limit_per_analysis"],
                plans_dir=plans_dir,
            ))

        for sample_id in options["sample"]:
            self.stdout.write(f"== Sample {sample_id} (synthetic) ==")
            rows.extend(self._profile_sample_synthetic(
                sample_id,
                rerun=options["rerun"],
                explain=options["explain"],
                plans_dir=plans_dir,
            ))

        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
            writer.writeheader()
            for row in rows:
                writer.writerow({k: row.get(k, "") for k in CSV_FIELDS})

        self._write_meta(meta_path, options)
        self.stdout.write(self.style.SUCCESS(f"Wrote {len(rows)} rows -> {csv_path}"))
        self.stdout.write(f"Bundle: {out_dir}")

    # ---- Analysis-mode profiling ----------------------------------------

    def _profile_analysis(self, analysis_id, *, rerun, explain, node_type_filter, limit, plans_dir):
        from analysis.models import Analysis
        try:
            analysis = Analysis.objects.get(pk=analysis_id)
        except Analysis.DoesNotExist:
            self.stderr.write(f"Analysis {analysis_id} not found")
            return [{
                "source": "analysis", "analysis_id": analysis_id, "error": "Analysis not found",
            }]

        nodes_qs = AnalysisNode.objects.filter(analysis=analysis).select_subclasses().order_by("pk")
        nodes = list(nodes_qs)
        if limit:
            nodes = nodes[:limit]

        rows = []
        for node in nodes:
            node_type = type(node).__name__
            if node_type_filter and node_type not in node_type_filter:
                continue
            row = self._profile_node(
                node,
                source="analysis",
                analysis_id=analysis_id,
                rerun=rerun,
                explain=explain,
                plans_dir=plans_dir,
            )
            rows.append(row)
            self.stdout.write(self._row_summary(row))
        return rows

    def _profile_node(self, node, *, source, analysis_id, rerun, explain, plans_dir):
        node_type = type(node).__name__
        row = {
            "source": source,
            "analysis_id": analysis_id,
            "node_id": node.pk,
            "node_type": node_type,
            "node_name": getattr(node, "name", "") or "",
            "config_summary": _config_summary(node),
            "count": node.count,
            "cached_load_seconds": node.load_seconds,
        }

        try:
            row["parent_input_count"] = _sum_parent_counts(node)
        except Exception as e:
            row["parent_input_count"] = ""
            row["error"] = f"parent_count: {e}"

        try:
            qs = node.get_queryset()
        except Exception as e:
            row["error"] = (row.get("error") or "") + f" get_queryset: {e}"
            row["sql_truncated"] = ""
            return row

        try:
            sql, params = _qs_sql_with_params(qs)
            row["sql_truncated"] = _truncate(sql, SQL_TRUNCATE)
        except Exception as e:
            sql, params = None, None
            row["sql_truncated"] = ""
            row["error"] = (row.get("error") or "") + f" sql: {e}"

        if rerun:
            t0 = time.perf_counter()
            try:
                qs.count()
                row["rerun_count_seconds"] = round(time.perf_counter() - t0, 4)
            except Exception as e:
                row["error"] = (row.get("error") or "") + f" rerun: {e}"

        if explain and sql is not None:
            plan_file = os.path.join(plans_dir, f"{source}_a{analysis_id}_n{node.pk}.json")
            try:
                planning, execution = _run_explain(sql, params, plan_file)
                row["explain_planning_ms"] = planning
                row["explain_execution_ms"] = execution
                row["explain_plan_file"] = os.path.relpath(plan_file, os.path.dirname(plans_dir))
            except Exception as e:
                row["error"] = (row.get("error") or "") + f" explain: {e}"

        return row

    # ---- Synthetic-mode profiling ---------------------------------------

    def _profile_sample_synthetic(self, sample_id, *, rerun, explain, plans_dir):
        from snpdb.models import Sample, Variant
        from annotation.models import VariantAnnotationVersion, VariantGeneOverlap
        from genes.models import GeneList

        try:
            sample = Sample.objects.get(pk=sample_id)
        except Sample.DoesNotExist:
            self.stderr.write(f"Sample {sample_id} not found")
            return [{
                "source": "synthetic", "node_id": sample_id, "error": "Sample not found",
            }]

        rows = []
        cgc = sample.cohort_genotype_collection
        cohort_alias = cgc.cohortgenotype_alias

        # Pattern 1 - sample any-variant (HET or HOM_ALT) via regex on packed zygosity string
        ann_kwargs_cgc = cgc.get_annotation_kwargs()
        sample_idx = cgc.get_sql_index_for_sample_id(sample.pk)
        regex = f"^.{{{sample_idx - 1}}}[EO]" if sample_idx > 1 else "^[EO]"
        qs_variant_regex = (Variant.objects
                            .annotate(**ann_kwargs_cgc)
                            .filter(Q(**{f"{cohort_alias}__samples_zygosity__regex": regex})))
        rows.append(self._profile_synthetic_pattern(
            "sample_variant_regex", sample_id, qs_variant_regex,
            f"sample={sample_id} cgc={cgc.pk} idx={sample_idx}",
            rerun=rerun, explain=explain, plans_dir=plans_dir))

        # Pattern 2 - sample any-variant via Substr + IN
        ann_kwargs_sample = sample.get_annotation_kwargs()
        zyg_alias = sample.zygosity_alias
        qs_variant_substr = (Variant.objects
                             .annotate(**ann_kwargs_cgc)
                             .annotate(**ann_kwargs_sample)
                             .filter(Q(**{f"{zyg_alias}__in": ["E", "O"]})))
        rows.append(self._profile_synthetic_pattern(
            "sample_variant_substr_in", sample_id, qs_variant_substr,
            f"sample={sample_id} alias={zyg_alias}",
            rerun=rerun, explain=explain, plans_dir=plans_dir))

        # Pattern 3 - gene list via VariantGeneOverlap (issue #1542)
        vav_field_names = {f.name for f in VariantAnnotationVersion._meta.get_fields()}
        if "status" in vav_field_names:
            active_vav_filter = {"status": "ACTIVE"}
        else:
            active_vav_filter = {"active": True}
        vav = (VariantAnnotationVersion.objects
               .filter(genome_build=sample.genome_build, **active_vav_filter)
               .order_by("-pk").first())
        gene_list = (GeneList.objects
                     .filter(import_status="S")
                     .order_by("-pk").first())
        if vav and gene_list:
            gene_ids = list(gene_list.get_genes(vav.gene_annotation_release).values_list("pk", flat=True))
            qs_vgo = Variant.objects.filter(
                pk__in=VariantGeneOverlap.objects
                    .filter(version=vav, gene__in=gene_ids)
                    .values_list("variant_id", flat=True)
            )
            rows.append(self._profile_synthetic_pattern(
                "gene_list_via_vgo", sample_id, qs_vgo,
                f"vav={vav.pk} gene_list={gene_list.pk} ({len(gene_ids)} genes)",
                rerun=rerun, explain=explain, plans_dir=plans_dir))
        else:
            rows.append({
                "source": "synthetic", "node_type": "gene_list_via_vgo",
                "analysis_id": "", "node_id": sample_id,
                "config_summary": f"sample={sample_id}",
                "error": f"missing prerequisite: vav={vav} gene_list={gene_list}",
            })

        # Pattern 4 - population gnomAD AF filter on sample variants (subquery)
        try:
            sample_variants_qs = sample.get_variant_qs()
            qs_af_subq = sample_variants_qs.filter(variantannotation__gnomad_af__lte=0.01)
            rows.append(self._profile_synthetic_pattern(
                "population_af_subquery", sample_id, qs_af_subq,
                f"sample={sample_id} gnomad_af<=0.01",
                rerun=rerun, explain=explain, plans_dir=plans_dir))
        except Exception as e:
            rows.append({
                "source": "synthetic", "node_type": "population_af_subquery",
                "analysis_id": "", "node_id": sample_id,
                "config_summary": f"sample={sample_id}",
                "error": f"build qs: {e}",
            })

        # Pattern 5 - population AF filter via materialized pk__in list (the known-fast path)
        try:
            variant_ids = list(sample.get_variant_qs().values_list("pk", flat=True)[:100_000])
            qs_af_list = Variant.objects.filter(
                pk__in=variant_ids,
                variantannotation__gnomad_af__lte=0.01,
            )
            rows.append(self._profile_synthetic_pattern(
                "population_af_pk_in_list", sample_id, qs_af_list,
                f"sample={sample_id} gnomad_af<=0.01 list_size={len(variant_ids)}",
                rerun=rerun, explain=explain, plans_dir=plans_dir))
        except Exception as e:
            rows.append({
                "source": "synthetic", "node_type": "population_af_pk_in_list",
                "analysis_id": "", "node_id": sample_id,
                "config_summary": f"sample={sample_id}",
                "error": f"build qs: {e}",
            })

        for row in rows:
            self.stdout.write(self._row_summary(row))
        return rows

    def _profile_synthetic_pattern(self, name, sample_id, qs, config, *, rerun, explain, plans_dir):
        row = {
            "source": "synthetic",
            "analysis_id": "",
            "node_id": sample_id,
            "node_type": name,
            "node_name": "",
            "config_summary": config,
        }
        try:
            sql, params = _qs_sql_with_params(qs)
            row["sql_truncated"] = _truncate(sql, SQL_TRUNCATE)
        except Exception as e:
            row["error"] = f"sql: {e}"
            return row

        if rerun:
            t0 = time.perf_counter()
            try:
                row["count"] = qs.count()
                row["rerun_count_seconds"] = round(time.perf_counter() - t0, 4)
            except Exception as e:
                row["error"] = (row.get("error") or "") + f" rerun: {e}"

        if explain:
            plan_file = os.path.join(plans_dir, f"synthetic_s{sample_id}_{name}.json")
            try:
                planning, execution = _run_explain(sql, params, plan_file)
                row["explain_planning_ms"] = planning
                row["explain_execution_ms"] = execution
                row["explain_plan_file"] = os.path.relpath(plan_file, os.path.dirname(plans_dir))
            except Exception as e:
                row["error"] = (row.get("error") or "") + f" explain: {e}"

        return row

    # ---- Helpers --------------------------------------------------------

    @staticmethod
    def _row_summary(row):
        bits = [
            row.get("source", ""),
            f"a{row.get('analysis_id') or '-'}",
            f"n{row.get('node_id') or '-'}",
            row.get("node_type", ""),
        ]
        if row.get("count") not in (None, ""):
            bits.append(f"count={row['count']}")
        if row.get("cached_load_seconds") not in (None, ""):
            bits.append(f"cached={row['cached_load_seconds']}s")
        if row.get("rerun_count_seconds") not in (None, ""):
            bits.append(f"rerun={row['rerun_count_seconds']}s")
        if row.get("explain_execution_ms") not in (None, ""):
            bits.append(f"exec={row['explain_execution_ms']}ms")
        if row.get("error"):
            bits.append(f"ERR={row['error']}")
        return "  ".join(bits)

    @staticmethod
    def _write_meta(path, options):
        with connection.cursor() as cursor:
            cursor.execute("SELECT version()")
            pg_version = cursor.fetchone()[0]
        with open(path, "w") as fh:
            fh.write(f"timestamp: {datetime.now().isoformat()}\n")
            fh.write(f"hostname: {socket.gethostname()}\n")
            fh.write(f"platform: {platform.platform()}\n")
            fh.write(f"python: {sys.version}\n")
            fh.write(f"pg_version: {pg_version}\n")
            fh.write(f"argv: {' '.join(sys.argv)}\n")
            fh.write("options:\n")
            for k, v in sorted(options.items()):
                fh.write(f"  {k}: {v}\n")


def _sum_parent_counts(node):
    parents = list(AnalysisNode.objects.filter(children=node.id, children__isnull=False))
    if not parents:
        return ""
    counts = [p.count for p in parents if p.count is not None]
    if not counts:
        return ""
    return sum(counts)


def _config_summary(node):
    bits = [str(node)] if str(node) else []
    skip = {
        "id", "x", "y", "version", "appearance_version", "auto_node_name",
        "output_node", "hide_node_and_descendants_upon_template_configuration_error",
        "ready", "valid", "visible", "count", "errors", "shadow_color",
        "load_seconds", "cloned_from", "status", "name", "analysis", "modified",
        "created", "node_ptr",
    }
    extras = []
    for field in node._meta.get_fields():
        if not getattr(field, "concrete", False):
            continue
        if field.name in skip:
            continue
        try:
            val = getattr(node, field.attname, None)
        except Exception:
            continue
        if val in (None, "", False, 0):
            continue
        extras.append(f"{field.name}={val}")
    if extras:
        bits.append(" ".join(extras[:8]))
    return " | ".join(b for b in bits if b)[:500]


def _qs_sql_with_params(qs):
    sql, params = qs.query.sql_with_params()
    return sql, params


def _truncate(s, n):
    s = " ".join(str(s).split())
    return s if len(s) <= n else s[: n - 3] + "..."


def _run_explain(sql, params, plan_path):
    """Run EXPLAIN (ANALYZE, BUFFERS, FORMAT JSON) and write the plan to plan_path.

    Returns (planning_ms, execution_ms) extracted from the plan."""
    explain_sql = f"EXPLAIN (ANALYZE, BUFFERS, FORMAT JSON) {sql}"
    with connection.cursor() as cursor:
        cursor.execute("SET statement_timeout = 0")
        cursor.execute(explain_sql, params)
        rows = cursor.fetchall()

    plan = rows[0][0]
    if isinstance(plan, str):
        plan = json.loads(plan)
    if isinstance(plan, list):
        plan_root = plan[0]
    else:
        plan_root = plan

    with open(plan_path, "w") as fh:
        json.dump(plan, fh, indent=2)

    planning = plan_root.get("Planning Time")
    execution = plan_root.get("Execution Time")
    return planning, execution
