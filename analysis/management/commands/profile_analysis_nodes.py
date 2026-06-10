"""Profile AnalysisNode querysets and write a CSV + EXPLAIN plans bundle.

See https://github.com/SACGF/variantgrid/issues/1494

Usage:
    python3 manage.py profile_analysis_nodes \\
        --analysis 12345 67890 \\
        --sample 555 \\
        --rerun --explain \\
        --out /tmp/prof_$(date +%Y%m%d_%H%M%S)

    # A/B the issue #546 explicit-PK substitution on a real analysis:
    python3 manage.py profile_analysis_nodes \\
        --analysis 12345 --rerun --explain --pk-substitution both \\
        --out /tmp/prof_546_$(date +%Y%m%d_%H%M%S)
    # then compare rerun_count_seconds / explain_execution_ms between the
    # pk_substitution=on and pk_substitution=off rows.

Output bundle:
    nodes.csv                           single index row per profiled node/pattern
    explain_plans/<source>_<id>.json    full EXPLAIN ANALYZE plans (when --explain)
    meta.txt                            run metadata (host, args, PG version, VAVs)
"""
import csv
import json
import os
import platform
import random
import socket
import sys
import time
from datetime import datetime

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from django.db.models import Q
from django.db.models.functions import Substr as DjSubstr

from analysis.models import Analysis
from analysis.models.nodes.analysis_node import AnalysisNode
from annotation.models import VariantAnnotationVersion, VariantGeneOverlap
from genes.models import GeneList
from snpdb.models import Cohort, Sample, Trio, Variant, VariantCollection
from snpdb.models.models_enums import ProcessingStatus

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
    "pk_substitution",     # #546 explicit-PK substitution mode this row was profiled under (on/off)
    "build_seconds",       # only set by cohort_exclude_vc_join — one-time pre-compute cost
    "vc_record_count",     # only set by cohort_exclude_vc_join — size of the cached set
    "analyze_seconds",     # only set by cohort_exclude_vc_join_postanalyze — ANALYZE cost
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
        parser.add_argument("--trio", type=int, nargs="*", default=[],
                            help="Trio ID(s) — multi-sample regex vs substr-AND comparison (de novo pattern)")
        parser.add_argument("--cohort", type=int, nargs="*", default=[],
                            help="Cohort ID(s) — multi-sample regex vs substr-AND comparison (random ~half carriers, plus exclude lookahead vs substr-OR)")
        parser.add_argument("--cohort-seed", type=int, default=42,
                            help="Seed for deterministic cohort sample-half selection (default 42)")
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
        parser.add_argument("--pk-substitution", choices=["on", "off", "both"], default="on",
                            help="Issue #546 explicit-PK substitution: on (default — small parents "
                                 "collapse to a literal Q(pk__in=[...]) via "
                                 "ANALYSIS_NODE_STORE_ID_SIZE_MAX), off (threshold forced to 0, "
                                 "parents always run as a full subquery), both (run each profile twice "
                                 "for A/B comparison). Only affects --analysis real nodes; the synthetic "
                                 "patterns build querysets directly and never use the substitution.")

        parser.add_argument("--planner-diagnostic", action="store_true",
                            help="For --cohort runs only: also re-run cohort_exclude_lookahead "
                                 "and cohort_exclude_vc_join under work_mem/random_page_cost "
                                 "tuning variants and (for VC) post-ANALYZE — answers 'is the "
                                 "variant-hash slowdown a planner stats issue or a RAM issue?' "
                                 "without changing cluster-wide settings (#1546).")

        parser.add_argument("--merge-threshold-bumped", type=int, default=10000,
                            help="For --analysis MergeNodes: per-parent survey classifies each parent's count "
                                 "against the current ANALYSIS_NODE_STORE_ID_SIZE_MAX setting and this "
                                 "hypothetical bumped value. Parents in (current, bumped] are the candidates "
                                 "that would shift from path C (subquery form) to path A (literal IN list) "
                                 "under the bump (#546). Default 10000. Survey is read-only — no SQL re-execution.")

    def handle(self, *args, **options):
        if not (options["analysis"] or options["sample"] or options["trio"] or options["cohort"]):
            raise CommandError("Provide at least one of --analysis / --sample / --trio / --cohort")

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

        subst_modes = {"on": ["on"], "off": ["off"], "both": ["on", "off"]}[options["pk_substitution"]]

        # Issue #546 explicit-PK substitution is gated on ANALYSIS_NODE_STORE_ID_SIZE_MAX
        # (read live by AnalysisNode.get_small_parent_arg_q_dict). "on" uses the configured
        # threshold (falling back to 1000 if it's unset/0 so "on" is genuinely enabled); "off"
        # forces it to 0 so every parent runs as a full subquery.
        original_threshold = getattr(settings, "ANALYSIS_NODE_STORE_ID_SIZE_MAX", 1000)
        subst_threshold = {"on": original_threshold or 1000, "off": 0}

        try:
            for subst_mode in subst_modes:
                settings.ANALYSIS_NODE_STORE_ID_SIZE_MAX = subst_threshold[subst_mode]
                self._pk_substitution_mode = subst_mode
                if len(subst_modes) > 1:
                    self.stdout.write(
                        f"##### pk_substitution={subst_mode} "
                        f"(ANALYSIS_NODE_STORE_ID_SIZE_MAX={subst_threshold[subst_mode]}) #####")

                def tag(profiled_rows, subst_mode=subst_mode):
                    return self._tag_subst(profiled_rows, subst_mode)

                for analysis_id in options["analysis"]:
                    self.stdout.write(f"== Analysis {analysis_id} ==")
                    rows.extend(tag(self._profile_analysis(
                        analysis_id,
                        rerun=options["rerun"],
                        explain=options["explain"],
                        node_type_filter=node_type_filter,
                        limit=options["limit_per_analysis"],
                        plans_dir=plans_dir,
                        merge_threshold_bumped=options["merge_threshold_bumped"],
                    )))

                for sample_id in options["sample"]:
                    self.stdout.write(f"== Sample {sample_id} (synthetic) ==")
                    rows.extend(tag(self._profile_sample_synthetic(
                        sample_id,
                        rerun=options["rerun"],
                        explain=options["explain"],
                        plans_dir=plans_dir,
                    )))

                for trio_id in options["trio"]:
                    self.stdout.write(f"== Trio {trio_id} (synthetic) ==")
                    rows.extend(tag(self._profile_trio_synthetic(
                        trio_id,
                        rerun=options["rerun"],
                        explain=options["explain"],
                        plans_dir=plans_dir,
                    )))

                for cohort_id in options["cohort"]:
                    self.stdout.write(f"== Cohort {cohort_id} (synthetic) ==")
                    rows.extend(tag(self._profile_cohort_synthetic(
                        cohort_id,
                        rerun=options["rerun"],
                        explain=options["explain"],
                        plans_dir=plans_dir,
                        seed=options["cohort_seed"],
                        planner_diagnostic=options["planner_diagnostic"],
                    )))
        finally:
            settings.ANALYSIS_NODE_STORE_ID_SIZE_MAX = original_threshold

        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
            writer.writeheader()
            for row in rows:
                writer.writerow({k: row.get(k, "") for k in CSV_FIELDS})

        self._write_meta(meta_path, options)
        self.stdout.write(self.style.SUCCESS(f"Wrote {len(rows)} rows -> {csv_path}"))
        self.stdout.write(f"Bundle: {out_dir}")

    # ---- Analysis-mode profiling ----------------------------------------

    def _profile_analysis(self, analysis_id, *, rerun, explain, node_type_filter, limit, plans_dir,
                          merge_threshold_bumped):
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

        threshold_current = getattr(settings, "ANALYSIS_NODE_STORE_ID_SIZE_MAX", 1000) or 0

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

            if node_type == "MergeNode":
                survey_rows = self._merge_parent_survey(
                    node, analysis_id, threshold_current, merge_threshold_bumped)
                for sr in survey_rows:
                    rows.append(sr)
                    self.stdout.write(self._row_summary(sr))
        return rows

    @staticmethod
    def _merge_parent_survey(merge_node, analysis_id, threshold_current, threshold_bumped):
        """Per-parent band classification for a MergeNode (#546).

        Read-only walk: emits one row per non-empty parent reporting the parent's count and
        which path-band it falls into under the current vs hypothetical bumped threshold.
        Path A (literal IN) is taken when count <= threshold; otherwise path C (subquery form).
        Parents in (threshold_current, threshold_bumped] are the experiment candidates.

        Backwards-compat: only touches stable APIs (parent.count, get_non_empty_parents).
        Doesn't reference get_parent_pks / cache_memoize / the substitution helper, which
        only exist on post-#546 code.
        """
        rows = []
        try:
            parents = list(merge_node.get_non_empty_parents())
        except Exception as e:
            return [{
                "source": "merge_parent",
                "analysis_id": analysis_id,
                "node_id": merge_node.pk,
                "node_type": "MergeNode_parent",
                "node_name": f"merge {merge_node.pk}",
                "error": f"get_non_empty_parents: {e}",
            }]

        for parent in parents:
            count = getattr(parent, "count", None)
            if count is None:
                band = "unknown"
            elif count <= threshold_current:
                band = "A_current"
            elif count <= threshold_bumped:
                band = "A_bumped"
            else:
                band = "C"
            parent_type = type(parent).__name__
            parent_name = (getattr(parent, "name", "") or "").replace("\n", " ")[:40]
            config = (f"parent_id={parent.pk} parent_type={parent_type} count={count} "
                      f"threshold_current={threshold_current} threshold_bumped={threshold_bumped} "
                      f"band={band}")
            rows.append({
                "source": "merge_parent",
                "analysis_id": analysis_id,
                "node_id": merge_node.pk,
                "node_type": "MergeNode_parent",
                "node_name": f"merge {merge_node.pk} <- parent {parent.pk} {parent_type} {parent_name}".strip(),
                "config_summary": config,
                "parent_input_count": count if count is not None else "",
                "count": count if count is not None else "",
            })
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
            "pk_substitution": getattr(self, "_pk_substitution_mode", ""),
        }

        try:
            row["parent_input_count"] = _sum_parent_counts(node)
        except Exception as e:
            row["parent_input_count"] = ""
            row["error"] = f"parent_count: {e}"

        # Regenerate the arg_q_dict fresh rather than reading a (possibly production-warmed)
        # cached one keyed only on node version - otherwise toggling the #546 threshold between
        # passes wouldn't change the emitted SQL.
        node._cache_node_q = False

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
            subst_mode = getattr(self, "_pk_substitution_mode", "on")
            plan_file = os.path.join(
                plans_dir, f"{source}_a{analysis_id}_n{node.pk}_subst{subst_mode}.json")
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

        # Pattern 4 - population gnomAD AF filter on sample variants (subquery, single VAV partition)
        if vav:
            try:
                sample_variants_qs = sample.get_variant_qs()
                qs_af_subq = sample_variants_qs.filter(
                    variantannotation__version=vav,
                    variantannotation__gnomad_af__lte=0.01,
                )
                rows.append(self._profile_synthetic_pattern(
                    "population_af_subquery", sample_id, qs_af_subq,
                    f"sample={sample_id} vav={vav.pk} gnomad_af<=0.01",
                    rerun=rerun, explain=explain, plans_dir=plans_dir))
            except Exception as e:
                rows.append({
                    "source": "synthetic", "node_type": "population_af_subquery",
                    "analysis_id": "", "node_id": sample_id,
                    "config_summary": f"sample={sample_id}",
                    "error": f"build qs: {e}",
                })
        else:
            rows.append({
                "source": "synthetic", "node_type": "population_af_subquery",
                "analysis_id": "", "node_id": sample_id,
                "config_summary": f"sample={sample_id}",
                "error": "missing prerequisite: no active VariantAnnotationVersion",
            })

        # Pattern 5 - population AF filter via materialized pk__in list (single VAV partition)
        if vav:
            try:
                variant_ids = list(sample.get_variant_qs().values_list("pk", flat=True))
                qs_af_list = Variant.objects.filter(
                    pk__in=variant_ids,
                    variantannotation__version=vav,
                    variantannotation__gnomad_af__lte=0.01,
                )
                rows.append(self._profile_synthetic_pattern(
                    "population_af_pk_in_list", sample_id, qs_af_list,
                    f"sample={sample_id} vav={vav.pk} gnomad_af<=0.01 list_size={len(variant_ids)}",
                    rerun=rerun, explain=explain, plans_dir=plans_dir))
            except Exception as e:
                rows.append({
                    "source": "synthetic", "node_type": "population_af_pk_in_list",
                    "analysis_id": "", "node_id": sample_id,
                    "config_summary": f"sample={sample_id}",
                    "error": f"build qs: {e}",
                })
        else:
            rows.append({
                "source": "synthetic", "node_type": "population_af_pk_in_list",
                "analysis_id": "", "node_id": sample_id,
                "config_summary": f"sample={sample_id}",
                "error": "missing prerequisite: no active VariantAnnotationVersion",
            })

        for row in rows:
            self.stdout.write(self._row_summary(row))
        return rows

    def _profile_synthetic_pattern(self, name, entity_id, qs, config, *,
                                   rerun, explain, plans_dir, file_prefix="synthetic_s"):
        row = {
            "source": "synthetic",
            "analysis_id": "",
            "node_id": entity_id,
            "node_type": name,
            "node_name": "",
            "config_summary": config,
            "pk_substitution": getattr(self, "_pk_substitution_mode", ""),
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
            subst_mode = getattr(self, "_pk_substitution_mode", "on")
            plan_file = os.path.join(
                plans_dir, f"{file_prefix}{entity_id}_{name}_subst{subst_mode}.json")
            try:
                planning, execution = _run_explain(sql, params, plan_file)
                row["explain_planning_ms"] = planning
                row["explain_execution_ms"] = execution
                row["explain_plan_file"] = os.path.relpath(plan_file, os.path.dirname(plans_dir))
            except Exception as e:
                row["error"] = (row.get("error") or "") + f" explain: {e}"

        return row

    # ---- Trio synthetic -------------------------------------------------

    def _profile_trio_synthetic(self, trio_id, *, rerun, explain, plans_dir):
        try:
            trio = Trio.objects.get(pk=trio_id)
        except Trio.DoesNotExist:
            self.stderr.write(f"Trio {trio_id} not found")
            return [{"source": "synthetic", "node_type": "trio_denovo",
                     "analysis_id": "", "node_id": trio_id, "error": "Trio not found"}]

        cohort = trio.cohort.get_base_cohort()
        cgc = cohort.cohort_genotype_collection
        ann_kwargs = cgc.get_annotation_kwargs()
        cohort_alias = cgc.cohortgenotype_alias
        zyg_field = f"{cohort_alias}__samples_zygosity"
        n = cgc.num_samples

        # Pull packed indices from the underlying CohortSamples (the trio's CohortSamples
        # belong to its own cohort; the base cohort packing uses the same index).
        proband_idx1 = trio.proband.cohort_genotype_packed_field_index + 1
        mother_idx1 = trio.mother.cohort_genotype_packed_field_index + 1
        father_idx1 = trio.father.cohort_genotype_packed_field_index + 1

        # De novo: proband HET/HOM, parents REF/MISSING (treat '.' and 'R' both as not-carrier)
        # Regex: build a length-n string with '[EO]' at proband, '[R.]' at parents, '.' elsewhere
        regex_parts = ["."] * n
        regex_parts[proband_idx1 - 1] = "[EO]"
        regex_parts[mother_idx1 - 1] = "[R.]"
        regex_parts[father_idx1 - 1] = "[R.]"
        regex = "^" + "".join(regex_parts)

        rows = []
        qs_regex = (Variant.objects
                    .annotate(**ann_kwargs)
                    .filter(Q(**{f"{cohort_alias}__samples_zygosity__regex": regex})))
        rows.append(self._profile_synthetic_pattern(
            "trio_denovo_regex", trio_id, qs_regex,
            f"trio={trio_id} cgc={cgc.pk} proband={proband_idx1} mother={mother_idx1} father={father_idx1}",
            rerun=rerun, explain=explain, plans_dir=plans_dir, file_prefix="synthetic_t"))

        # Substr-AND: 3 independent substring conditions
        qs_substr = (Variant.objects
                     .annotate(**ann_kwargs)
                     .annotate(_proband_z=DjSubstr(zyg_field, proband_idx1, 1))
                     .annotate(_mother_z=DjSubstr(zyg_field, mother_idx1, 1))
                     .annotate(_father_z=DjSubstr(zyg_field, father_idx1, 1))
                     .filter(_proband_z__in=["E", "O"])
                     .filter(_mother_z__in=["R", "."])
                     .filter(_father_z__in=["R", "."]))
        rows.append(self._profile_synthetic_pattern(
            "trio_denovo_substr_and", trio_id, qs_substr,
            f"trio={trio_id} cgc={cgc.pk} proband={proband_idx1} mother={mother_idx1} father={father_idx1}",
            rerun=rerun, explain=explain, plans_dir=plans_dir, file_prefix="synthetic_t"))

        for row in rows:
            self.stdout.write(self._row_summary(row))
        return rows

    # ---- Cohort synthetic -----------------------------------------------

    def _profile_cohort_synthetic(self, cohort_id, *, rerun, explain, plans_dir, seed,
                                   planner_diagnostic=False):
        try:
            cohort = Cohort.objects.get(pk=cohort_id)
        except Cohort.DoesNotExist:
            self.stderr.write(f"Cohort {cohort_id} not found")
            return [{"source": "synthetic", "node_type": "cohort",
                     "analysis_id": "", "node_id": cohort_id, "error": "Cohort not found"}]

        base = cohort.get_base_cohort()
        cgc = base.cohort_genotype_collection
        ann_kwargs = cgc.get_annotation_kwargs()
        cohort_alias = cgc.cohortgenotype_alias
        zyg_field = f"{cohort_alias}__samples_zygosity"
        n = cgc.num_samples
        if n < 2:
            return [{"source": "synthetic", "node_type": "cohort",
                     "analysis_id": "", "node_id": cohort_id,
                     "error": f"cohort num_samples={n} is too small"}]

        rng = random.Random(seed)
        half_count = max(1, n // 2)
        chosen = sorted(rng.sample(range(n), half_count))  # 0-based indices
        chosen_1based = [i + 1 for i in chosen]

        # Pattern A: half-carriers regex — '[EO]' at chosen positions, '.' elsewhere
        regex_parts = ["."] * n
        for i in chosen:
            regex_parts[i] = "[EO]"
        regex_pos = "^" + "".join(regex_parts)

        rows = []
        qs_regex = (Variant.objects
                    .annotate(**ann_kwargs)
                    .filter(Q(**{f"{cohort_alias}__samples_zygosity__regex": regex_pos})))
        rows.append(self._profile_synthetic_pattern(
            "cohort_half_carriers_regex", cohort_id, qs_regex,
            f"cohort={cohort_id} cgc={cgc.pk} n={n} chosen={half_count} (seed={seed})",
            rerun=rerun, explain=explain, plans_dir=plans_dir, file_prefix="synthetic_c"))

        # Pattern B: half-carriers substr-AND — N independent substring IN ('E','O') predicates
        qs_substr_and = Variant.objects.annotate(**ann_kwargs)
        for k, idx1 in enumerate(chosen_1based):
            alias = f"_zc{k}"
            qs_substr_and = qs_substr_and.annotate(**{alias: DjSubstr(zyg_field, idx1, 1)})
            qs_substr_and = qs_substr_and.filter(**{f"{alias}__in": ["E", "O"]})
        rows.append(self._profile_synthetic_pattern(
            "cohort_half_carriers_substr_and", cohort_id, qs_substr_and,
            f"cohort={cohort_id} cgc={cgc.pk} n={n} chosen={half_count} (seed={seed})",
            rerun=rerun, explain=explain, plans_dir=plans_dir, file_prefix="synthetic_c"))

        # Pattern C: exclude lookahead — current production anti-pattern
        # Exclude variants where ALL chosen samples are HET/HOM (i.e. a "rare in cohort" filter)
        # Production rewrite: regex_string = f"^((?!{regex_pos[1:]}))"   (drop our '^' since it re-anchors)
        regex_excl = f"^((?!{''.join(regex_parts)}))"
        qs_excl_regex = (Variant.objects
                         .annotate(**ann_kwargs)
                         .filter(Q(**{f"{cohort_alias}__samples_zygosity__regex": regex_excl})))
        rows.append(self._profile_synthetic_pattern(
            "cohort_exclude_lookahead", cohort_id, qs_excl_regex,
            f"cohort={cohort_id} cgc={cgc.pk} n={n} chosen={half_count} (seed={seed})",
            rerun=rerun, explain=explain, plans_dir=plans_dir, file_prefix="synthetic_c"))

        # Pattern D: exclude via substr-OR — equivalent positive form
        # NOT (substring(.,i_0,1) IN E,O AND substring(.,i_1,1) IN E,O ...)
        #   ==  substring(.,i_0,1) NOT IN E,O OR substring(.,i_1,1) NOT IN E,O ...
        qs_excl_or = Variant.objects.annotate(**ann_kwargs)
        for k, idx1 in enumerate(chosen_1based):
            alias = f"_zc{k}"
            qs_excl_or = qs_excl_or.annotate(**{alias: DjSubstr(zyg_field, idx1, 1)})
        or_q = Q()
        for k in range(len(chosen_1based)):
            or_q |= ~Q(**{f"_zc{k}__in": ["E", "O"]})
        qs_excl_or = qs_excl_or.filter(or_q)
        rows.append(self._profile_synthetic_pattern(
            "cohort_exclude_substr_or", cohort_id, qs_excl_or,
            f"cohort={cohort_id} cgc={cgc.pk} n={n} chosen={half_count} (seed={seed})",
            rerun=rerun, explain=explain, plans_dir=plans_dir, file_prefix="synthetic_c"))

        # Pattern E: exclude via cached VariantCollection JOIN
        # Models the "pre-compute once, join for every subsequent query" strategy: build a
        # transient VC populated by the same regex-excl predicate, then measure ONLY the
        # JOIN-query cost. Build cost is recorded separately as `build_seconds` so it's
        # visible but not counted as the rerun number (it's amortised across queries).
        # Mirrors the design of caching sub-cohort filters for the production EXCLUDE hot
        # path (#1546 alternative to v2 GIN dispatch).
        rows.append(self._profile_cohort_exclude_vc_join(
            cohort_id, cgc, regex_excl, chosen,
            f"cohort={cohort_id} cgc={cgc.pk} n={n} chosen={half_count} (seed={seed})",
            rerun=rerun, explain=explain, plans_dir=plans_dir))

        if planner_diagnostic:
            # Diagnostic: re-run cohort_exclude_lookahead under tuning variants — answers
            # "is the variant-hash bottleneck a planner / RAM issue?" without changing
            # cluster-wide settings. SET runs at session level via the context manager;
            # RESET in __exit__ restores defaults so other cohorts in the same invocation
            # are not affected.
            tuning_variants = [
                ("wm4gb",       {"work_mem":         "'4GB'"}),
                ("rpc11",       {"random_page_cost": "1.1"}),
                ("wm4gb_rpc11", {"work_mem":         "'4GB'", "random_page_cost": "1.1"}),
            ]
            for label_suffix, pg_settings in tuning_variants:
                with _PgSessionSettings(**pg_settings):
                    rows.append(self._profile_synthetic_pattern(
                        f"cohort_exclude_lookahead_{label_suffix}", cohort_id, qs_excl_regex,
                        f"cohort={cohort_id} cgc={cgc.pk} n={n} chosen={half_count} "
                        f"(seed={seed}) settings={'+'.join(f'{k}={v}' for k,v in pg_settings.items())}",
                        rerun=rerun, explain=explain, plans_dir=plans_dir,
                        file_prefix="synthetic_c"))

            # VC diagnostic: build VC once, run ANALYZE on the partition, then profile under
            # baseline + tuning variants. The post-ANALYZE row tells us whether bad stats on
            # the freshly-INSERTed VCR partition are why the planner picked Hash(variants) on
            # the prod EXPLAIN — see #1546 closing comment.
            rows.extend(self._profile_cohort_exclude_vc_join_diagnostic(
                cohort_id, cgc, regex_excl, chosen,
                f"cohort={cohort_id} cgc={cgc.pk} n={n} chosen={half_count} (seed={seed})",
                tuning_variants=tuning_variants,
                rerun=rerun, explain=explain, plans_dir=plans_dir))

        for row in rows:
            self.stdout.write(self._row_summary(row))
        return rows

    def _profile_cohort_exclude_vc_join(self, cohort_id, cgc, regex_excl, chosen, config, *,
                                        rerun, explain, plans_dir):
        """ Build a transient VariantCollection with the regex-excl result set, then profile
        the JOIN-query cost. Cleans up the transient VC after profiling so re-runs stay
        idempotent. """
        partition_table = cgc.get_partition_table()
        common_partition_table = None
        if cgc.common_collection_id:
            common_partition_table = cgc.common_collection.get_partition_table()
        # Match Pattern C's regex; sub_cohort production also runs this against both the
        # uncommon and common partitions, so include both when present.
        partitions = [partition_table] + ([common_partition_table] if common_partition_table else [])

        vc_name = f"synthetic_cohort_{cohort_id}_excl_chosen{len(chosen)}"
        # Idempotency: drop any leftover VC from a previous synthetic run on the same cohort.
        for stale in VariantCollection.objects.filter(name=vc_name):
            try:
                stale.delete_related_objects()
            except Exception:  # pragma: no cover — partition may already be gone
                pass
            stale.delete()

        vc = VariantCollection.objects.create(name=vc_name, status=ProcessingStatus.CREATED)
        vc_partition_table = vc.get_partition_table()

        # Build: INSERT FROM SELECT, mirrors what a production cache-build would do (one
        # SQL statement, no row-by-row Python overhead). Time it.
        union_sql_parts = [
            f'SELECT DISTINCT {vc.pk} AS variant_collection_id, variant_id '
            f'FROM "{p}" '
            f'WHERE collection_id = ANY(%s) AND samples_zygosity ~ %s'
            for p in partitions
        ]
        build_sql = (f'INSERT INTO "{vc_partition_table}" (variant_collection_id, variant_id) '
                     f'{" UNION ".join(union_sql_parts)};')
        coll_ids = [cgc.pk] + ([cgc.common_collection_id] if cgc.common_collection_id else [])
        # Each partition predicate needs its own (collection_ids_array, regex) pair.
        params = []
        for _ in partitions:
            params.extend([coll_ids, regex_excl])

        build_seconds = None
        try:
            t0 = time.perf_counter()
            with connection.cursor() as cur:
                cur.execute(build_sql, params)
                vc_record_count = cur.rowcount
            build_seconds = round(time.perf_counter() - t0, 4)
            vc.count = vc_record_count
            vc.status = ProcessingStatus.SUCCESS
            vc.save()
        except Exception as e:
            row = {
                "source": "synthetic", "analysis_id": "", "node_id": cohort_id,
                "node_type": "cohort_exclude_vc_join", "config_summary": config,
                "error": f"vc build: {e}",
            }
            try:
                vc.delete_related_objects()
            except Exception:  # pragma: no cover
                pass
            vc.delete()
            return row

        # Profile the JOIN query — variant ⋈ variantcollectionrecord ⋈ cohortgenotype.
        # cohortgenotype is still joined because production analysis nodes need its
        # cohortgenotype_alias annotation for downstream filters (filters, AF, etc.); the
        # change is that the EXCLUDE predicate becomes a VC membership instead of a regex.
        ann_kwargs = cgc.get_annotation_kwargs()
        ann_kwargs.update(vc.get_annotation_kwargs())
        cohort_alias = cgc.cohortgenotype_alias
        qs = (Variant.objects
              .annotate(**ann_kwargs)
              .filter(Q(**{f"{cohort_alias}__isnull": False}),
                      Q(**{f"{vc.variant_collection_alias}__isnull": False})))
        row = self._profile_synthetic_pattern(
            "cohort_exclude_vc_join", cohort_id, qs, config,
            rerun=rerun, explain=explain, plans_dir=plans_dir,
            file_prefix="synthetic_c")
        row["build_seconds"] = build_seconds
        row["vc_record_count"] = vc.count

        try:
            vc.delete_related_objects()
        except Exception:  # pragma: no cover
            pass
        vc.delete()
        return row

    def _profile_cohort_exclude_vc_join_diagnostic(self, cohort_id, cgc, regex_excl, chosen,
                                                    config, *, tuning_variants,
                                                    rerun, explain, plans_dir):
        """ Diagnostic counterpart to _profile_cohort_exclude_vc_join — builds the VC ONCE,
        runs ANALYZE on the freshly-INSERTed VCR partition, then profiles the JOIN under
        the baseline session settings plus each tuning variant. Tells us whether the prod
        variant-hash slowdown (#1546 closing comment) is fixed by stats alone, or also
        needs work_mem / random_page_cost tuning. """
        partition_table = cgc.get_partition_table()
        common_partition_table = None
        if cgc.common_collection_id:
            common_partition_table = cgc.common_collection.get_partition_table()
        partitions = [partition_table] + ([common_partition_table] if common_partition_table else [])

        vc_name = f"synthetic_cohort_{cohort_id}_excl_diag_chosen{len(chosen)}"
        for stale in VariantCollection.objects.filter(name=vc_name):
            try:
                stale.delete_related_objects()
            except Exception:  # pragma: no cover
                pass
            stale.delete()

        vc = VariantCollection.objects.create(name=vc_name, status=ProcessingStatus.CREATED)
        vc_partition_table = vc.get_partition_table()

        union_sql_parts = [
            f'SELECT DISTINCT {vc.pk} AS variant_collection_id, variant_id '
            f'FROM "{p}" '
            f'WHERE collection_id = ANY(%s) AND samples_zygosity ~ %s'
            for p in partitions
        ]
        build_sql = (f'INSERT INTO "{vc_partition_table}" (variant_collection_id, variant_id) '
                     f'{" UNION ".join(union_sql_parts)};')
        coll_ids = [cgc.pk] + ([cgc.common_collection_id] if cgc.common_collection_id else [])
        params = []
        for _ in partitions:
            params.extend([coll_ids, regex_excl])

        try:
            t0 = time.perf_counter()
            with connection.cursor() as cur:
                cur.execute(build_sql, params)
                vc_record_count = cur.rowcount
            build_seconds = round(time.perf_counter() - t0, 4)
            vc.count = vc_record_count
            vc.status = ProcessingStatus.SUCCESS
            vc.save()
        except Exception as e:
            err_row = {
                "source": "synthetic", "analysis_id": "", "node_id": cohort_id,
                "node_type": "cohort_exclude_vc_join_postanalyze", "config_summary": config,
                "error": f"vc build: {e}",
            }
            try:
                vc.delete_related_objects()
            except Exception:  # pragma: no cover
                pass
            vc.delete()
            return [err_row]

        # ANALYZE the VCR partition — gives the planner accurate row counts and value
        # distribution for the just-INSERTed rows. Without this PG falls back to defaults
        # and may pick Hash(snpdb_variant 40M) instead of NestedLoop+IndexScan.
        old_autocommit = connection.get_autocommit()
        connection.set_autocommit(True)
        analyze_seconds = None
        try:
            t0 = time.perf_counter()
            with connection.cursor() as cur:
                cur.execute(f'ANALYZE "{vc_partition_table}";')
            analyze_seconds = round(time.perf_counter() - t0, 4)
        finally:
            connection.set_autocommit(old_autocommit)

        ann_kwargs = cgc.get_annotation_kwargs()
        ann_kwargs.update(vc.get_annotation_kwargs())
        cohort_alias = cgc.cohortgenotype_alias
        qs = (Variant.objects
              .annotate(**ann_kwargs)
              .filter(Q(**{f"{cohort_alias}__isnull": False}),
                      Q(**{f"{vc.variant_collection_alias}__isnull": False})))

        rows = []
        # Baseline post-ANALYZE: same SQL as cohort_exclude_vc_join above, but the planner
        # now has fresh stats — usually the dominant difference vs the non-diagnostic row.
        baseline_row = self._profile_synthetic_pattern(
            "cohort_exclude_vc_join_postanalyze", cohort_id, qs,
            f"{config} post_analyze",
            rerun=rerun, explain=explain, plans_dir=plans_dir,
            file_prefix="synthetic_c")
        baseline_row["build_seconds"] = build_seconds
        baseline_row["vc_record_count"] = vc.count
        baseline_row["analyze_seconds"] = analyze_seconds
        rows.append(baseline_row)

        # Tuning variants: same query, same fresh stats, additional session settings.
        for label_suffix, pg_settings in tuning_variants:
            with _PgSessionSettings(**pg_settings):
                variant_row = self._profile_synthetic_pattern(
                    f"cohort_exclude_vc_join_postanalyze_{label_suffix}", cohort_id, qs,
                    f"{config} post_analyze settings="
                    f"{'+'.join(f'{k}={v}' for k,v in pg_settings.items())}",
                    rerun=rerun, explain=explain, plans_dir=plans_dir,
                    file_prefix="synthetic_c")
                variant_row["vc_record_count"] = vc.count
            rows.append(variant_row)

        try:
            vc.delete_related_objects()
        except Exception:  # pragma: no cover
            pass
        vc.delete()
        return rows

    # ---- Helpers --------------------------------------------------------

    @staticmethod
    def _tag_subst(rows, subst_mode):
        for row in rows:
            row["pk_substitution"] = subst_mode
        return rows

    @staticmethod
    def _row_summary(row):
        bits = [
            row.get("source", ""),
            f"a{row.get('analysis_id') or '-'}",
            f"n{row.get('node_id') or '-'}",
            row.get("node_type", ""),
        ]
        if row.get("pk_substitution"):
            bits.append(f"subst={row['pk_substitution']}")
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


class _PgSessionSettings:
    """ Context manager — SET session-level PG settings around a block, RESET on exit.
    Used by --planner-diagnostic to test work_mem / random_page_cost variants on the
    cohort-exclude queries without touching cluster config (#1546). """

    def __init__(self, **settings):
        self.settings = settings

    def __enter__(self):
        if self.settings:
            with connection.cursor() as cur:
                for k, v in self.settings.items():
                    cur.execute(f"SET {k} = {v}")
        return self

    def __exit__(self, *_):
        if self.settings:
            with connection.cursor() as cur:
                for k in self.settings:
                    cur.execute(f"RESET {k}")


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
