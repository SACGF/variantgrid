"""
Read-only aggregate data-mining for the VariantGrid paper.

>>> SAFETY CONTRACT (read before running) <<<
  * READ-ONLY: this command only counts/aggregates. It performs NO writes/migrations.
  * AGGREGATE ONLY: every output CSV is counts / rates / distributions. No patient-level rows.
  * NO IDENTIFIERS: never emits patient IDs, sample names, usernames, condition free-text, or
    person-linked HGVS/coordinates. Users appear only as *counts of distinct users*.
  * SMALL-CELL SUPPRESSION: any aggregate cell with count < --min-cell (default 5) is dropped.
  * Prefer running against a READ REPLICA or a restored backup, not live clinical prod.
  * Claude did not and will not run this; review every query yourself before running.

Install: copy this file to  <some_app>/management/commands/paper_data_mining.py
         (e.g. snpdb/management/commands/) then:
             python manage.py paper_data_mining --output-dir /tmp/vg_paper_stats
         Options:
             --min-cell N        suppression threshold (default 5)
             --only a,b,c        run only named collectors (default: all "safe" ones)
             --skip-heavy        skip slow collectors (genotype ratio, recurrence)
Targets CURRENT (VG4 / master) models. Run after migrating vg3_sapath_prod data up.
Lines marked `# VERIFY` use a field that may differ across schema versions — confirm in-context.
"""
import csv
import os
import traceback

from django.apps import apps
from django.core.management.base import BaseCommand
from django.db.models import Count
from django.db.models.functions import TruncMonth


def get_model(app_label, model_name):
    """Return a model or None (so a missing model just skips its collector)."""
    try:
        return apps.get_model(app_label, model_name)
    except LookupError:
        return None


class Command(BaseCommand):
    help = "Read-only AGGREGATE usage stats for the VariantGrid paper (no PII, small cells suppressed)."

    def add_arguments(self, parser):
        parser.add_argument("--output-dir", required=True)
        parser.add_argument("--min-cell", type=int, default=5)
        parser.add_argument("--only", default="", help="comma-separated collector names")
        parser.add_argument("--skip-heavy", action="store_true")

    def handle(self, *args, **opts):
        self.out_dir = opts["output_dir"]
        self.min_cell = opts["min_cell"]
        os.makedirs(self.out_dir, exist_ok=True)
        self.manifest = []

        collectors = {
            "ingestion_over_time": self.ingestion_over_time,
            "annotation_versions": self.annotation_versions,
            "annotation_version_diffs": self.annotation_version_diffs,
            "classifications_over_time": self.classifications_over_time,
            "significance_change_flags": self.significance_change_flags,
            "analysis_usage": self.analysis_usage,
            "node_type_usage": self.node_type_usage,
            "template_usage": self.template_usage,
            "variant_tagging": self.variant_tagging,
        }
        heavy = {
            "genotype_variant_ratio": self.genotype_variant_ratio,
            "classified_variant_recurrence": self.classified_variant_recurrence,
        }
        if not opts["skip_heavy"]:
            collectors.update(heavy)

        only = {c.strip() for c in opts["only"].split(",") if c.strip()}
        for name, fn in collectors.items():
            if only and name not in only:
                continue
            try:
                rows = fn()
                self.manifest.append((name, "ok", rows))
                self.stdout.write(self.style.SUCCESS(f"{name}: {rows} rows"))
            except Exception as e:  # noqa: BLE001 — one collector failing must not kill the rest
                self.manifest.append((name, f"FAILED: {e}", 0))
                self.stderr.write(self.style.ERROR(f"{name} FAILED: {e}"))
                self.stderr.write(traceback.format_exc())

        self._write("_run_manifest.csv", ["collector", "status", "rows_emitted"], self.manifest)
        self.stdout.write(self.style.SUCCESS(f"Done. CSVs in {self.out_dir} (min_cell={self.min_cell})."))

    # ---- helpers -----------------------------------------------------------
    def _write(self, filename, header, rows):
        path = os.path.join(self.out_dir, filename)
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(header)
            w.writerows(rows)
        return path

    def _suppress(self, count):
        """Return count, or '<min' if below the suppression threshold (and non-zero)."""
        return count if (count == 0 or count >= self.min_cell) else f"<{self.min_cell}"

    @staticmethod
    def _month(dt):
        return dt.strftime("%Y-%m") if dt else "unknown"

    # ---- collectors --------------------------------------------------------
    def ingestion_over_time(self):
        """VCFs and samples ingested per month. Source: snpdb.VCF.date, Sample via vcf."""
        VCF = get_model("snpdb", "VCF")
        Sample = get_model("snpdb", "Sample")
        rows = []
        vcf_by_month = (VCF.objects.annotate(m=TruncMonth("date"))  # VERIFY: VCF.date populated on import
                        .values("m").annotate(n=Count("id")).order_by("m"))
        # samples per month keyed off their VCF's date (Sample has no own timestamp)
        sample_by_month = (Sample.objects.annotate(m=TruncMonth("vcf__date"))
                           .values("m").annotate(n=Count("id")).order_by("m"))
        sbm = {self._month(r["m"]): r["n"] for r in sample_by_month}
        for r in vcf_by_month:
            mth = self._month(r["m"])
            rows.append([mth, self._suppress(r["n"]), self._suppress(sbm.get(mth, 0))])
        self._write("ingestion_by_month.csv", ["month", "vcfs", "samples"], rows)
        return len(rows)

    def genotype_variant_ratio(self):
        """Total variants vs total genotype observations (the 2.4B-vs-46M style figure). HEAVY."""
        Variant = get_model("snpdb", "Variant")
        CGC = get_model("snpdb", "CohortGenotypeCollection")
        CG = get_model("snpdb", "CohortGenotype")
        num_variants = Variant.objects.count()
        # genotype observations = sum over collections of (rows in collection * num_samples)
        obs = 0
        for cgc in CGC.objects.all().only("id", "num_samples"):
            n_rows = CG.objects.filter(collection_id=cgc.id).count()  # VERIFY: CG.collection FK name
            obs += n_rows * (cgc.num_samples or 0)
        ratio = round(obs / num_variants, 2) if num_variants else 0
        rows = [[num_variants, obs, ratio]]
        self._write("genotype_variant_ratio.csv",
                    ["num_variants", "num_genotype_observations", "genotypes_per_variant"], rows)
        return len(rows)

    def annotation_versions(self):
        VAV = get_model("annotation", "VariantAnnotationVersion")
        rows = []
        for v in VAV.objects.all().order_by("created"):
            rows.append([v.id, self._month(getattr(v, "created", None)),
                         self._month(getattr(v, "annotation_date", None)),
                         getattr(v, "columns_version", ""), getattr(v, "status", "")])
        self._write("annotation_versions.csv",
                    ["version_id", "created_month", "annotation_month", "columns_version", "status"], rows)
        return len(rows)

    def annotation_version_diffs(self):
        """Added/modified/removed/unchanged between annotation versions — backs reanalysis claim."""
        VD = get_model("annotation", "VariantAnnotationVersionDiff")
        rows = []
        for d in VD.objects.all():
            rows.append([getattr(d, "version_from_id", ""), getattr(d, "version_to_id", ""),
                         getattr(d, "num_added", ""), getattr(d, "num_modified", ""),
                         getattr(d, "num_removed", ""), getattr(d, "num_unchanged", "")])
        self._write("annotation_version_diffs.csv",
                    ["version_from", "version_to", "num_added", "num_modified",
                     "num_removed", "num_unchanged"], rows)
        # Per-column change detail (e.g. ClinVar significance churn). Best-effort: introspect the
        # FromToResult model and emit aggregate change counts per column. # VERIFY field names.
        FTR = get_model("annotation", "VersionDiffFromToResult")
        if FTR is not None:
            fields = {f.name for f in FTR._meta.get_fields()}
            col_field = next((c for c in ("column", "column_name", "field", "key") if c in fields), None)
            if col_field:
                detail = (FTR.objects.values(col_field)
                          .annotate(n=Count("id")).order_by("-n"))
                drows = [[r[col_field], r["n"]] for r in detail]
                self._write("annotation_version_diff_columns.csv", ["column", "change_rows"], drows)
        return len(rows)

    def classifications_over_time(self):
        """LIGHT (Shariant's turf): counts by month + clinical significance. Suppressed."""
        C = get_model("classification", "Classification")
        agg = (C.objects.annotate(m=TruncMonth("created"))
               .values("m", "clinical_significance").annotate(n=Count("id")).order_by("m"))
        rows = [[self._month(r["m"]), r["clinical_significance"] or "none", self._suppress(r["n"])]
                for r in agg]
        self._write("classifications_by_month.csv", ["month", "clinical_significance", "count"], rows)
        return len(rows)

    def significance_change_flags(self):
        """Reclassification-over-time evidence: significance-change flags per month."""
        Flag = get_model("flags", "Flag")
        agg = (Flag.objects.filter(flag_type__id="classification_significance_change")  # VERIFY id
               .annotate(m=TruncMonth("created")).values("m").annotate(n=Count("id")).order_by("m"))
        rows = [[self._month(r["m"]), self._suppress(r["n"])] for r in agg]
        self._write("significance_changes_by_month.csv", ["month", "count"], rows)
        return len(rows)

    def analysis_usage(self):
        """Analyses created per month + distinct users (COUNT only, never identities)."""
        A = get_model("analysis", "Analysis")
        agg = (A.objects.annotate(m=TruncMonth("created")).values("m")
               .annotate(n=Count("id"), users=Count("user", distinct=True)).order_by("m"))
        rows = [[self._month(r["m"]), self._suppress(r["n"]), self._suppress(r["users"])] for r in agg]
        self._write("analyses_by_month.csv", ["month", "analyses", "distinct_users"], rows)
        return len(rows)

    def node_type_usage(self):
        """Histogram of AnalysisNode subclass usage (MTI). Backs 'which nodes people use'."""
        AnalysisNode = get_model("analysis", "AnalysisNode")
        rows = []
        for m in apps.get_app_config("analysis").get_models():
            if issubclass(m, AnalysisNode) and m is not AnalysisNode and not m._meta.abstract:
                try:
                    n = m.objects.count()
                except Exception:  # noqa: BLE001 — proxy/abstract edge cases
                    continue
                if n:
                    rows.append([m.__name__, self._suppress(n)])
        rows.sort(key=lambda r: (r[1] == f"<{self.min_cell}", r[1]), reverse=True)
        self._write("node_type_usage.csv", ["node_type", "count"], rows)
        return len(rows)

    def template_usage(self):
        """Template runs per month — supports 'exploration vs fixed pipeline' thesis."""
        TR = get_model("analysis", "AnalysisTemplateRun")
        if TR is None:
            return 0
        agg = (TR.objects.annotate(m=TruncMonth("created"))  # VERIFY: AnalysisTemplateRun has created
               .values("m").annotate(n=Count("id")).order_by("m"))
        rows = [[self._month(r["m"]), self._suppress(r["n"])] for r in agg]
        self._write("template_runs_by_month.csv", ["month", "template_runs"], rows)
        return len(rows)

    def variant_tagging(self):
        """Tagging behaviour: counts by tag type (suppressed) and by month."""
        VT = get_model("analysis", "VariantTag")
        by_tag = VT.objects.values("tag_id").annotate(n=Count("id")).order_by("-n")
        rows = [[r["tag_id"], self._suppress(r["n"])] for r in by_tag]
        self._write("variant_tags_by_type.csv", ["tag", "count"], rows)
        by_month = (VT.objects.annotate(m=TruncMonth("created")).values("m")
                    .annotate(n=Count("id")).order_by("m"))
        mrows = [[self._month(r["m"]), self._suppress(r["n"])] for r in by_month]
        self._write("variant_tags_by_month.csv", ["month", "count"], mrows)
        return len(rows) + len(mrows)

    def classified_variant_recurrence(self):
        """How often a classified allele is carried by samples (internal-population-DB value). HEAVY.
        Outputs only a single aggregate distribution; verify performance on a replica first."""
        C = get_model("classification", "Classification")
        VZC = get_model("snpdb", "VariantZygosityCount")
        # distinct alleles that have at least one classification
        classified_alleles = (C.objects.exclude(allele__isnull=True)
                              .values_list("allele_id", flat=True).distinct())
        n_classified = classified_alleles.count()
        # of those, how many are carried (het+hom) by >=1 sample, via zygosity counts on their variants
        # VERIFY: linkage allele -> variant -> VariantZygosityCount.variant; this is a coarse proxy.
        carried = 0
        # iterate in batches to bound memory; this is intentionally simple/auditable
        for allele_id in classified_alleles.iterator(chunk_size=1000):
            has = (VZC.objects.filter(variant__variantallele__allele_id=allele_id)  # VERIFY relation
                   .filter(het_count__gt=0).exists() or
                   VZC.objects.filter(variant__variantallele__allele_id=allele_id)
                   .filter(hom_count__gt=0).exists())
            if has:
                carried += 1
        rows = [[n_classified, carried,
                 round(100.0 * carried / n_classified, 1) if n_classified else 0]]
        self._write("classified_variant_recurrence.csv",
                    ["classified_alleles", "carried_by_a_sample", "percent_recurrent"], rows)
        return len(rows)
