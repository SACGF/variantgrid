import json
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError

from annotation.vep_columns import VEP_COLUMNS


class Command(BaseCommand):
    help = ("One-off: load column_vep_field.json and verify the new VEP_COLUMNS registry "
            "describes an identical set of (source_field, vep_plugin, vep_custom, "
            "source_field_has_custom_prefix, category, processing_description, summary_stats, "
            "version-bounds, genome_build, pipeline_type, variant_grid_column) tuples. "
            "Delete me + column_vep_field.json once green.")

    def add_arguments(self, parser):
        parser.add_argument("--json", default="column_vep_field.json",
                            help="Path to the JSON dump of the old ColumnVEPField table.")

    def handle(self, *args, **opts):
        json_path = Path(opts["json"])
        if not json_path.exists():
            raise CommandError(f"{json_path} not found")

        rows = json.loads(json_path.read_text())

        json_tuples = set()
        for r in rows:
            json_tuples.add((
                r["source_field"],
                r["vep_plugin"],
                r["vep_custom"],
                r["source_field_has_custom_prefix"],
                r["category"],
                r["source_field_processing_description"],
                r["summary_stats"],
                r["min_columns_version"],
                r["max_columns_version"],
                r["min_vep_version"],
                r["max_vep_version"],
                r["genome_build_id"],
                r["pipeline_type"],
                r["variant_grid_column_id"],
            ))

        registry_tuples = set()
        for c in VEP_COLUMNS:
            builds = sorted(c.genome_builds) or [None]
            pipelines = [p.value for p in sorted(c.pipeline_types, key=lambda p: p.value)] or [None]
            for build in builds:
                for pipeline in pipelines:
                    for vgc in c.variant_grid_columns:
                        registry_tuples.add((
                            c.source_field,
                            c.vep_plugin.value if c.vep_plugin else None,
                            c.vep_custom.value if c.vep_custom else None,
                            c.source_field_has_custom_prefix,
                            c.category.value,
                            c.source_field_processing_description,
                            c.summary_stats,
                            c.min_columns_version,
                            c.max_columns_version,
                            c.min_vep_version,
                            c.max_vep_version,
                            build,
                            pipeline,
                            vgc,
                        ))

        only_in_json = json_tuples - registry_tuples
        only_in_registry = registry_tuples - json_tuples

        if not only_in_json and not only_in_registry:
            self.stdout.write(self.style.SUCCESS(
                f"OK: {len(json_tuples)} mappings match across JSON and registry."
            ))
            return

        if only_in_json:
            self.stdout.write(self.style.ERROR(
                f"\n{len(only_in_json)} mapping(s) in JSON but missing from registry:"))
            for t in sorted(only_in_json, key=lambda x: tuple(str(v) for v in x)):
                self.stdout.write(f"  {t}")
        if only_in_registry:
            self.stdout.write(self.style.ERROR(
                f"\n{len(only_in_registry)} mapping(s) in registry but missing from JSON:"))
            for t in sorted(only_in_registry, key=lambda x: tuple(str(v) for v in x)):
                self.stdout.write(f"  {t}")

        raise CommandError("Registry does not match JSON. See diff above.")
