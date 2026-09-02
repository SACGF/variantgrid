#!/usr/bin/env python3
"""
Backfill VariantAnnotation columns from an annotated VCF (#1675).

--dump writes the variants already annotated in a version to a VCF carrying variant_id in INFO; the
operator annotates that file however suits the column (bcftools annotate against a COSMIC VCF, or VEP);
--import reads it back and updates the named columns.

See claude/plans/1675_backfill_annotation_column_plan.md.
"""
import os

from django.core.management.base import BaseCommand, CommandError

from annotation.backfill_columns import (
    DEFAULT_BATCH_SIZE,
    BackfillColumnError,
    dump_annotated_variants,
    import_backfill_vcf,
    resolve_backfill_columns,
)
from annotation.models.models import VariantAnnotationVersion
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    help = "Dump/import a VCF to backfill VariantAnnotation columns (#1675)"

    def add_arguments(self, parser):
        mode = parser.add_mutually_exclusive_group(required=True)
        mode.add_argument("--dump", action="store_true",
                          help="Write the version's annotated variants to a VCF for external annotation")
        mode.add_argument("--import", action="store_true", dest="import_mode",
                          help="Read an annotated VCF back, updating --columns")

        parser.add_argument("--genome-build", required=True)
        parser.add_argument("--variant-annotation-version", type=int,
                            help="VariantAnnotationVersion pk (default: the build's ACTIVE version)")
        parser.add_argument("--columns", required=True,
                            help="Comma separated VariantGrid column ids, eg 'cosmic_count,cosmic_legacy_id'")
        parser.add_argument("--info", action="append", default=[], metavar="COLUMN=INFO_FIELD",
                            help="Read a column from this INFO field instead of the registry default "
                                 "(the source file's naming varies by release). Repeatable")
        parser.add_argument("--csq", action="append", default=[], metavar="COLUMN=CSQ_FIELD",
                            help="Read a column off the PICK'd VEP CSQ entry rather than a top level "
                                 "INFO field. Repeatable")
        parser.add_argument("--output", help="--dump: VCF to write (.gz writes bgzip, which tabix can index)")
        parser.add_argument("--only-missing", action="store_true",
                            help="--dump: only variants where every --columns value is null")
        parser.add_argument("--not-null", dest="not_null",
                            help="--dump: only variants that already have a value for these comma "
                                 "separated columns, eg a sibling column the same source writes - "
                                 "keeps the dump to the variants that source can say anything about")
        parser.add_argument("--min-variant-id", type=int,
                            help="--dump: split a large dump into pieces that can be annotated in parallel")
        parser.add_argument("--max-variant-id", type=int, help="--dump: see --min-variant-id")
        parser.add_argument("--vcf", help="--import: the externally annotated VCF")
        parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
                            help="--import: records per COPY/UPDATE cycle (default %(default)s)")
        parser.add_argument("--clear-missing", action="store_true",
                            help="--import: write NULL where the source has no value, making the file "
                                 "authoritative for the dumped variants (default: leave it alone)")
        parser.add_argument("--dry-run", action="store_true",
                            help="--import: report how many rows would change, changing nothing")

    def handle(self, *args, **options):
        genome_build = GenomeBuild.get_name_or_alias(options["genome_build"])
        variant_annotation_version = self._get_variant_annotation_version(genome_build, options)
        column_ids = [c.strip() for c in options["columns"].split(",") if c.strip()]
        if not column_ids:
            raise CommandError("--columns needs at least one VariantGrid column id")

        try:
            targets = resolve_backfill_columns(variant_annotation_version, column_ids,
                                               info_overrides=self._parse_overrides(options["info"], "--info"),
                                               csq_overrides=self._parse_overrides(options["csq"], "--csq"))
        except BackfillColumnError as e:
            raise CommandError(str(e))

        for target in targets:
            tables = ", ".join(target.base_table_names)
            self.stdout.write(f"{target.column} <- {target.source_description} ({tables})")

        if options["dump"]:
            self._run_dump(variant_annotation_version, targets, options)
        else:
            self._run_import(variant_annotation_version, targets, options)

    @staticmethod
    def _get_variant_annotation_version(genome_build, options) -> VariantAnnotationVersion:
        if pk := options["variant_annotation_version"]:
            try:
                variant_annotation_version = VariantAnnotationVersion.objects.get(pk=pk)
            except VariantAnnotationVersion.DoesNotExist:
                raise CommandError(f"No VariantAnnotationVersion with pk={pk}")
            if variant_annotation_version.genome_build != genome_build:
                raise CommandError(f"VariantAnnotationVersion pk={pk} is "
                                   f"{variant_annotation_version.genome_build}, not {genome_build}")
        else:
            variant_annotation_version = VariantAnnotationVersion.latest(genome_build)
            if variant_annotation_version is None:
                raise CommandError(f"No ACTIVE VariantAnnotationVersion for {genome_build} - name one "
                                   f"with --variant-annotation-version")
        return variant_annotation_version

    @staticmethod
    def _parse_overrides(values, arg_name) -> dict[str, str]:
        overrides = {}
        for value in values:
            column, sep, source_field = value.partition("=")
            if not (sep and column.strip() and source_field.strip()):
                raise CommandError(f"{arg_name} '{value}' must be COLUMN=FIELD")
            overrides[column.strip()] = source_field.strip()
        return overrides

    def _run_dump(self, variant_annotation_version, targets, options):
        output_filename = options.get("output")
        if not output_filename:
            raise CommandError("--dump requires --output")

        not_null_columns = [c.strip() for c in (options["not_null"] or "").split(",") if c.strip()]
        try:
            count = dump_annotated_variants(variant_annotation_version, output_filename, targets=targets,
                                            only_missing=options["only_missing"],
                                            not_null_columns=not_null_columns,
                                            min_variant_id=options["min_variant_id"],
                                            max_variant_id=options["max_variant_id"])
        except BackfillColumnError as e:
            raise CommandError(str(e))
        self.stdout.write(f"Dumped {count} variants annotated by {variant_annotation_version} "
                          f"to '{output_filename}'")

    def _run_import(self, variant_annotation_version, targets, options):
        vcf_filename = options.get("vcf")
        if not vcf_filename:
            raise CommandError("--import requires --vcf")
        if not os.path.exists(vcf_filename):
            raise CommandError(f"--vcf '{vcf_filename}' does not exist")

        def emit(message):
            self.stdout.write(message)
            self.stdout.flush()

        try:
            report = import_backfill_vcf(variant_annotation_version, vcf_filename, targets,
                                         batch_size=options["batch_size"],
                                         clear_missing=options["clear_missing"],
                                         dry_run=options["dry_run"], emit=emit)
        except BackfillColumnError as e:
            raise CommandError(str(e))

        verb = "would change" if options["dry_run"] else "changed"
        self.stdout.write(f"{variant_annotation_version}: read {report['records_read']} records, "
                          f"{verb} {report['rows_changed']} rows")
