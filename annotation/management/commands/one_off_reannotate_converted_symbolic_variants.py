"""
Re-annotates, through the STRUCTURAL_VARIANT pipeline, variants whose stored annotation came from the STANDARD
pipeline although their version now routes them to STRUCTURAL_VARIANT (#2786). These are the explicit del/dup/inv
that one_off_fix_symbolic_variants converted in place to a symbolic <DEL>/<DUP>/<INV> at or over the version's
structural_variant_min_size: the conversion kept the pk, and so the annotation VEP wrote for the explicit
sequence, and the scheduler only ever fills new pk ranges.

For each NEW/ACTIVE version, the variants are dumped to one VCF and run through VEP with the STRUCTURAL_VARIANT
pipeline's command line (get_vep_command), here rather than as an AnnotationRun. Then, one transaction per range lock,
their STANDARD rows (VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap) are deleted and VEP's
output is inserted by BulkVEPVCFAnnotationInserter, attributed to the lock's existing finished STRUCTURAL_VARIANT
run - so every row matches its run's pipeline_type again, no run or range lock is added, and a second run of this
command finds nothing. Where the lock's AnnotSV run has finished at the ACTIVE AnnotSV version, AnnotSV is run over
the same variants and its TSV imported through that run as well; an AnnotSV run that is still to come takes them
with every other SV in its range. The run rows themselves are never saved: they keep the history of the run that
annotated the rest of the range.

A lock with no finished STRUCTURAL_VARIANT run is reported and left alone. HISTORICAL versions are counted but not
touched: their VEP no longer matches the installed one. --dry-run prints the counts and the runs the rows would go to.
Registered as a deploy-time ManualOperation (annotation/migrations/0017_reannotate_converted_symbolic_variants.py).
"""
import os
import shutil
import tempfile
from dataclasses import dataclass, field
from typing import Optional

import cyvcf2
from django.conf import settings
from django.core.management.base import BaseCommand
from django.db import transaction

from annotation.annotation_run_files import write_qs_to_vcf
from annotation.annotation_version_querysets import (
    filter_vep_sv_max_size,
    get_queryset_for_annotation_version,
    get_variants_qs_for_annotation,
)
from annotation.annotsv_annotation import get_annotsv_command, get_annotsv_tsv_filename
from annotation.enums import AnnotationStatus, VariantAnnotationPipelineType, VEPSkippedReason
from annotation.models.models import (
    VariantAnnotation,
    VariantGeneOverlap,
    VariantTranscriptAnnotation,
)
from annotation.models.models_annotation_run import AnnotationPipelineVersion, AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion
from annotation.pipelines import get_pipeline, get_runner
from annotation.sv_conservation import (
    conservation_sidecar_filename,
    get_sv_conservation_tracks,
    score_sv_vcf,
    write_conservation_sidecar,
)
from annotation.vcf_files.bulk_annotsv_tsv_inserter import import_annotsv_tsv
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import BulkVEPVCFAnnotationInserter
from annotation.vcf_files.import_vcf_annotations import (
    _add_vep_skipped_variant,
    _insert_vep_skipped_variants,
)
from annotation.vep_annotation import (
    get_vep_command,
    get_vep_skipped_variants_filename,
    vep_check_annotated_file_version_match,
    vep_check_command_line_version_match,
)
from annotation.vep_command_check import check_vep_command
from annotation.vep_warnings import parse_skipped_variants_file
from library.django_utils.data_archive_mixin import DataArchivedError
from library.utils import execute_cmd
from snpdb.models import Variant
from snpdb.vcf_utils import cyvcf2_header_types

STANDARD = VariantAnnotationPipelineType.STANDARD
STRUCTURAL_VARIANT = VariantAnnotationPipelineType.STRUCTURAL_VARIANT
ANNOTSV = VariantAnnotationPipelineType.ANNOTSV


@dataclass
class LockGroup:
    """ Misrouted variants annotated by one range lock's STANDARD run, and the lock's runs their rows go to """
    range_lock_id: int
    variant_ids: list[int] = field(default_factory=list)
    standard_run_ids: set[int] = field(default_factory=set)
    sv_run: Optional[AnnotationRun] = None
    annotsv_run: Optional[AnnotationRun] = None

    def __str__(self):
        s = f"range lock {self.range_lock_id}: {len(self.variant_ids)} variants"
        if self.sv_run:
            s += f" -> STRUCTURAL_VARIANT run {self.sv_run.pk}"
            if self.annotsv_run:
                s += f", AnnotSV run {self.annotsv_run.pk}"
        else:
            s += " - no finished STRUCTURAL_VARIANT run, left as they are"
        return s


def find_misrouted_variants(annotation_version) -> list[LockGroup]:
    """ Variants the version routes to STRUCTURAL_VARIANT that hold a STANDARD run's annotation """
    qs = get_variants_qs_for_annotation(annotation_version, pipeline_type=STRUCTURAL_VARIANT, annotated=True)
    qs = qs.filter(variantannotation__annotation_run__pipeline_type=STANDARD)
    rows = qs.order_by("pk").values_list("pk", "variantannotation__annotation_run_id",
                                         "variantannotation__annotation_run__annotation_range_lock_id")
    groups: dict[int, LockGroup] = {}
    for variant_id, run_id, range_lock_id in rows:
        group = groups.setdefault(range_lock_id, LockGroup(range_lock_id))
        group.variant_ids.append(variant_id)
        group.standard_run_ids.add(run_id)

    finished_runs_qs = AnnotationRun.objects.filter(annotation_range_lock__in=list(groups),
                                                    status=AnnotationStatus.FINISHED)
    for run in finished_runs_qs.filter(pipeline_type=STRUCTURAL_VARIANT).select_related("annotation_range_lock__version"):
        groups[run.annotation_range_lock_id].sv_run = run
    annotsv_runs_qs = finished_runs_qs.filter(pipeline_type=ANNOTSV,
                                              pipeline_version__status=AnnotationPipelineVersion.Status.ACTIVE)
    for run in annotsv_runs_qs.select_related("annotation_range_lock__version", "pipeline_version"):
        groups[run.annotation_range_lock_id].annotsv_run = run
    return list(groups.values())


def run_sv_vep(vav: VariantAnnotationVersion, vcf_dump_filename: str, vcf_annotated_filename: str):
    """ VEP over the dump with the STRUCTURAL_VARIANT pipeline's arguments, plus the conservation sidecar its
        import reads, as VEPRunner.annotate does for a run """
    vep_check_command_line_version_match(vav)
    cmd = get_vep_command(vcf_dump_filename, vcf_annotated_filename, vav.genome_build, vav.annotation_consortium,
                          STRUCTURAL_VARIANT, variant_annotation_version=vav)
    check_vep_command(vav, STRUCTURAL_VARIANT, cmd)
    return_code, _std_out, std_err = execute_cmd(cmd)
    if return_code != 0:
        raise RuntimeError(f"VEP returned {return_code}: {std_err}")
    if get_runner(STRUCTURAL_VARIANT).sv_conservation_enabled:
        genome_build = vav.genome_build
        write_conservation_sidecar(conservation_sidecar_filename(vcf_annotated_filename),
                                   score_sv_vcf(vcf_dump_filename, genome_build),
                                   get_sv_conservation_tracks(genome_build))


def run_annotsv(annotsv_run: AnnotationRun, vcf_dump_filename: str, output_dir: str) -> str:
    """ AnnotSV over the dump, as AnnotSVRunner.annotate does for a run. Returns its TSV """
    get_runner(ANNOTSV).check_tool_version(annotsv_run)
    cmd = get_annotsv_command(vcf_dump_filename, output_dir, annotsv_run.genome_build,
                              annotsv_run.annotation_consortium)
    return_code, _std_out, std_err = execute_cmd(cmd, timeout=settings.ANNOTATION_ANNOTSV_TIMEOUT_SECONDS)
    tsv_filename = get_annotsv_tsv_filename(vcf_dump_filename, output_dir)
    if return_code != 0 or not os.path.exists(tsv_filename):
        raise RuntimeError(f"AnnotSV returned {return_code}: {std_err}")
    return tsv_filename


def check_vep_output(vcf_annotated_filename: str, dump_count: int):
    """ Every dumped record must come back annotated or listed as skipped - VEP can exit 0 having written a short
        VCF (see import_vcf_annotations.check_all_dumped_records_accounted_for) """
    records_read = sum(1 for _ in cyvcf2.VCF(vcf_annotated_filename))
    skipped = len(parse_skipped_variants_file(get_vep_skipped_variants_filename(vcf_annotated_filename)))
    if records_read + skipped != dump_count:
        raise ValueError(f"VEP output '{vcf_annotated_filename}' has {records_read} records and {skipped} skipped, "
                         f"but {dump_count} were dumped - truncated")


def reannotate_group(annotation_version, group: LockGroup, vcf_annotated_filename: Optional[str],
                     annotsv_tsv_filename: Optional[str]):
    """ Swap the group's STANDARD rows for VEP's STRUCTURAL_VARIANT output, attributed to the lock's SV run """
    variant_ids = set(group.variant_ids)
    sv_run = group.sv_run
    with transaction.atomic():
        for klass in [VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap]:
            qs = get_queryset_for_annotation_version(klass, annotation_version)
            qs.filter(variant_id__in=variant_ids, annotation_run_id__in=group.standard_run_ids).delete()

        # The inserter reads the annotated VCF (and its conservation sidecar) off the run - this instance is never saved
        sv_run.vcf_annotated_filename = vcf_annotated_filename
        annotated_ids = set()
        if vcf_annotated_filename:
            vcf_reader = cyvcf2.VCF(vcf_annotated_filename)
            bulk_inserter = BulkVEPVCFAnnotationInserter(sv_run, cyvcf2_header_types(vcf_reader)["INFO"])
            bulk_inserter.remove_processing_files()  # Its scratch dir is the run's - clear what the run left
            for v in vcf_reader:
                variant_id = v.INFO["variant_id"]
                if variant_id in variant_ids:
                    bulk_inserter.process_entry(v)
                    if v.INFO.get("CSQ"):
                        annotated_ids.add(variant_id)
            bulk_inserter.finish()
        else:
            bulk_inserter = BulkVEPVCFAnnotationInserter(sv_run, validate_columns=False, vep_skipped_only=True)
            bulk_inserter.remove_processing_files()

        sv_max_size = settings.ANNOTATION_VEP_SV_MAX_SIZE
        skipped_qs = Variant.objects.filter(pk__in=variant_ids - annotated_ids)
        for v in skipped_qs.select_related("locus__contig", "locus__ref", "alt"):
            too_long = sv_max_size and v.svlen is not None and abs(v.svlen) > sv_max_size  # never dumped
            reason = VEPSkippedReason.TOO_LONG if too_long else VEPSkippedReason.UNKNOWN
            _add_vep_skipped_variant(sv_run, bulk_inserter, v, reason)
        _insert_vep_skipped_variants(bulk_inserter)
        bulk_inserter.remove_processing_files()

        if annotsv_tsv_filename and group.annotsv_run:
            group.annotsv_run.vcf_annotated_filename = annotsv_tsv_filename  # never saved, as above
            import_annotsv_tsv(group.annotsv_run)


def reannotate_version(annotation_version, groups: list[LockGroup], work_dir: str) -> int:
    """ Returns the number of variants AnnotSV was run over """
    vav = annotation_version.variant_annotation_version
    genome_build = vav.genome_build
    variant_ids = [variant_id for group in groups for variant_id in group.variant_ids]
    variants_qs = Variant.objects.filter(pk__in=variant_ids)

    dump_qs = variants_qs
    if settings.ANNOTATION_VEP_SV_MAX_SIZE:
        dump_qs = filter_vep_sv_max_size(dump_qs, settings.ANNOTATION_VEP_SV_MAX_SIZE, too_long=False)
    vcf_dump_filename = os.path.join(work_dir, "structural_variant.vcf")
    dump_count = write_qs_to_vcf(vcf_dump_filename, genome_build, dump_qs,
                                 structural_variant_min_size=vav.structural_variant_min_size)
    vcf_annotated_filename = None
    if dump_count:
        vcf_annotated_filename = os.path.join(work_dir, f"structural_variant.vep_annotated_{genome_build.name}.vcf.gz")
        run_sv_vep(vav, vcf_dump_filename, vcf_annotated_filename)
        vep_check_annotated_file_version_match(vav, vcf_annotated_filename)
        check_vep_output(vcf_annotated_filename, dump_count)

    annotsv_tsv_filename = None
    annotsv_groups = [group for group in groups if group.annotsv_run]
    num_annotsv = sum(len(group.variant_ids) for group in annotsv_groups)
    if annotsv_groups:
        annotsv_dir = os.path.join(work_dir, "annotsv")
        os.makedirs(annotsv_dir)
        annotsv_dump_filename = os.path.join(annotsv_dir, "annotsv.vcf")
        annotsv_ids = [variant_id for group in annotsv_groups for variant_id in group.variant_ids]
        write_qs_to_vcf(annotsv_dump_filename, genome_build, Variant.objects.filter(pk__in=annotsv_ids),
                        samples=get_runner(ANNOTSV).dump_samples,
                        structural_variant_min_size=vav.structural_variant_min_size)
        annotsv_tsv_filename = run_annotsv(annotsv_groups[0].annotsv_run, annotsv_dump_filename, annotsv_dir)

    for group in groups:
        reannotate_group(annotation_version, group, vcf_annotated_filename, annotsv_tsv_filename)
    return num_annotsv


class Command(BaseCommand):
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument("--dry-run", action="store_true",
                            help="Print the counts and the runs the rows would go to, without running VEP or writing")
        parser.add_argument("--variant-annotation-version", type=int,
                            help="Only this version (default: every version whose data isn't archived)")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        qs = VariantAnnotationVersion.objects.filter(data_archived_date__isnull=True).select_related("genome_build")
        if vav_id := options["variant_annotation_version"]:
            qs = qs.filter(pk=vav_id)
        for vav in qs.order_by("pk"):
            prefix = "DRY RUN " if dry_run else ""
            self.stdout.write(f"{prefix}{vav}: {self._handle_version(vav, dry_run)}")

    def _handle_version(self, vav: VariantAnnotationVersion, dry_run: bool) -> str:
        annotation_version = vav.get_any_annotation_version()
        if annotation_version is None:
            return "skipped: no AnnotationVersion"
        try:
            groups = find_misrouted_variants(annotation_version)
        except DataArchivedError as e:
            return f"skipped: {e}"
        num_variants = sum(len(g.variant_ids) for g in groups)
        summary = f"{num_variants} variants with STANDARD annotation now routed to STRUCTURAL_VARIANT"
        if not num_variants:
            return summary
        if vav.status == VariantAnnotationVersion.Status.HISTORICAL:
            return f"{summary} - left as they are: HISTORICAL versions aren't re-annotated"
        if blocker := vav.get_annotation_run_blocker():
            return f"{summary} - left as they are until annotation is unblocked: {blocker}"

        annotsv = get_pipeline(ANNOTSV)
        if not (annotsv.enabled and annotsv.runner.supports_genome_build(vav.genome_build)):
            for group in groups:
                group.annotsv_run = None
        ready = [group for group in groups if group.sv_run]
        if dry_run:
            for group in groups:
                self.stdout.write(f"  {group}")
        num_ready = sum(len(g.variant_ids) for g in ready)
        summary += f", {num_ready} in {len(ready)} range locks with a finished STRUCTURAL_VARIANT run"
        if dry_run or not ready:
            return summary

        os.makedirs(settings.ANNOTATION_VCF_DUMP_DIR, exist_ok=True)
        work_dir = tempfile.mkdtemp(prefix=f"reannotate_converted_symbolic_vav{vav.pk}_",
                                    dir=settings.ANNOTATION_VCF_DUMP_DIR)
        num_annotsv = reannotate_version(annotation_version, ready, work_dir)
        shutil.rmtree(work_dir)  # Left in place for a look when anything above raised
        return f"{summary}: re-annotated (AnnotSV over {num_annotsv})"
