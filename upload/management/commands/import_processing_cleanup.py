"""
Sweeps settings.IMPORT_PROCESSING_DIR of scratch directories whose owner has finished with them (#928).

Everything under that directory is disposable working space written by an import, but until #928 several
paths never removed theirs, so an existing deployment has years of them. Each owner now cleans up as it
finishes (UploadPipeline.success and its post_delete, the annotation run reset, the SQL COPY commands),
so this is mostly a one-off catch-up - registered as a ManualOperation in
upload/migrations/0042_one_off_import_processing_cleanup.py - but it is worth keeping: a deployment
running with IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS off still accumulates them.

An entry is only removed when its owner says it is done with it. Anything still owned (a PROCESSING
pipeline, an ERROR annotation run kept for investigation), anything modified within --min-age-days, and
anything whose name this does not recognise is reported and left alone.
"""
import os
import shutil
import time

from django.conf import settings
from django.core.management.base import BaseCommand

from annotation.models import (
    AnnotationRun,
    GeneAnnotationVersion,
    HumanProteinAtlasAnnotationVersion,
)
from annotation.models.models_enums import AnnotationStatus
from snpdb.models import SomalierVCFExtract
from snpdb.models.models_enums import ProcessingStatus
from upload.models import (
    UploadedClassificationImport,
    UploadedLiftover,
    UploadedManualVariantEntryCollection,
    UploadPipeline,
)

DAY_SECS = 24 * 60 * 60


def _dir_size(path: str) -> int:
    total = 0
    for dir_path, _dir_names, filenames in os.walk(path):
        for filename in filenames:
            try:
                total += os.lstat(os.path.join(dir_path, filename)).st_size
            except OSError:
                pass
    return total


def _pipeline_finished(file_upload_id=None, pk=None) -> bool:
    """ A pipeline that's gone, or succeeded, has no further use for its files. ERROR/PROCESSING keeps
        them - that is what 'Retry import' reloads from, and what a failure is investigated with """
    if pk is not None:
        qs = UploadPipeline.objects.filter(pk=pk)
    else:
        qs = UploadPipeline.objects.filter(file_upload_id=file_upload_id)
    status = qs.values_list("status", flat=True).first()
    return status is None or status == ProcessingStatus.SUCCESS


def _upload_pipeline_finished(pk: int) -> bool:
    """ pipeline_<pk> - split VCF chunks, per-step logs and COPY CSVs of a VCF import """
    return _pipeline_finished(pk=pk)


def _annotation_run_finished(pk: int) -> bool:
    """ annotation_run_<pk> - the CSVs BulkVEPVCFAnnotationInserter SQL COPYs from """
    if annotation_run := AnnotationRun.objects.filter(pk=pk).first():
        return annotation_run.get_status() == AnnotationStatus.FINISHED
    return True


def _generated_input_finished(upload_data) -> bool:
    """ liftover_<pk> / manual_variants_<pk> / classification_import_<pk> hold a VCF we generated as the
        input to a pipeline, so they live as long as that pipeline wants to be able to re-read it """
    if upload_data is None:
        return True
    return _pipeline_finished(file_upload_id=upload_data.file_upload_id)


def _liftover_finished(pk: int) -> bool:
    return _generated_input_finished(UploadedLiftover.objects.filter(liftover_id=pk).first())


def _manual_variants_finished(pk: int) -> bool:
    return _generated_input_finished(UploadedManualVariantEntryCollection.objects.filter(collection_id=pk).first())


def _classification_import_finished(pk: int) -> bool:
    return _generated_input_finished(UploadedClassificationImport.objects.filter(classification_import_id=pk).first())


def _somalier_vcf_extract_finished(pk: int) -> bool:
    status = SomalierVCFExtract.objects.filter(pk=pk).values_list("status", flat=True).first()
    return status is None or status == ProcessingStatus.SUCCESS


def _somalier_relate_finished(_pk: int) -> bool:
    """ snpdb.tasks.somalier_tasks._somalier_relate removes its own dir the moment somalier returns, so
        one that survived belongs to a run that died. SomalierRelate is abstract (a table per subclass)
        so the pk in the name doesn't identify a row anyway - the age guard is what keeps a live run safe """
    return True


def _gene_annotation_finished(pk: int) -> bool:
    """ The COPY runs after the version row is created, so a row means the import reached the insert """
    return GeneAnnotationVersion.objects.filter(pk=pk).exists()


def _human_protein_atlas_finished(pk: int) -> bool:
    return HumanProteinAtlasAnnotationVersion.objects.filter(pk=pk).exists()


# <prefix>_<pk> directory name -> is its owner finished with it? Longest prefix wins when matching
PK_OWNERS = {
    "pipeline": _upload_pipeline_finished,
    "annotation_run": _annotation_run_finished,
    "liftover": _liftover_finished,
    "manual_variants": _manual_variants_finished,
    "classification_import": _classification_import_finished,
    "somalier_vcf_extract": _somalier_vcf_extract_finished,
    "somalier_relate": _somalier_relate_finished,
    "gene_annotation": _gene_annotation_finished,
    "human_protein_atlas": _human_protein_atlas_finished,
}

CLINGEN_PREFIX = "clingen_allele_registry_"
UNIT_TEST_DIR_NAME = "test"
# genes.models.models_gene_coverage nests a dir per collection under here rather than at the top level -
# a sequencing run loads coverage for every sample, so there would be thousands of top level entries
GENE_COVERAGE_DIR_NAME = "gene_coverage"


class Command(BaseCommand):
    help = "Remove import_processing scratch directories whose owner has finished with them"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help="Report what would be removed, removing nothing")
        parser.add_argument('--min-age-days', type=float, default=1,
                            help="Leave entries modified more recently than this alone (default: 1)")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        cutoff = time.time() - options["min_age_days"] * DAY_SECS
        import_processing_dir = settings.IMPORT_PROCESSING_DIR
        if not os.path.isdir(import_processing_dir):
            self.stdout.write(f"'{import_processing_dir}' does not exist - nothing to do\n")
            return

        removed = {}
        kept = {}
        for name in sorted(os.listdir(import_processing_dir)):
            path = os.path.join(import_processing_dir, name)
            keep_reason = self._keep_reason(path, name, cutoff)
            self._tally(kept if keep_reason else removed, path, name, keep_reason)
            if keep_reason:
                continue
            if not dry_run:
                if os.path.isdir(path):
                    shutil.rmtree(path, ignore_errors=True)
                else:
                    os.remove(path)

        self._report("Would remove" if dry_run else "Removed", removed)
        self._report("Kept", kept)

    def _keep_reason(self, path: str, name: str, cutoff: float) -> str:
        """ Why this entry stays; empty string means it can go """
        try:
            modified = os.lstat(path).st_mtime
        except OSError:
            return "unreadable"
        if modified > cutoff:
            return "modified recently"

        if name == UNIT_TEST_DIR_NAME:
            return ""  # Unit test scratch - @see annotation.fake_annotation
        if name == GENE_COVERAGE_DIR_NAME:
            return "gene coverage collections clean up their own dirs"
        if name.startswith(CLINGEN_PREFIX):
            # Empty ones were minted by ClinGenAlleleRegistryAPI merely being constructed (#928);
            # a non-empty one holds an API failure dump someone may still want to read
            if os.path.isdir(path) and not os.listdir(path):
                return ""
            return "holds ClinGen API failure dumps"

        prefix, pk = self._split_prefix_pk(name)
        if prefix is None:
            return "unrecognised"
        if PK_OWNERS[prefix](pk):
            return ""
        return "owner not finished"

    @staticmethod
    def _split_prefix_pk(name: str) -> tuple:
        for prefix in sorted(PK_OWNERS, key=len, reverse=True):
            if name.startswith(prefix + "_"):
                pk = name[len(prefix) + 1:]
                if pk.isdigit():
                    return prefix, int(pk)
        return None, None

    def _tally(self, counts: dict, path: str, name: str, keep_reason: str):
        prefix, _pk = self._split_prefix_pk(name)
        label = prefix or (CLINGEN_PREFIX.rstrip("_") if name.startswith(CLINGEN_PREFIX) else name)
        if keep_reason:
            label = f"{label} ({keep_reason})"
        try:
            size = _dir_size(path) if os.path.isdir(path) else os.lstat(path).st_size
        except OSError:
            size = 0
        num, total_size = counts.get(label, (0, 0))
        counts[label] = (num + 1, total_size + size)

    def _report(self, heading: str, counts: dict):
        if not counts:
            return
        self.stdout.write(f"{heading}:\n")
        for label, (num, size) in sorted(counts.items()):
            self.stdout.write(f"  {label}: {num} ({size / 1024 / 1024:.1f} MB)\n")
