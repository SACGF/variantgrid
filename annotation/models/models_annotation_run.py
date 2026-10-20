"""
One batch of variants through an annotation pipeline: AnnotationRangeLock (a variant pk range per
VariantAnnotationVersion) and AnnotationRun (one per pipeline type over that range, status derived from its
timestamps), plus AnnotationPipelineVersion for the non-VEP tools and VCFAnnotationStats. Scheduling and the
lease are in annotation/tasks/annotation_scheduler_task.py; AnnotationRun.get_for_variant finds a variant's run.
"""

import os
from datetime import date, datetime, time
from typing import Optional

from django.conf import settings
from django.db import models, transaction
from django.db.models import Q
from django.db.models.deletion import CASCADE, PROTECT
from django.urls import reverse
from django.utils import timezone
from django.utils.text import slugify
from django_extensions.db.models import TimeStampedModel

from annotation.annotation_pipeline_routing import pipeline_type_for_variant
from annotation.annotation_run_files import ANNOTATION_RUN_IMPORT_PROCESSING_PREFIX
from annotation.enums.annotation_enums import (
    AnnotationStatus,
    VariantAnnotationPipelineType,
)
from annotation.models.models import (
    VariantAnnotation,
    VariantGeneOverlap,
    VariantTranscriptAnnotation,
)
from annotation.models.models_version import VariantAnnotationVersion
from library.django_utils.django_file_utils import remove_import_processing_dir
from snpdb.models import (
    VCF,
    GenomeBuild,
    Variant,
)


class VCFAnnotationStats(models.Model):
    """ This is produced in calculate_sample_stats """
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)
    vep_skipped_count = models.IntegerField()

    class Meta:
        unique_together = ('vcf', 'variant_annotation_version')


class AnnotationPipelineVersion(TimeStampedModel):
    """ #720: version of a non-VEP pipeline's tool + data bundle, and what its AnnotationRuns are
        scheduled against.

        Deliberately separate from VariantAnnotationVersion. Pinning eg AnnotSV on the VEP version is what
        makes rolling its bundle cost a full re-annotation of everything - the thing splitting these
        pipelines out exists to avoid. VEP pipelines don't use this: their version IS the
        VariantAnnotationVersion.

        Same lifecycle as VariantAnnotationVersion, for the same reason: upgrading a tool creates a NEW
        row, and promoting it to ACTIVE is what makes the scheduler create runs against it - so a bundle
        swapped on disk never silently re-annotates the genome. The scheduler then backfills through the
        machinery it already has, since "this lock has no run at the ACTIVE version" is the same query as
        "this lock has no run of this type at all".

        HISTORICAL means less here than it does for a VariantAnnotationVersion, which keeps its own
        partition and so keeps its data. A supplementary pipeline writes into columns of the VEP
        partition, one set per VariantAnnotationVersion, so a new version overwrites its predecessor's
        values: HISTORICAL is provenance on the runs it produced, not data still queryable. """

    class Status(models.TextChoices):
        NEW = "NEW", "New"
        ACTIVE = "ACTIVE", "Active"
        HISTORICAL = "HISTORICAL", "Historical"

    pipeline_type = models.CharField(max_length=1, choices=VariantAnnotationPipelineType.choices)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    code_version = models.TextField()                       # eg AnnotSV "3.5.8", from the binary itself
    data_version = models.TextField(null=True, blank=True)  # eg the AnnotSV bundle release string
    status = models.CharField(max_length=10, choices=Status.choices, default=Status.NEW)

    class Meta:
        unique_together = ('pipeline_type', 'genome_build', 'code_version', 'data_version')

    @staticmethod
    def get_active(pipeline_type, genome_build) -> Optional['AnnotationPipelineVersion']:
        """ The version the scheduler creates runs against, or None if the tool has never been registered
            (see the create_new_annotation_pipeline_version command). """
        return AnnotationPipelineVersion.objects.filter(
            pipeline_type=pipeline_type, genome_build=genome_build,
            status=AnnotationPipelineVersion.Status.ACTIVE).first()

    @staticmethod
    def get_new(pipeline_type, genome_build) -> Optional['AnnotationPipelineVersion']:
        """ A registered upgrade awaiting promotion. Latest first - re-registering twice before anyone
            promotes leaves the most recently created one as the candidate. """
        return AnnotationPipelineVersion.objects.filter(
            pipeline_type=pipeline_type, genome_build=genome_build,
            status=AnnotationPipelineVersion.Status.NEW).order_by("-created").first()

    def promote_to_active(self):
        """ Make this the version the scheduler schedules against, retiring the one it replaces. """
        with transaction.atomic():
            AnnotationPipelineVersion.objects.filter(
                pipeline_type=self.pipeline_type, genome_build=self.genome_build,
                status=AnnotationPipelineVersion.Status.ACTIVE,
            ).exclude(pk=self.pk).update(status=AnnotationPipelineVersion.Status.HISTORICAL)
            self.status = AnnotationPipelineVersion.Status.ACTIVE
            self.save()

    @property
    def is_active(self) -> bool:
        return self.status == AnnotationPipelineVersion.Status.ACTIVE

    def __str__(self):
        data = f"/{self.data_version}" if self.data_version else ""
        return f"{self.get_pipeline_type_display()} {self.code_version}{data} ({self.genome_build})"


class AnnotationRangeLock(models.Model):
    MIN_SIZE_FOR_SUBDIVISION = 1000  # Do we need this?

    version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)
    min_variant = models.ForeignKey(Variant, related_name='min_variant', on_delete=PROTECT)
    max_variant = models.ForeignKey(Variant, related_name='max_variant', on_delete=PROTECT)
    count = models.IntegerField(null=True)

    def __str__(self):
        min_v = self.min_variant_id
        max_v = self.max_variant_id
        s = f"AnnotationRangeLock: (v. {self.version}) {min_v} - {max_v}"
        if self.count is not None:
            s += f" (count={self.count})"
        return s

    @staticmethod
    def release_variant(variant: Variant):
        """ Sometimes we want to delete a variant that is either min/max of an AnnotationRangeLock
            This moves the min/max ID in so we can delete it.
        """

        if arl := AnnotationRangeLock.objects.filter(min_variant=variant).first():
            if arl.max_variant == variant:
                # Just contained this variant
                arl.delete()
            else:
                q_contig = Variant.get_contigs_q(arl.version.genome_build)
                if new_min := Variant.objects.filter(q_contig, pk__gt=arl.min_variant_id, pk__lte=arl.max_variant_id).first():
                    arl.min_variant = new_min
                    arl.save()
        elif arl := AnnotationRangeLock.objects.filter(max_variant=variant).first():
            q_contig = Variant.get_contigs_q(arl.version.genome_build)
            if new_max := Variant.objects.filter(q_contig, pk__gte=arl.min_variant_id, pk__lt=arl.max_variant_id).first():
                arl.max_variant = new_max
                arl.save()

    def can_subdivide(self) -> bool:
        """ This is approx (as could be interleaved w/structural etc) """
        return int(self.max_variant_id) - int(self.min_variant_id) >= self.MIN_SIZE_FOR_SUBDIVISION


class AnnotationRun(TimeStampedModel):
    status = models.CharField(max_length=1, choices=AnnotationStatus.choices, default=AnnotationStatus.CREATED)
    annotation_range_lock = models.ForeignKey(AnnotationRangeLock, null=True, on_delete=CASCADE)
    pipeline_type = models.CharField(max_length=1, choices=VariantAnnotationPipelineType.choices,
                                     default=VariantAnnotationPipelineType.STANDARD)
    # task_id is used as a lock to prevent multiple Celery jobs from executing same job
    task_id = models.CharField(max_length=36, null=True)
    # Lease fields (#2667) - mirror analysis NodeTask. The dispatcher (dispatch_annotation_runs)
    # is the single authority that leases a run before launching it; lease_expires lets a dead
    # worker's run be reclaimed. task_id remains the in-annotate_variants execution lock.
    leased_by = models.CharField(max_length=64, null=True)  # worker/dispatch id holding the lease
    lease_expires = models.DateTimeField(null=True)  # for dead-worker reclaim
    attempt_count = models.IntegerField(default=0)  # executions that reached a worker; bounds retries
    # dispatches (lease + task launch). Unbounded and diagnostic only: dispatch_count pulling far ahead
    # of attempt_count means tasks are being queued but never consumed - a broken or starved queue.
    dispatch_count = models.IntegerField(default=0)
    # External annotation (#1568): set by the annotation_external --dump command. The normal scheduler /
    # annotate_variants skip these so VEP is never auto-run on a run the operator is managing externally.
    external = models.BooleanField(default=False)
    # #1657: our pyBigWig stage, not VEP, wrote this run's SV conservation (phastCons/phyloP) _max
    # columns. Recorded rather than inferred from the conservation --custom args in pipeline_command,
    # which an external run doesn't have at all - "which code produced these values" is a property of
    # the run worth keeping once the command that shows it is gone.
    sv_conservation_pybigwig = models.BooleanField(default=False)
    dump_start = models.DateTimeField(null=True)
    dump_end = models.DateTimeField(null=True)
    annotation_start = models.DateTimeField(null=True)
    annotation_end = models.DateTimeField(null=True)
    upload_start = models.DateTimeField(null=True)
    upload_end = models.DateTimeField(null=True)
    upload_attempts = models.IntegerField(default=0)

    pipeline_command = models.TextField(null=True)
    pipeline_stdout = models.TextField(null=True)
    pipeline_stderr = models.TextField(null=True)
    error_exception = models.TextField(null=True)
    vep_skipped_count = models.IntegerField(null=True)
    vep_warnings = models.TextField(null=True)
    vcf_dump_filename = models.TextField(null=True)
    vcf_annotated_filename = models.TextField(null=True)
    # #1701: VEP's --skipped_variants_file for this attempt - one row per record VEP deliberately dropped,
    # which the import lane counts against dump_count to catch a silently truncated annotated VCF
    vep_skipped_variants_filename = models.TextField(null=True)
    # #1646: variants to process, pre-counted off-thread by count_annotation_runs (null until counted;
    # 0 finishes an empty run without a dump). Reset to null when a merge grows the range lock.
    count = models.IntegerField(null=True)
    dump_count = models.IntegerField(null=True)
    annotated_count = models.IntegerField(null=True)
    celery_task_logs = models.JSONField(null=False, default=dict)  # Key=task_id, so we keep logs from multiple runs

    # #720: which version of a non-VEP tool this run is scheduled against. NULL for VEP pipelines, whose
    # version is the VariantAnnotationVersion on the range lock.
    pipeline_version = models.ForeignKey(AnnotationPipelineVersion, null=True, blank=True, on_delete=PROTECT)

    class Meta:
        # #720: one run per lock per pipeline type - except for versioned pipelines, which get one per
        # version, because that is how they re-run. A new VEP version brings its own range locks, so VEP
        # never needed this; a supplementary pipeline hangs off the locks the VEP run already has, and a
        # tool upgrade has to produce a second run over the same lock. Split in two rather than one
        # 3-column constraint because NULLs compare distinct in Postgres, which would drop the VEP
        # guarantee entirely.
        constraints = [
            models.UniqueConstraint(fields=["annotation_range_lock", "pipeline_type"],
                                    condition=Q(pipeline_version__isnull=True),
                                    name="one_unversioned_run_per_range_lock_and_type"),
            models.UniqueConstraint(fields=["annotation_range_lock", "pipeline_type", "pipeline_version"],
                                    condition=Q(pipeline_version__isnull=False),
                                    name="one_versioned_run_per_range_lock_type_and_version"),
        ]

    def get_absolute_url(self):
        return reverse('view_annotation_run', kwargs={'annotation_run_id': self.pk})

    @staticmethod
    def count_not_successful_runs_for_version(version: VariantAnnotationVersion) -> int:
        qs = AnnotationRun.objects.filter(annotation_range_lock__version=version)
        return qs.exclude(status=AnnotationStatus.FINISHED).count()

    @staticmethod
    def cancel_incomplete_runs(variant_annotation_version_ids, reason: str) -> int:
        """ ERROR the incomplete runs on these versions, for a version the dispatcher no longer leases
            (HISTORICAL). A run a worker holds a live lease on is left to finish or fail on its own: the
            worker's next save would overwrite the error with its in-memory state. """
        qs = AnnotationRun.objects.filter(annotation_range_lock__version__in=variant_annotation_version_ids)
        qs = qs.exclude(status__in=AnnotationStatus.get_completed_states())
        qs = qs.exclude(lease_expires__gte=timezone.now())
        # update() skips save(), so set the status get_status() derives from error_exception
        return qs.update(error_exception=reason, status=AnnotationStatus.ERROR)

    @staticmethod
    def get_for_variant(variant: Variant, genome_build) -> Optional['AnnotationRun']:
        # For newly created variants, there will only be one per build for latest annotation version.
        # Each version routes by its own cut-off, so the pipeline is decided per run
        qs = AnnotationRun.objects.filter(annotation_range_lock__version__genome_build=genome_build,
                                          annotation_range_lock__min_variant__lte=variant.pk,
                                          annotation_range_lock__max_variant__gte=variant.pk)
        for ar in qs.select_related("annotation_range_lock__version").order_by("pk"):
            sv_min_size = ar.annotation_range_lock.version.structural_variant_min_size
            if ar.pipeline_type == pipeline_type_for_variant(variant, sv_min_size):
                return ar
        return None

    @staticmethod
    def get_active_runs(genome_build):
        qs = AnnotationRun.objects.filter(annotation_range_lock__version__genome_build=genome_build)
        return qs.exclude(status__in=AnnotationStatus.get_completed_states())

    @property
    def variant_annotation_version(self):
        return self.annotation_range_lock.version

    def save(self, *args, **kwargs):
        self.status = self.get_status()
        super().save(*args, **kwargs)

    def get_status(self):
        status = AnnotationStatus.CREATED
        if self.error_exception:
            status = AnnotationStatus.ERROR
        else:
            if self.dump_start:
                status = AnnotationStatus.DUMP_STARTED
            if self.dump_end:
                status = AnnotationStatus.DUMP_COMPLETED
            if self.dump_count == 0:
                status = AnnotationStatus.FINISHED
            elif self.external and self.vcf_annotated_filename is None:
                # External annotation (#1568): an operator-managed run parks awaiting external annotation for
                # its whole life (from being claimed external, through the local dump, until --import sets
                # vcf_annotated_filename and the normal flow resumes) rather than surfacing the transient
                # local dump states (DUMP_STARTED/DUMP_COMPLETED). This keeps `status` in sync with the
                # `external` flag - the two never disagree - so a run whose dump is in progress or died
                # mid-dump still reads "Awaiting external annotation" instead of getting stuck at DUMP_STARTED.
                status = AnnotationStatus.EXTERNAL_DUMP_COMPLETED
            else:
                if self.annotation_start:
                    status = AnnotationStatus.ANNOTATION_STARTED
                if self.annotation_end:
                    status = AnnotationStatus.ANNOTATION_COMPLETED
                if self.upload_start:
                    status = AnnotationStatus.UPLOAD_STARTED
                if self.upload_end:
                    status = AnnotationStatus.FINISHED
        return status

    def _lease_is_live(self, now=None) -> bool:
        if self.lease_expires is None:
            return False
        if now is None:
            now = timezone.now()
        return self.lease_expires >= now

    def is_dispatchable(self, now=None) -> bool:
        """ #2667: Ready for the dispatcher to lease + launch - pending state with no live lease.
            External runs (#1568) are operator-managed and never auto-dispatched.

            #1646: strictly CREATED - a resume-upload run (is_upload_resumable) is dispatched via a
            separate lane. Keeping this CREATED-only means merge (_range_lock_is_dispatchable) never
            grows a lock whose run already owns committed VEP output for its exact range. """
        if self.external:
            return False
        if self.task_id is not None:
            return False
        if self.status != AnnotationStatus.CREATED:
            return False
        return not self._lease_is_live(now)

    @property
    def is_empty_finished(self) -> bool:
        """ #1646: a run the count task found empty (no variants of its pipeline_type in range) and
            finished without a dump. Owns no annotation rows, so a merge that grows its range lock can
            safely reopen it (_absorb_range_lock) for a re-count over the larger range. """
        return self.status == AnnotationStatus.FINISHED and self.dump_count == 0

    def reopen_to_created(self):
        """ #1646: return an empty-finished run to an un-counted CREATED state so it is re-counted over
            a now-larger range (a merge grew its lock). Only valid for empty runs - they own no rows. """
        self.dump_start = None
        self.dump_end = None
        self.dump_count = None
        self.count = None
        self.save()  # get_status() -> CREATED

    def reset_for_retry(self):
        """ #1654: return this run to a clean, un-counted CREATED state so the dispatcher re-counts +
            re-annotates it from scratch, keeping its existing annotation_range_lock. Reusing the run
            in place (rather than delete + recreate) means no rangeless AnnotationRun is ever committed
            and there's no unique_together (annotation_range_lock, pipeline_type) window to work around.

            Clears any partially-imported annotation rows first - which can take a minute or two for a
            large range - so callers run this on a worker rather than in a request. The flip to CREATED
            happens only after that delete, so the dispatcher never re-annotates over stale rows; an
            interrupted delete leaves the run in its current (visible, retryable) state, never orphaned. """
        with transaction.atomic():
            self.delete_related_objects()
            # Overwrite the row in place with a fresh run carrying the same pk + range lock, so every other
            # field returns to its default (no field-by-field reset to keep in sync as the model grows).
            # pk is already set -> save() does an UPDATE of all columns; get_status() -> CREATED.
            # pipeline_version (#720) is part of what this run was scheduled against, like its lock and
            # type - a retry has to stay on the same version or check_tool_version has nothing to compare
            # the installed tool against.
            fresh = AnnotationRun(annotation_range_lock=self.annotation_range_lock,
                                  pipeline_type=self.pipeline_type,
                                  pipeline_version=self.pipeline_version)
            fresh.pk = self.pk
            fresh.created = self.created  # auto_now_add isn't re-applied on the UPDATE below
            fresh.save()
        self.refresh_from_db()  # keep this instance consistent with the reset row
        # #1596: the DB rows are gone, so drop the matching on-disk import-processing scratch dir too.
        # Otherwise the previous attempt's TSVs survive the reset, and because upload_attempts is back to 0
        # the upload_attempts>1 cleanup in import_vcf_annotations is skipped - so the leftover files trip
        # write_sql_copy_csv's "don't want to overwrite" guard, which is meant only for genuinely out-of-sync
        # dirs (moved dump / double launch). Done outside the transaction (filesystem op) and after the DB
        # reset commits, so a rolled-back reset leaves the scratch dir intact.
        remove_import_processing_dir(self.pk, prefix=ANNOTATION_RUN_IMPORT_PROCESSING_PREFIX)

    def revert_external_to_local(self):
        """ #1568: return an external run to the normal local pipeline. Clears the external flag and dump
            state (like reopen_to_created) so get_status() -> CREATED and the dispatcher - which sweeps the
            latest NEW/ACTIVE versions filtering external=False - re-dumps + annotates it locally. Used when
            a dump is too small to be worth the external round-trip (see
            external_annotation.DEFAULT_MIN_EXTERNAL_VARIANTS). """
        self.external = False
        self.dump_start = None
        self.dump_end = None
        self.dump_count = None
        self.count = None
        self.save()  # get_status() -> CREATED

    def is_upload_resumable(self, now=None) -> bool:
        """ #1646: a stalled run that already has an annotated VCF (past VEP) and only needs the quick
            DB upload. reclaim_stalled_annotation_runs keeps the annotated file and scrubs partial upload
            rows, leaving status ANNOTATION_COMPLETED; the dispatcher re-launches it upload-only and
            annotate_variants skips straight to import_vcf_annotations. """
        if self.external:
            return False
        if self.task_id is not None:
            return False
        if self.status != AnnotationStatus.ANNOTATION_COMPLETED:
            return False
        if not self.vcf_annotated_filename:
            return False
        return not self._lease_is_live(now)

    def is_in_flight(self, now=None) -> bool:
        """ #2667: A run that currently occupies a worker slot - live lease, holding the task_id
            execution lock, or part-way through the pipeline (a running, non-completed status). A
            completed run never occupies a slot even if its lease has not been cleared yet. """
        if self.status in AnnotationStatus.get_completed_states():
            return False
        # #1646: a reclaimed resume-upload run waits for the dispatcher - it holds no worker slot.
        if self.is_upload_resumable(now):
            return False
        if self._lease_is_live(now):
            return True
        if self.task_id is not None:
            return True
        return self.status != AnnotationStatus.CREATED

    @property
    def annotation_consortium(self):
        return self.annotation_range_lock.version.annotation_consortium

    @property
    def genome_build(self):
        return self.annotation_range_lock.version.genome_build

    def delete_related_objects(self):
        from annotation.annotation_version_querysets import get_queryset_for_annotation_version

        if self.annotation_range_lock_id is None:
            return  # no version to name a partition by - delete()'s cascade over the base tables covers it
        annotation_version = self.variant_annotation_version.get_any_annotation_version()
        for klass in [VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap]:
            qs = get_queryset_for_annotation_version(klass, annotation_version)
            qs.filter(annotation_run=self).delete()

    PIPELINE_TYPE_DESC = {
        VariantAnnotationPipelineType.STANDARD: "standard",
        VariantAnnotationPipelineType.STRUCTURAL_VARIANT: "structural_variant",
        VariantAnnotationPipelineType.GENE_LEVEL: "gene_level",
        VariantAnnotationPipelineType.ANNOTSV: "annotsv",
    }

    def _get_dump_path_stem(self, dump_dir=None, task_token=None) -> str:
        """ Self-describing dump path without extension (#1568): the stem carries site/build/version/run
            identity so dumps copied to other machines remain self-explanatory. Matching is driven by the
            sidecar metadata (get_dump_metadata), not by parsing this name. dump_dir defaults to
            ANNOTATION_VCF_DUMP_DIR; the annotation_external --dump command passes its --output-dir instead.

            #1658: the local VEP pipeline passes task_token (the executing Celery task_id) so a reclaimed
            run's new attempt writes to a path distinct from the losing zombie worker's - two concurrent
            executions always have distinct task_ids by construction, so each VEP owns its own file and
            neither can corrupt the other. External dumps pass no token, keeping their stable name. """
        type_desc = self.PIPELINE_TYPE_DESC.get(self.pipeline_type, str(self.pipeline_type))
        vav = self.variant_annotation_version
        site = slugify(settings.SITE_NAME) or "site"
        stem = (f"{site}__{vav.genome_build.name}__{vav.get_annotation_consortium_display()}"
                f"__vep{vav.vep}__cv{vav.columns_version}__gar{vav.gene_annotation_release_id}"
                f"__run{self.pk}__{type_desc}")
        if task_token:
            stem += f"__task{slugify(task_token)}"
        return os.path.join(dump_dir or settings.ANNOTATION_VCF_DUMP_DIR, stem)

    def get_dump_filename(self, dump_dir=None, task_token=None) -> str:
        # External dumps (#1568) are gzipped - smaller to ship to the compute box, and VEP reads gzipped
        # input natively. In-VM dumps stay plain (fed straight to VEP then deleted).
        suffix = ".vcf.gz" if self.external else ".vcf"
        return self._get_dump_path_stem(dump_dir=dump_dir, task_token=task_token) + suffix

    def get_dump_metadata_filename(self, dump_dir=None) -> str:
        """ Sidecar metadata path written next to the dump VCF (#1568). """
        return self._get_dump_path_stem(dump_dir=dump_dir) + ".meta.json"

    def delete(self, using=None, keep_parents=False):
        self.delete_related_objects()
        super().delete(using=using, keep_parents=keep_parents)

    def set_task_log(self, key, value):
        """ #1658: append to this task's entry in celery_task_logs (keyed by task_id so a reclaimed run
            keeps every attempt's log - the audit trail for diagnosing a reclaim).

            Values are coerced to JSON-native types here rather than at the call sites: celery_task_logs is
            a plain JSONField with no encoder=, and every real caller logs a timezone.now(). Coercing here
            keeps the field's stored representation unchanged for every other value (adding
            encoder=DjangoJSONEncoder would need a migration and would alter how they all round-trip). """
        assert self.task_id is not None
        if isinstance(value, (datetime, date, time)):
            value = value.isoformat()
        self.celery_task_logs.setdefault(self.task_id, {})[key] = value

    def save_if_owner(self, task_id) -> bool:
        """ #1658: ownership-guarded counterpart to save(), for the writes that happen after a stage this
            worker may have lost the run during (VEP / the bulk import). A stalled attempt whose lease
            expired has its run reclaimed and handed to a fresh attempt, but it still holds the ORM
            instance it loaded before the reclaim - saving that whole instance writes its status, dump
            path and traceback onto a row the new attempt now owns, and blanks the new attempt's lease.

            Issues a single conditional UPDATE of all concrete fields, matched on the task_id we believe we
            hold; returns whether it matched a row. Being one statement it is atomic, so unlike a
            read-then-write ownership check there is no TOCTOU window. status is recomputed here exactly as
            save() does, so it can't desynchronise from the dump_*/annotation_* timestamps. """
        self.status = self.get_status()
        self.modified = timezone.now()  # auto_now isn't applied by .update()
        update_kwargs = {f.attname: getattr(self, f.attname)
                         for f in self._meta.concrete_fields if not f.primary_key}
        return bool(AnnotationRun.objects.filter(pk=self.pk, task_id=task_id).update(**update_kwargs))

    def __str__(self):
        return f"AnnotationRun: {self.pk}/{self.get_pipeline_type_display()}: ({self.get_status_display()})"
