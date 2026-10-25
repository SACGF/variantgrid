"""
The per-version partitioned annotation sources: SubVersionPartition is the base for every sub-version
(ClinVarVersion, GeneAnnotationVersion, HumanProteinAtlasAnnotationVersion, VariantAnnotationVersion) and
bumps AnnotationVersion when one is created. VariantAnnotationVersion is one VEP install plus its data files
for a build; VariantAnnotationVersion.latest and promote_to_active are the entry points.
"""

import logging
from collections.abc import Callable
from contextlib import contextmanager
from functools import cached_property
from typing import Optional

from django.apps import apps
from django.conf import settings
from django.db import models, transaction
from django.db.models import OuterRef, Q, QuerySet, Subquery
from django.db.models.deletion import CASCADE
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.utils import timezone

from annotation.enums import (
    FATHMMPrediction,
    MutationAssessorPrediction,
    MutationTasterPrediction,
    Polyphen2Prediction,
    SIFTPrediction,
    VariantAnnotationPipelineType,
)
from annotation.gene_annotation_release_matching import (
    GeneAnnotationReleaseMatch,
    VepGeneSetVersions,
    match_gene_annotation_release,
)
from annotation.pathogenicity_predictions import (
    TOOLS,
    pred_pathogenic_funcs,
    raw_score_pathogenic_funcs,
)
from annotation.vep_columns import visible_columns_for
from annotation.vep_config import VEPConfig
from genes.enums import AnnotationConsortium
from genes.models import (
    GeneAnnotationRelease,
)
from library.django_utils.data_archive_mixin import DataArchiveMixin
from library.django_utils.django_partition import RelatedModelsPartitionModel
from snpdb.models import (
    GenomeBuild,
)


class InvalidAnnotationVersionError(Exception):
    pass


class SubVersionPartition(RelatedModelsPartitionModel):
    RECORDS_FK_FIELD_TO_THIS_MODEL = "version_id"
    PARTITION_LABEL_TEXT = "version"
    _defer_new_sub_version = False
    created = models.DateTimeField(auto_now_add=True)  # Inserted into DB
    annotation_date = models.DateTimeField(auto_now_add=True)  # Date of annotation (what we sort by)

    class Meta:
        abstract = True

    @staticmethod
    @contextmanager
    def defer_new_sub_version():
        """ Create a batch of sub-versions without bumping AnnotationVersion on each save.
            The caller is responsible for creating the AnnotationVersion afterwards. """
        SubVersionPartition._defer_new_sub_version = True
        try:
            yield
        finally:
            SubVersionPartition._defer_new_sub_version = False

    def save(self, *args, **kwargs):
        created = not self.pk
        super().save(*args, **kwargs)
        if created and not SubVersionPartition._defer_new_sub_version:
            genome_build = getattr(self, "genome_build", None)
            # AnnotationVersion has an FK to every sub-version, so importing it here would be a cycle
            annotation_version_model = apps.get_model("annotation", "AnnotationVersion")
            annotation_version_model.new_sub_version(genome_build)

    def __str__(self):
        date_str = self.annotation_date.strftime("%d %B %Y")
        return f"v{self.pk}. ({date_str})"


class VariantAnnotationVersion(DataArchiveMixin, SubVersionPartition):
    REPRESENTATIVE_TRANSCRIPT_ANNOTATION = "annotation_variantannotation"
    TRANSCRIPT_ANNOTATION = "annotation_varianttranscriptannotation"
    VARIANT_GENE_OVERLAP = "annotation_variantgeneoverlap"
    # Returned by suggested_gene_annotation_release_name when the VEP versions don't tell us the gene set
    UNKNOWN_GENE_ANNOTATION_RELEASE_NAME = "TODO_RELEASE_NAME"
    RECORDS_BASE_TABLE_NAMES = [REPRESENTATIVE_TRANSCRIPT_ANNOTATION, TRANSCRIPT_ANNOTATION, VARIANT_GENE_OVERLAP]

    class Status(models.TextChoices):
        NEW = "NEW", "New"
        ACTIVE = "ACTIVE", "Active"
        HISTORICAL = "HISTORICAL", "Historical"

    class GencodeSubset(models.TextChoices):
        BASIC = "basic", "GENCODE Basic"
        PRIMARY = "primary", "GENCODE Primary"

    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    # GeneAnnotationRelease - imported GTF we can use to get gene/transcript versions that match VEP
    gene_annotation_release = models.ForeignKey(GeneAnnotationRelease, null=True, on_delete=CASCADE)
    last_checked_date = models.DateTimeField(null=True)
    # Lifecycle: NEW (being populated by scheduler) → ACTIVE (current working version)
    # → HISTORICAL (kept for old analyses). @see private issue #577
    status = models.CharField(max_length=10, choices=Status.choices, default=Status.NEW)
    # #1679: latched by annotate_variants when the installed VEP no longer matches this version's pins.
    # Set once (first failure reports to Rollbar), cleared only by the recheck_vep_version task. While set
    # the scheduler creates no range locks and the dispatcher launches no VEP runs against this version.
    vep_version_mismatch_reason = models.TextField(null=True, blank=True)
    vep_version_mismatch_date = models.DateTimeField(null=True, blank=True)

    vep = models.IntegerField()  # code version
    vep_cache = models.IntegerField(null=True)  # May need to have diff code/cache versions eg T2T
    columns_version = models.IntegerField(default=1)
    ensembl = models.TextField()
    # can be eg: ensembl=97.378db18 ensembl-variation=97.26a059c ensembl-io=97.dc917e1 ensembl-funcgen=97.24f4d3c
    ensembl_funcgen = models.TextField()
    ensembl_variation = models.TextField()
    ensembl_io = models.TextField()
    thousand_genomes = models.TextField(blank=True, null=True)  # 37/38 only
    cosmic = models.IntegerField(blank=True, null=True)  # 37/38 only
    hgmd = models.TextField(blank=True, null=True)  # 37/38 only
    assembly = models.TextField()
    dbsnp = models.IntegerField(blank=True, null=True)  # 37/38 only
    gencode = models.TextField(blank=True, null=True)  # 37/38 only
    gencode_subset = models.TextField(
        null=True, blank=True, choices=GencodeSubset.choices,
    )  # None = full Ensembl set; ignored on RefSeq VAVs
    genebuild = models.TextField()
    gnomad = models.TextField(blank=True, null=True)  # 37/38 only
    refseq = models.TextField(blank=True)
    regbuild = models.TextField(blank=True, null=True)  # 37/38 only
    sift = models.TextField(blank=True, null=True)  # 37/38 only
    dbnsfp = models.TextField(blank=True, null=True)  # 37/38 only
    denovo_db = models.TextField(blank=True, null=True)  # 37/38 only
    mave_db = models.TextField(blank=True, null=True)  # GRCh38 only - date string parsed from MaveDB filename
    open_targets = models.TextField(blank=True, null=True)  # release parsed from Open Targets filename
    popeve = models.TextField(blank=True, null=True)  # GRCh38 only - date string parsed from popEVE filename
    spliceai = models.TextField(blank=True, null=True)
    distance = models.IntegerField(default=5000)  # VEP --distance parameter

    # --- remaining VEP command line component pins (#462) ------------------
    # Data files and settings that go into the VEP command line but aren't reported in the ##VEP=
    # header. Files whose name carries a version store the parsed version, the rest store the
    # basename, so installing different data registers as a change - as does removing it, which
    # drops that data's columns. All derived by vep_config.vep_component_version_kwargs().
    # NULL on versions created before #462 other than the ones migration 0161 backfilled.
    alphagenome_splicing = models.TextField(null=True, blank=True)  # non-commercial custom VCF - GRCh38 only
    avi = models.TextField(null=True, blank=True)  # AlphaGenome Variant Impact scores - GRCh38 + VEP >= 116 only
    conservation = models.TextField(null=True, blank=True)  # phastCons/phyloP bigwig basenames
    dbscsnv = models.TextField(null=True, blank=True)
    eve = models.TextField(null=True, blank=True)  # GRCh38 + VEP >= 116 only
    fasta = models.TextField(null=True, blank=True)  # whatever _get_vep_fasta passes to --fasta
    gnomad_sv = models.TextField(null=True, blank=True)
    mastermind = models.TextField(null=True, blank=True)
    maxentscan = models.TextField(null=True, blank=True)
    promoter_ai = models.TextField(null=True, blank=True)  # GRCh38 + VEP >= 116 only
    protvar = models.TextField(null=True, blank=True)
    repeat_masker = models.TextField(null=True, blank=True)
    splicevault = models.TextField(null=True, blank=True)  # converted to RefSeq for RefSeq builds - GRCh37/38 only
    topmed = models.TextField(null=True, blank=True)
    transcript_blocklist = models.TextField(null=True, blank=True)  # --transcript_filter
    uk10k = models.TextField(null=True, blank=True)
    utr_annotator = models.TextField(null=True, blank=True)  # UTRAnnotator uORF file - GRCh37/38 only

    pick_order = models.TextField(null=True, blank=True)  # ANNOTATION_VEP_PICK_ORDER
    sift_enabled = models.BooleanField(null=True)  # vep_config "sift" - whether we pass --sift b
    sv_max_size = models.IntegerField(null=True, blank=True)  # ANNOTATION_VEP_SV_MAX_SIZE
    # ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE - a symbolic del/dup/inv shorter than this went to the STANDARD
    # pipeline (@see annotation.annotation_pipeline_routing)
    structural_variant_min_size = models.IntegerField()
    # ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE - an explicit (sequence) variant whose ref or alt is at least this long
    # goes to the STRUCTURAL_VARIANT pipeline. Null (versions made before #1708) routes them all to STANDARD
    explicit_structural_variant_min_size = models.IntegerField(null=True, blank=True)
    sv_overlap_min_fraction = models.FloatField(null=True, blank=True)  # gnomAD SV overlap_cutoff
    vep_args = models.TextField(null=True, blank=True)  # ANNOTATION_VEP_ARGS, space joined

    # Strategy used to map RefSeq <-> Ensembl transcripts when picking dbNSFP
    # per-transcript scores at insert time. Populated on RefSeq pipelines for
    # columns_version >= 4; NULL otherwise (older versions only consumed
    # variant-level dbNSFP fields, so no mapping was needed). Value is the
    # `name` of a RefSeqEnsemblTranscriptResolver implementation, e.g.
    # "dbnsfp_gene". @see annotation/refseq_ensembl_resolver.py
    transcript_resolver = models.TextField(null=True, blank=True)

    # --- backfill-completion flags ---------------------------------------
    # Whether a VAV's partition has the derived columns the query layer optimises on.
    # New VAVs default True because the bulk VEP inserter populates them at insert
    # time. Every ACTIVE VAV was backfilled before the #736 squash; HISTORICAL ones
    # may still be False, so the consuming node falls back (and warns) for analyses
    # left on them. Once no row is False, the field, the fallback branch and this
    # comment block can go.

    # #574 / commit c6f3c4e6f — DamageNode.spliceai filtering. Source:
    # annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:_add_spliceai_max_ds
    backfilled_spliceai_max_ds = models.BooleanField(default=True)

    # columns_version 2 / 4 rollouts — DamageNode.damage_predictions_min uses
    # predictions_num_pathogenic / predictions_num_benign aggregates which
    # default to 0 on the model, so legacy rows look "no pathogenic predictions"
    # until backfilled. Source:
    # annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:_add_pathogenicity_prediction_counts
    backfilled_damage_counts = models.BooleanField(default=True)

    # #579 - DamageNode.ptc_nmd_escaping filters on nmd_escape_status, which is null on
    # rows annotated before the PTC calculation existed. Source:
    # annotation/vcf_files/bulk_vep_vcf_annotation_inserter.py:_add_calculated_ptc
    backfilled_ptc = models.BooleanField(default=True)

    # #1657 - the SV conservation _max columns our pyBigWig stage wrote are off by one base at the
    # start of the span, and a run whose sidecar stage failed left them null. ConservationNode filters
    # on those columns. Source:
    # annotation/sv_conservation.py:sv_conservation_window
    backfilled_sv_conservation = models.BooleanField(default=True)

    # #2426 - UTRAnnotator, AVI and SpliceVault (and AlphaGenome splicing, #2433) were added to versions
    # already in use: their columns on variants annotated before then stay null until
    # `annotation_backfill_columns --mark-backfilled` imports them, one VEP run with all of them.
    # annotation/vep_plugin_backfill.py:stamp_plugin_pins sets it False on the NEW/ACTIVE versions one applies to.
    backfilled_utr_avi_splicevault = models.BooleanField(default=True)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["genome_build"],
                condition=Q(status="ACTIVE"),
                name="one_active_vav_per_build",
            ),
        ]

    @property
    def gnomad_major_version(self) -> int:
        # Pre-migration GRCh37 versions stored the VEP-style "r2.1"
        major = self.gnomad.split(".", maxsplit=1)[0]
        return int(major.removeprefix("r"))

    @property
    def uses_raw_spliceai(self) -> bool:
        # SpliceAI ships "raw" and "masked" precomputed score files; the masked
        # files zero out splice-site changes Illumina's README flags as typically
        # less pathogenic (strengthening annotated / weakening unannotated) and
        # are the version recommended for variant interpretation. The `spliceai`
        # VEP-pin string is set to e.g. "raw 1.3" or "masked 1.3" by
        # `manage.py vep_install` (see annotation_pipeline_data). See
        # SACGF/variantgrid_sapath#410.
        return bool(self.spliceai) and self.spliceai.strip().lower().startswith("raw")

    @property
    def is_active(self) -> bool:
        return self.status == VariantAnnotationVersion.Status.ACTIVE

    @staticmethod
    def latest(genome_build, status: 'VariantAnnotationVersion.Status' = None) -> 'VariantAnnotationVersion':
        if status is None:
            status = VariantAnnotationVersion.Status.ACTIVE
        qs = VariantAnnotationVersion.objects.filter(genome_build=genome_build, status=status)
        return qs.order_by("annotation_date").last()

    @staticmethod
    def latest_pipeline_routing_sizes(genome_build) -> tuple[int, Optional[int]]:
        """ (structural_variant_min_size, explicit_structural_variant_min_size) for routing variants not
            annotated yet - the setting, until a build has an active version """
        if vav := VariantAnnotationVersion.latest(genome_build):
            return vav.structural_variant_min_size, vav.explicit_structural_variant_min_size
        return settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE, settings.ANNOTATION_STRUCTURAL_VARIANT_MIN_SIZE

    @staticmethod
    def latest_gene_annotation_releases(status: 'VariantAnnotationVersion.Status' = None) -> list[GeneAnnotationRelease]:
        """ The GeneAnnotationRelease of each annotated build's latest VAV (by build name), for builds that have one """
        gene_annotation_releases = []
        for genome_build in GenomeBuild.builds_with_annotation().order_by("name"):
            if vav := VariantAnnotationVersion.latest(genome_build, status=status):
                if vav.gene_annotation_release:
                    gene_annotation_releases.append(vav.gene_annotation_release)
        return gene_annotation_releases

    @staticmethod
    def latest_for_all_builds(status: 'VariantAnnotationVersion.Status' = None) -> QuerySet:
        if status is None:
            status = VariantAnnotationVersion.Status.ACTIVE
        latest_qs = VariantAnnotationVersion.objects.filter(
            genome_build=OuterRef('genome_build'),
            status=status,
        ).order_by('-annotation_date')

        return VariantAnnotationVersion.objects.filter(
            pk__in=Subquery(latest_qs.values('pk')[:1]),
            genome_build__in=GenomeBuild.builds_with_annotation()
        )

    def get_gene_annotation_promote_blocker(self) -> Optional[str]:
        """ Reason gene annotation prevents promoting this VAV to ACTIVE, or None if OK. Activating a VAV whose
            AnnotationVersion is inconsistent (most commonly: gene annotation not run for the current OntologyVersion,
            e.g. ontology was re-imported after this VAV was created) breaks variant pages. The annotation page uses
            this to disable the Promote button and explain why; promote_to_active() enforces it. """
        if not (settings.ANNOTATION_GENE_ANNOTATION_VERSION_ENABLED and self.gene_annotation_release_id):
            return None
        av = self.annotationversion_set.order_by("annotation_date").last()
        if av is None:
            return None
        remediation = "Run 'python3 manage.py gene_annotation --new-releases'."
        if av.gene_annotation_version_id is None:
            return (f"No gene annotation has been created for this GeneAnnotationRelease / current "
                    f"OntologyVersion. {remediation}")
        try:
            av.validate_gene_annotation()
        except InvalidAnnotationVersionError as e:
            return f"{e} {remediation}"
        return None

    def get_annotation_run_blocker(self) -> Optional[str]:
        """ Reason annotation runs against this VAV will fail at pipeline time, or None if OK. Without a
            GeneAnnotationRelease we can't get transcript/versions into the DB for FK linking, and a
            latched VEP mismatch (#1679) fails every run at the version check, so every AnnotationRun
            errors (see annotation.tasks.annotate_variants). While set, the scheduler creates no range
            locks and the dispatcher's VEP lane launches nothing against this version; it is also shown
            as a warning on the annotation runs page. """
        if self.gene_annotation_release_id is None:
            return "missing GeneAnnotationRelease"
        if self.vep_version_mismatch_date:
            mismatch_date = timezone.localtime(self.vep_version_mismatch_date)
            return f"VEP mismatch since {mismatch_date:%Y-%m-%d %H:%M}: {self.vep_version_mismatch_reason}"
        return None

    def latch_vep_version_mismatch(self, reason: str) -> bool:
        """ #1679: latch a VEP mismatch. Returns whether this call made the transition - runs leased
            before the latch fail in the same wave on parallel workers, and the conditional update makes
            exactly one of them the reporter. """
        now = timezone.now()
        updated = VariantAnnotationVersion.objects.filter(
            pk=self.pk, vep_version_mismatch_date__isnull=True).update(vep_version_mismatch_reason=reason,
                                                                       vep_version_mismatch_date=now)
        if updated:
            self.vep_version_mismatch_reason = reason
            self.vep_version_mismatch_date = now
        return bool(updated)

    def clear_vep_version_mismatch(self):
        VariantAnnotationVersion.objects.filter(pk=self.pk).update(vep_version_mismatch_reason=None,
                                                                   vep_version_mismatch_date=None)
        self.vep_version_mismatch_reason = None
        self.vep_version_mismatch_date = None

    @transaction.atomic
    def promote_to_active(self):
        """ Demote prior ACTIVE for this genome_build → HISTORICAL, set self → ACTIVE. """
        if self.status == VariantAnnotationVersion.Status.ACTIVE:
            raise ValueError(f"VariantAnnotationVersion pk={self.pk} is already ACTIVE")
        if self.status == VariantAnnotationVersion.Status.HISTORICAL:
            raise ValueError(f"VariantAnnotationVersion pk={self.pk} is HISTORICAL — cannot promote")

        if blocker := self.get_gene_annotation_promote_blocker():
            raise InvalidAnnotationVersionError(
                f"Cannot promote VariantAnnotationVersion pk={self.pk} to ACTIVE - {blocker}")

        prior_qs = VariantAnnotationVersion.objects.select_for_update().filter(
            genome_build=self.genome_build,
            status=VariantAnnotationVersion.Status.ACTIVE,
        ).exclude(pk=self.pk)
        prior_pks = list(prior_qs.values_list("pk", flat=True))
        prior_qs.update(status=VariantAnnotationVersion.Status.HISTORICAL)
        # The dispatcher never leases runs on a HISTORICAL version, so theirs would sit incomplete forever (#1940)
        annotation_run_model = apps.get_model("annotation", "AnnotationRun")
        num_cancelled = annotation_run_model.cancel_incomplete_runs(
            prior_pks, f"Cancelled: VariantAnnotationVersion went HISTORICAL when pk={self.pk} was promoted")

        self.status = VariantAnnotationVersion.Status.ACTIVE
        self.save(update_fields=["status"])
        logging.info("Promoted VariantAnnotationVersion pk=%s (build=%s) to ACTIVE; demoted prior ACTIVE pks=%s to "
                     "HISTORICAL, cancelling %d incomplete annotation runs", self.pk, self.genome_build, prior_pks,
                     num_cancelled)

    def get_any_annotation_version(self):
        """ Often you don't care what annotation version you use, only that variant annotation version is this one """
        return self.annotationversion_set.last()

    def get_rankscore_pathogenic_prediction_funcs(self) -> dict[str, Callable]:
        if self.columns_version == 1:
            return {
                'sift': lambda d: d in SIFTPrediction.get_damage_or_greater_levels(),
                'fathmm_pred_most_damaging': lambda d: d in FATHMMPrediction.get_damage_or_greater_levels(),
                'mutation_assessor_pred_most_damaging': lambda d: d in MutationAssessorPrediction.get_damage_or_greater_levels(),
                'mutation_taster_pred_most_damaging': lambda d: d in MutationTasterPrediction.get_damage_or_greater_levels(),
                'polyphen2_hvar_pred_most_damaging': lambda d: d in Polyphen2Prediction.get_damage_or_greater_levels(),
            }
        if self.columns_version >= 2:
            pathogenic_rankscore = settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE
            version_columns = visible_columns_for(columns_version=self.columns_version)
            return {t.rankscore_field: lambda d: float(d) >= pathogenic_rankscore
                    for t in TOOLS if t.rankscore_field in version_columns}
        raise ValueError(f"Don't know fields for {self.columns_version=}")

    def get_raw_score_pathogenic_prediction_funcs(self) -> dict[str, Callable]:
        """Raw-score + pred contributions to predictions_num_pathogenic at v4. Empty pre-v4."""
        if self.columns_version < 4:
            return {}
        return {**raw_score_pathogenic_funcs(), **pred_pathogenic_funcs()}

    @cached_property
    def damage_predictions_description(self) -> str:
        if self.columns_version == 1:
            columns = ", ".join(self.get_rankscore_pathogenic_prediction_funcs())
            return f"Count of {columns} at the most damaging level."
        if self.columns_version >= 4:
            raw_columns = ", ".join(self.get_raw_score_pathogenic_prediction_funcs())
            return (f"Count of {raw_columns} reaching their PP3-supporting threshold "
                    "(Pejaver 2022 / Bergquist 2024 ClinGen calibration; see annotation/pathogenicity_predictions.py).")
        rankscore_columns = ", ".join(self.get_rankscore_pathogenic_prediction_funcs())
        return f"Count of {rankscore_columns} that exceed {settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE}"

    @cached_property
    def _vep_config(self) -> dict:
        return self.genome_build.settings["vep_config"]

    @cached_property
    def has_phastcons_30_way_mammalian(self) -> bool:
        return self._vep_config.get("phastcons30way")

    @cached_property
    def has_phylop_30_way_mammalian(self) -> bool:
        return self._vep_config.get("phylop30way")

    @cached_property
    def has_phastcons_46_way_mammalian(self) -> bool:
        return self._vep_config.get("phastcons46way")

    @cached_property
    def has_phylop_46_way_mammalian(self) -> bool:
        return self._vep_config.get("phylop46way")

    def get_visible_columns(self, pipeline_type: Optional[VariantAnnotationPipelineType] = None) -> frozenset[str]:
        """ VariantGrid columns this version populates, off the same VEP_COLUMNS table that controls what
            annotation writes (#1148). Pass a pipeline_type to narrow to what that pipeline wrote. """
        return visible_columns_for(
            vep_config=VEPConfig(self.genome_build),
            genome_build_name=self.genome_build.name,
            pipeline_type=pipeline_type,
            columns_version=self.columns_version,
            vep_version=self.vep,
            cosmic_version=self.cosmic,
            gnomad4_minor_version=self.gnomad,
        )

    @cached_property
    def visible_columns(self) -> frozenset[str]:
        """ Columns populated by any pipeline on this version - what filter nodes can offer, as they
            work off the version rather than an individual annotation run. """
        return self.get_visible_columns()

    @cached_property
    def vep_gene_set_versions(self) -> VepGeneSetVersions:
        """ The VEP version strings describing which gene set was annotated against """
        return VepGeneSetVersions.from_variant_annotation_version(self)

    @property
    def cdot_gene_release_token(self) -> Optional[str]:
        """ The release identifier cdot uses in its per-GFF asset names, eg 'RS_2025_08' """
        return self.vep_gene_set_versions.release_token

    @property
    def suggested_gene_annotation_release_name(self) -> str:
        if release_token := self.cdot_gene_release_token:
            return f"{self.get_annotation_consortium_display()}_{release_token}"
        return self.UNKNOWN_GENE_ANNOTATION_RELEASE_NAME

    def match_gene_annotation_release(self) -> GeneAnnotationReleaseMatch:
        """ Existing GeneAnnotationReleases for our build/consortium built from our VEP's gene set """
        release_urls_by_id = dict(
            GeneAnnotationRelease.objects.filter(genome_build=self.genome_build,
                                                 annotation_consortium=self.annotation_consortium)
            .values_list("pk", "gene_annotation_import__url"))
        return match_gene_annotation_release(release_urls_by_id, self.vep_gene_set_versions)

    def link_gene_annotation_release(self) -> Optional[GeneAnnotationRelease]:
        """ Set gene_annotation_release from an existing matching release. Returns what it linked,
            None if there was no single match (a build's annotation often doesn't change between VEP
            versions, so the release the previous version used is usually still the right one) """
        match = self.match_gene_annotation_release()
        if release_id := match.unique_release_id:
            self.gene_annotation_release = GeneAnnotationRelease.objects.get(pk=release_id)
            self.save()
            logging.info("%s: linked GeneAnnotationRelease '%s' - matched on %s",
                         self, self.gene_annotation_release, match.description)
            return self.gene_annotation_release

        if match.release_ids:
            releases = GeneAnnotationRelease.objects.filter(pk__in=match.release_ids)
            logging.warning("%s: %d GeneAnnotationReleases match on %s (%s) - leaving unlinked",
                            self, len(match.release_ids), match.description,
                            ", ".join(str(r) for r in releases))
        return None

    def short_string(self) -> str:
        """ Short VEP description """
        return f"v{self.pk} {self.vep} {self.get_annotation_consortium_display()}"

    def __str__(self):
        super_str = super().__str__()
        return f"{super_str} VEP: {self.vep} / {self.get_annotation_consortium_display()}/{self.genome_build_id}"


@receiver(pre_delete, sender=VariantAnnotationVersion)
def variant_annotation_version_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    instance.delete_related_objects()
