import json
import logging
import os
import re
from collections import defaultdict
from datetime import datetime, timedelta
from typing import List, Optional, Dict

from Bio import Entrez
from Bio.Data.IUPACData import protein_letters_1to3
from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db import models, transaction
from django.db.models.deletion import PROTECT, CASCADE, SET_NULL
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.utils import timezone
from django.utils.timezone import localtime
from django_extensions.db.models import TimeStampedModel
from lazy import lazy

from annotation.external_search_terms import get_variant_search_terms, get_variant_pubmed_search_terms
from annotation.models.damage_enums import Polyphen2Prediction, FATHMMPrediction, MutationTasterPrediction, \
    SIFTPrediction, PathogenicityImpact, MutationAssessorPrediction
from annotation.models.models_enums import AnnotationStatus, CitationSource, \
    VariantClass, ColumnAnnotationCategory, VEPPlugin, VEPCustom, ClinVarReviewStatus, VEPSkippedReason, \
    ManualVariantEntryType, HumanProteinAtlasAbundance, LOFTEEConfidence
from genes.models import GeneSymbol, Gene, TranscriptVersion, Transcript, GeneAnnotationRelease
from genes.models_enums import AnnotationConsortium
from library.django_utils import object_is_referenced
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.utils import invert_dict
from ontology.models import OntologyImport
from patients.models_enums import GnomADPopulation
from snpdb.models import GenomeBuild, Variant, VariantGridColumn, Q, VCF, DBSNP_PATTERN, VARIANT_PATTERN
from snpdb.models.models_enums import ImportStatus


class SubVersionPartition(RelatedModelsPartitionModel):
    RECORDS_FK_FIELD_TO_THIS_MODEL = "version_id"
    PARTITION_LABEL_TEXT = "version"
    created = models.DateTimeField(auto_now_add=True)  # Inserted into DB
    annotation_date = models.DateTimeField(auto_now_add=True)  # Date of annotation (what we sort by)

    class Meta:
        abstract = True

    def save(self, **kwargs):
        created = not self.pk
        super().save(**kwargs)
        if created:
            genome_build = getattr(self, "genome_build", None)
            AnnotationVersion.new_sub_version(genome_build)

    def __str__(self):
        date_str = self.annotation_date.strftime("%d %B %Y")
        return f"v{self.pk}. ({date_str})"


class ClinVarVersion(SubVersionPartition):
    RECORDS_BASE_TABLE_NAMES = ["annotation_clinvar"]
    filename = models.TextField()
    md5_hash = models.CharField(max_length=32)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)

    @staticmethod
    def get_annotation_date_from_filename(filename):
        # file looks like: clinvar_20160302.vcf.gz
        base_name = os.path.basename(filename)
        # Allow a bit of lee-way as may want to call it eg clinvar_20190708_grch37.vcf.gz
        CLINVAR_PATTERN = r"^clinvar_(\d{8}).*\.vcf"
        if m := re.match(CLINVAR_PATTERN, base_name):
            date_time = m.group(1)
            annotation_date = datetime.strptime(date_time, "%Y%m%d")
            return annotation_date
        msg = f"File name '{base_name}' didn't match pattern {CLINVAR_PATTERN}"
        raise ValueError(msg)


@receiver(pre_delete, sender=ClinVarVersion)
def clinvar_version_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    instance.delete_related_objects()


class ClinVar(models.Model):
    ALLELE_ORIGIN = {0: 'unknown',
                     1: 'germline',
                     2: 'somatic',
                     4: 'inherited',
                     8: 'paternal',
                     16: 'maternal',
                     32: 'de-novo',
                     64: 'biparental',
                     128: 'uniparental',
                     256: 'not-tested',
                     512: 'tested-inconclusive',
                     1073741824: 'other'}

    SUSPECT_REASON_CODES = {0: 'unspecified',
                            1: 'Paralog',
                            2: 'byEST',
                            4: 'oldAlign',
                            8: 'Para_EST',
                            16: '1kg_failed',
                            1024: 'other'}

    version = models.ForeignKey(ClinVarVersion, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, on_delete=PROTECT)
    clinvar_variation_id = models.IntegerField()
    clinvar_allele_id = models.IntegerField()
    clinvar_preferred_disease_name = models.TextField(null=True, blank=True)
    clinvar_disease_database_name = models.TextField(null=True, blank=True)
    clinvar_review_status = models.CharField(max_length=1, null=True, choices=ClinVarReviewStatus.choices)
    clinical_significance = models.TextField(null=True, blank=True)
    # If clinical_significance = 'Conflicting_interpretations_of_pathogenicity'
    conflicting_clinical_significance = models.TextField(null=True, blank=True)
    highest_pathogenicity = models.IntegerField(default=0)  # Highest of clinical_significance
    clinvar_clinical_sources = models.TextField(null=True, blank=True)
    clinvar_origin = models.IntegerField(default=0)
    clinvar_suspect_reason_code = models.IntegerField(default=0)
    drug_response = models.BooleanField(default=False)

    @property
    def stars(self):
        return ClinVarReviewStatus(self.clinvar_review_status).stars()

    def get_origin_display(self):
        return ClinVar.ALLELE_ORIGIN.get(self.clinvar_origin)

    def get_suspect_reason_code_display(self):
        return ClinVar.SUSPECT_REASON_CODES.get(self.clinvar_suspect_reason_code)

    def get_citations(self):
        cvc_qs = ClinVarCitation.objects.filter(clinvar_variation_id=self.clinvar_variation_id,
                                                clinvar_allele_id=self.clinvar_allele_id)
        return Citation.objects.filter(clinvarcitation__in=cvc_qs)

    def __str__(self):
        return f"ClinVar: variant: {self.variant}, path: {self.highest_pathogenicity}"


class ClinVarCitationsCollection(models.Model):
    cached_web_resource = models.ForeignKey('CachedWebResource', null=True, on_delete=CASCADE)


class Citation(models.Model):
    citation_source = models.CharField(max_length=1, choices=CitationSource.choices)
    citation_id = models.TextField()

    class Meta:
        unique_together = ("citation_source", "citation_id")

    def unique_code(self):
        return f"{self.citation_source}_{self.citation_id}"

    def __str__(self):
        citation_source = self.get_citation_source_display()
        return f"{self.citation_id} ({citation_source})"

    def ref_id(self):
        citation_source = self.get_citation_source_display()
        return f"{citation_source}:{self.citation_id}"

    @property
    def sort_key(self):
        return self.citation_source, self.citation_id.rjust(10, '0')

    def __lt__(self, other):
        self._sort_key < other._sort_key

    @staticmethod
    def citations_from_text(text):
        """ returns a list of (unsaved) Citation objects from text """
        # TODO replace with code from dbregexes

        citation_source_codes = dict({k.lower(): v for k, v in CitationSource.CODES.items()})
        regex_pattern = r"(%s):\s*(\d+)" % '|'.join(citation_source_codes)
        pattern = re.compile(regex_pattern, flags=re.IGNORECASE)  # @UndefinedVariable

        citations_list = []
        if text:
            # Find all PUBMED
            for m in pattern.finditer(text):
                citation_source_string = m.group(1)
                citation_id = m.group(2)

                citation_source = citation_source_codes.get(citation_source_string.lower())
                if citation_source:
                    citation, _ = Citation.objects.get_or_create(citation_source=citation_source,
                                                                 citation_id=citation_id)
                    citations_list.append(citation)

        return citations_list


class ClinVarCitation(models.Model):
    clinvar_citations_collection = models.ForeignKey(ClinVarCitationsCollection, on_delete=CASCADE)
    clinvar_variation_id = models.IntegerField()
    clinvar_allele_id = models.IntegerField()
    citation = models.ForeignKey(Citation, null=True, on_delete=CASCADE)


class CitationException(Exception):
    pass


class CachedCitation(TimeStampedModel):
    citation = models.OneToOneField(Citation, null=True, on_delete=CASCADE)
    json_string = models.TextField()
    has_error = models.BooleanField(blank=True, null=False, default=False)

    def get_record_or_fail(self):
        if self.has_error:
            raise CitationException(self.json_string)
        return self.get_record()

    def get_record(self):
        """ You can cache this by setting _record """
        record = getattr(self, "_record", None)
        if not record:
            record = json.loads(self.json_string)
            self._record = record
        return record


class GeneAnnotationVersion(SubVersionPartition):
    RECORDS_BASE_TABLE_NAMES = ["annotation_geneannotation"]
    gene_annotation_release = models.ForeignKey(GeneAnnotationRelease, on_delete=CASCADE)
    last_ontology_import = models.ForeignKey(OntologyImport, on_delete=PROTECT)
    gnomad_import_date = models.DateTimeField()

    @property
    def genome_build(self):
        return self.gene_annotation_release.genome_build


@receiver(pre_delete, sender=GeneAnnotationVersion)
def gene_annotation_version_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    instance.delete_related_objects()


class GeneAnnotation(models.Model):
    """ This is generated against genes found via GeneAnnotationRelease
        so that data matches up in analyses """
    version = models.ForeignKey(GeneAnnotationVersion, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    hpo_terms = models.TextField(null=True)
    omim_terms = models.TextField(null=True)
    mondo_terms = models.TextField(null=True)
    gene_disease_moderate_or_above = models.TextField(null=True)
    gene_disease_supportive_or_below = models.TextField(null=True)
    gnomad_oe_lof = models.FloatField(null=True)  # Copied from genes.GnomADGeneConstraint

    class Meta:
        unique_together = ("version", "gene")

    @staticmethod
    def get_for_gene_version(gene_version):
        pass


class HumanProteinAtlasAnnotationVersion(SubVersionPartition):
    RECORDS_BASE_TABLE_NAMES = ["annotation_humanproteinatlasannotation"]
    filename = models.TextField()
    md5_hash = models.CharField(max_length=32)
    hpa_version = models.FloatField()
    unit = models.TextField()  # What unit "value" is in

    def get_minimum_for_abundance_level(self, abundance: HumanProteinAtlasAbundance):
        """ Attempt to remain backwards compatabile - see how this is calculated here:
            https://github.com/SACGF/variantgrid/issues/9#issuecomment-563510678 """
        abundance_mins = {
            HumanProteinAtlasAbundance.NOT_DETECTED: 0,
            HumanProteinAtlasAbundance.LOW: 1.1,
            HumanProteinAtlasAbundance.MEDIUM: 15.5,
            HumanProteinAtlasAbundance.HIGH: 83.6,
        }

        if self.hpa_version == 15:
            abundance_mins = {
                HumanProteinAtlasAbundance.NOT_DETECTED: 0,
                HumanProteinAtlasAbundance.LOW: 0.5,
                HumanProteinAtlasAbundance.MEDIUM: 10,
                HumanProteinAtlasAbundance.HIGH: 50,
            }
        return abundance_mins[abundance]


@receiver(pre_delete, sender=HumanProteinAtlasAnnotationVersion)
def human_protein_annotation_version_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    instance.delete_related_objects()


class HumanProteinAtlasTissueSample(models.Model):
    name = models.TextField(unique=True)

    def __str__(self):
        return self.name


class HumanProteinAtlasAnnotation(models.Model):
    version = models.ForeignKey(HumanProteinAtlasAnnotationVersion, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, null=True, on_delete=CASCADE)  # Always Ensembl
    tissue_sample = models.ForeignKey(HumanProteinAtlasTissueSample, on_delete=CASCADE)
    value = models.FloatField()


class ColumnVEPField(models.Model):
    """ For VariantAnnotation/Transcript columns derived from VEP fields """
    column = models.TextField(unique=True)
    variant_grid_column = models.ForeignKey(VariantGridColumn, on_delete=CASCADE)
    genome_build = models.ForeignKey(GenomeBuild, null=True, on_delete=CASCADE)  # null = all builds
    category = models.CharField(max_length=1, choices=ColumnAnnotationCategory.choices)
    source_field = models.TextField(null=True)  # @see use vep_info_field
    source_field_processing_description = models.TextField(null=True)
    vep_plugin = models.CharField(max_length=1, choices=VEPPlugin.choices, null=True)
    vep_custom = models.CharField(max_length=1, choices=VEPCustom.choices, null=True)
    source_field_has_custom_prefix = models.BooleanField(default=False)
    # We can use these min/max versions to turn on/off columns over time
    min_vep_columns_version = models.IntegerField(null=True)
    max_vep_columns_version = models.IntegerField(null=True)

    @property
    def vep_info_field(self):
        """ For VCFs, be sure to set source_field_has_custom_prefix=True
            Annotating with a VCF (short name = TopMed) brings in the DB column as "TopMed" and
            prefixes INFO fields, eg 'TOPMED' => TopMed_TOPMED.
            We need to adjust for this in BulkVEPVCFAnnotationInserter """

        vif = self.source_field
        if self.source_field_has_custom_prefix:
            vif = self.get_vep_custom_display() + "_" + vif
        return vif

    @staticmethod
    def get_columns_version_q(columns_version: int) -> Q:
        q_min = Q(min_vep_columns_version__isnull=True) | Q(min_vep_columns_version__lte=columns_version)
        q_max = Q(max_vep_columns_version__isnull=True) | Q(max_vep_columns_version__gte=columns_version)
        return q_min & q_max

    @staticmethod
    def filter_for_build(genome_build: GenomeBuild):
        """ genome_build = NULL (no build) or matches provided build """
        return ColumnVEPField.objects.filter(Q(genome_build=genome_build) | Q(genome_build__isnull=True))

    @staticmethod
    def get_source_fields(genome_build: GenomeBuild, *columnvepfield_args, **columnvepfield_kwargs):
        qs = ColumnVEPField.filter_for_build(genome_build)
        qs = qs.filter(*columnvepfield_args, **columnvepfield_kwargs).distinct("source_field")
        return list(qs.values_list("source_field", flat=True).order_by("source_field"))


class VariantAnnotationVersion(SubVersionPartition):
    REPRESENTATIVE_TRANSCRIPT_ANNOTATION = "annotation_variantannotation"
    TRANSCRIPT_ANNOTATION = "annotation_varianttranscriptannotation"
    VARIANT_GENE_OVERLAP = "annotation_variantgeneoverlap"
    RECORDS_BASE_TABLE_NAMES = [REPRESENTATIVE_TRANSCRIPT_ANNOTATION, TRANSCRIPT_ANNOTATION, VARIANT_GENE_OVERLAP]

    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.choices)
    # GeneAnnotationRelease - imported GTF we can use to get gene/transcript versions that match VEP
    gene_annotation_release = models.ForeignKey(GeneAnnotationRelease, null=True, on_delete=CASCADE)
    last_checked_date = models.DateTimeField(null=True)
    active = models.BooleanField(default=True)

    vep = models.IntegerField()
    columns_version = models.IntegerField(default=1)
    ensembl = models.TextField()
    # can be eg: ensembl=97.378db18 ensembl-variation=97.26a059c ensembl-io=97.dc917e1 ensembl-funcgen=97.24f4d3c
    ensembl_funcgen = models.TextField()
    ensembl_variation = models.TextField()
    ensembl_io = models.TextField()
    thousand_genomes = models.TextField()
    cosmic = models.IntegerField()
    hgmd = models.TextField()
    assembly = models.TextField()
    dbsnp = models.IntegerField()
    gencode = models.TextField()
    genebuild = models.TextField()
    gnomad = models.TextField()
    refseq = models.TextField(blank=True)
    regbuild = models.TextField()
    sift = models.TextField()
    dbnsfp = models.TextField()
    distance = models.IntegerField(default=5000)  # VEP --distance parameter

    @staticmethod
    def latest(genome_build, active=True):
        qs = VariantAnnotationVersion.objects.filter(genome_build=genome_build, active=active)
        return qs.order_by("annotation_date").last()

    def get_any_annotation_version(self):
        """ Often you don't care what annotation version you use, only that variant annotation version is this one """
        return self.annotationversion_set.last()

    def get_functional_prediction_pathogenic_levels(self) -> Dict:
        if self.columns_version == 1:
            return {
                'sift': SIFTPrediction.get_damage_or_greater_levels(),
                'fathmm_pred_most_damaging': FATHMMPrediction.get_damage_or_greater_levels(),
                'mutation_assessor_pred_most_damaging': MutationAssessorPrediction.get_damage_or_greater_levels(),
                'mutation_taster_pred_most_damaging': MutationTasterPrediction.get_damage_or_greater_levels(),
                'polyphen2_hvar_pred_most_damaging': Polyphen2Prediction.get_damage_or_greater_levels(),
            }
        elif self.columns_version == 2:
            return {}  # none used
        raise ValueError(f"Don't know fields for {self.columns_version=}")

    @lazy
    def _vep_config(self) -> Dict:
        return self.genome_build.settings["vep_config"]

    @lazy
    def has_phastcons_30_way_mammalian(self) -> bool:
        return self._vep_config.get("phastcons30way")

    @lazy
    def has_phylop_30_way_mammalian(self) -> bool:
        return self._vep_config.get("phylop30way")

    @lazy
    def has_phastcons_46_way_mammalian(self) -> bool:
        return self._vep_config.get("phastcons46way")

    @lazy
    def has_phylop_46_way_mammalian(self) -> bool:
        return self._vep_config.get("phylop46way")

    def short_string(self) -> str:
        """ Short VEP description """
        return f"v{self.pk} {self.vep} {self.get_annotation_consortium_display()}"

    def __str__(self):
        super_str = super().__str__()
        return f"{super_str} VEP: {self.vep} / {self.get_annotation_consortium_display()}/{self.genome_build_id}"


@receiver(pre_delete, sender=VariantAnnotationVersion)
def variant_annotation_version_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    instance.delete_related_objects()


class VCFAnnotationStats(models.Model):
    """ This is produced in calculate_sample_stats """
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)
    vep_skipped_count = models.IntegerField()

    class Meta:
        unique_together = ('vcf', 'variant_annotation_version')


class AnnotationRangeLock(models.Model):
    version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)
    min_variant = models.ForeignKey(Variant, related_name='min_variant', on_delete=PROTECT)
    max_variant = models.ForeignKey(Variant, related_name='max_variant', on_delete=PROTECT)
    count = models.IntegerField(null=True)

    def __str__(self):
        min_v = self.min_variant_id
        max_v = self.max_variant_id
        return f"AnnotationRangeLock: (v. {self.version}) {min_v} - {max_v}"


class AnnotationRun(TimeStampedModel):
    status = models.CharField(max_length=1, choices=AnnotationStatus.choices, default=AnnotationStatus.CREATED)
    annotation_range_lock = models.OneToOneField(AnnotationRangeLock, null=True, on_delete=CASCADE)
    # task_id is used as a lock to prevent multiple Celery jobs from executing same job
    task_id = models.CharField(max_length=36, null=True)
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
    annotated_count = models.IntegerField(null=True)
    celery_task_logs = models.JSONField(null=False, default=dict)  # Key=task_id, so we keep logs from multiple runs

    @property
    def variant_annotation_version(self):
        return self.annotation_range_lock.version

    def save(self, **kwargs):
        self.status = self.get_status()
        super().save(**kwargs)

    def get_status(self):
        status = AnnotationStatus.CREATED
        if self.error_exception:
            status = AnnotationStatus.ERROR
        else:
            if self.dump_start:
                status = AnnotationStatus.DUMP_STARTED
            if self.dump_end:
                status = AnnotationStatus.DUMP_COMPLETED
            if self.annotation_start:
                status = AnnotationStatus.ANNOTATION_STARTED
            if self.annotation_end:
                status = AnnotationStatus.ANNOTATION_COMPLETED
            if self.upload_start:
                status = AnnotationStatus.UPLOAD_STARTED
            if self.upload_end:
                status = AnnotationStatus.FINISHED
        return status

    @property
    def annotation_consortium(self):
        return self.annotation_range_lock.version.annotation_consortium

    @property
    def genome_build(self):
        return self.annotation_range_lock.version.genome_build

    def delete_related_objects(self):
        from annotation.annotation_version_querysets import get_queryset_for_annotation_version

        annotation_version = self.variant_annotation_version.get_any_annotation_version()
        for klass in [VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap]:
            qs = get_queryset_for_annotation_version(klass, annotation_version)
            qs.filter(annotation_run=self).delete()

    def delete(self, using=None, keep_parents=False):
        self.delete_related_objects()
        super().delete(using=using, keep_parents=keep_parents)

    def set_task_log(self, key, value):
        assert self.task_id is not None
        task_log = self.celery_task_logs.get(self.task_id, {})
        task_log[key] = value

    def __str__(self):
        return f"AnnotationRun: {localtime(self.modified)} ({self.status})"


class AbstractVariantAnnotation(models.Model):
    """ Common fields between VariantAnnotation and VariantTranscriptAnnotation
        These fields are PER-TRANSCRIPT """
    version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    annotation_run = models.ForeignKey(AnnotationRun, on_delete=CASCADE)

    # gene/transcript are set by VEP and don't have a version.
    # can be set when up/downstream of gene (in which case HGVS_C and thus transcript_version will be null)
    gene = models.ForeignKey(Gene, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, on_delete=SET_NULL)
    # Linked from HGVS transcript (to get version) @see BulkVEPVCFAnnotationInserter.get_transcript_version_id
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, on_delete=SET_NULL)

    # VEP Fields
    # The best way to see how these map to VEP fields is via the annotation details page
    amino_acids = models.TextField(null=True, blank=True)
    cadd_phred = models.FloatField(null=True, blank=True)
    canonical = models.BooleanField(null=True, blank=True)
    nmd_escaping_variant = models.BooleanField(null=True, blank=True)
    codons = models.TextField(null=True, blank=True)
    consequence = models.TextField(null=True, blank=True)
    distance = models.IntegerField(null=True, blank=True)
    domains = models.TextField(null=True, blank=True)
    ensembl_protein = models.TextField(null=True, blank=True)
    exon = models.TextField(null=True, blank=True)
    fathmm_pred_most_damaging = models.CharField(max_length=1, choices=FATHMMPrediction.CHOICES, null=True, blank=True)
    flags = models.TextField(null=True, blank=True)
    gerp_pp_rs = models.FloatField(null=True, blank=True)
    grantham = models.IntegerField(null=True, blank=True)
    hgvs_c = models.TextField(null=True, blank=True)
    hgvs_p = models.TextField(null=True, blank=True)
    impact = models.CharField(max_length=1, choices=PathogenicityImpact.CHOICES, null=True, blank=True)
    interpro_domain = models.TextField(null=True, blank=True)
    intron = models.TextField(null=True, blank=True)
    # LOFTEE
    lof = models.CharField(max_length=1, choices=LOFTEEConfidence.choices, null=True, blank=True)
    lof_filter = models.TextField(null=True, blank=True)
    lof_flags = models.TextField(null=True, blank=True)
    lof_info = models.TextField(null=True, blank=True)
    # LoF Tool
    loftool = models.FloatField(null=True, blank=True)
    maxentscan_alt = models.FloatField(null=True, blank=True)
    maxentscan_diff = models.FloatField(null=True, blank=True)
    maxentscan_ref = models.FloatField(null=True, blank=True)
    maxentscan_percent_diff_ref = models.FloatField(null=True, blank=True)
    mutation_assessor_pred_most_damaging = models.CharField(max_length=1, choices=MutationAssessorPrediction.CHOICES, null=True, blank=True)
    mutation_taster_pred_most_damaging = models.CharField(max_length=1, choices=MutationTasterPrediction.CHOICES, null=True, blank=True)
    polyphen2_hvar_pred_most_damaging = models.CharField(max_length=1, choices=Polyphen2Prediction.CHOICES, null=True, blank=True)
    # protein_position = text as it can be eg indel: "22-23" or splicing: "?-10" or "10-?"
    protein_position = models.TextField(null=True, blank=True)
    revel_score = models.FloatField(null=True, blank=True)
    sift = models.CharField(max_length=1, choices=SIFTPrediction.CHOICES, null=True, blank=True)
    splice_region = models.TextField(null=True, blank=True)
    symbol = models.TextField(null=True, blank=True)

    class Meta:
        abstract = True

    @property
    def transcript_accession(self):
        """ Get transcript_id (with version if possible) """
        if self.transcript_version:
            return self.transcript_version.accession
        if self.transcript:
            return self.transcript.identifier
        return None

    @staticmethod
    def get_domains_components(domains):
        domain_dict = defaultdict(list)
        for d in domains.split("&"):
            key, *values = d.split(":")
            for v in values:
                domain_dict[key].append(v)
        return domain_dict

    @staticmethod
    def amino_acid_3_to_1(protein_string: str) -> str:
        """ p.Ala2351Thr -> p.A2351T

            We can't use pyhgvs.HGVSName(hgvs_string) as pyhgvs gives:
            InvalidHGVSName: Invalid HGVS protein allele "His1753_Lys1780del"
        """
        aa_3_to_1 = invert_dict(protein_letters_1to3)
        aa_3_regex_or = "|".join(aa_3_to_1)

        protein_aa1 = protein_string
        for aa3 in re.findall(f"({aa_3_regex_or})", protein_string):
            protein_aa1 = protein_aa1.replace(aa3, aa_3_to_1[aa3])
        return protein_aa1

    @staticmethod
    def protein_position_to_int(protein_position: str) -> int:
        """ VEP protein position can be:
            185, 185-187, ?-187, 185-? Takes 1st integer it finds """

        for i in protein_position.split("-"):
            try:
                return int(i)
            except ValueError:
                pass  # Skip "?"
        raise ValueError(f"Unable to handle protein_position={protein_position}")


class VariantAnnotation(AbstractVariantAnnotation):
    """ This is the "representative transcript" chosen (1 per variant/annotation version) """
    GENE_COLUMN = "variantannotation__gene"

    # Population frequency
    af_1kg = models.FloatField(null=True, blank=True)
    af_uk10k = models.FloatField(null=True, blank=True)
    gnomad_af = models.FloatField(null=True, blank=True)
    gnomad2_liftover_af = models.FloatField(null=True, blank=True)
    gnomad_ac = models.IntegerField(null=True, blank=True)
    gnomad_an = models.IntegerField(null=True, blank=True)
    gnomad_hom_alt = models.IntegerField(null=True, blank=True)
    gnomad_afr_af = models.FloatField(null=True, blank=True)
    gnomad_amr_af = models.FloatField(null=True, blank=True)
    gnomad_asj_af = models.FloatField(null=True, blank=True)
    gnomad_eas_af = models.FloatField(null=True, blank=True)
    gnomad_fin_af = models.FloatField(null=True, blank=True)
    gnomad_nfe_af = models.FloatField(null=True, blank=True)
    gnomad_oth_af = models.FloatField(null=True, blank=True)
    gnomad_sas_af = models.FloatField(null=True, blank=True)
    gnomad_popmax_af = models.FloatField(null=True, blank=True)
    gnomad_popmax_ac = models.IntegerField(null=True, blank=True)
    gnomad_popmax_an = models.IntegerField(null=True, blank=True)
    gnomad_popmax_hom_alt = models.IntegerField(null=True, blank=True)
    topmed_af = models.FloatField(null=True, blank=True)
    gnomad_filtered = models.BooleanField(null=True, blank=True)
    gnomad_popmax = models.CharField(max_length=3, choices=GnomADPopulation.choices, null=True, blank=True)

    # From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267638/
    # "optimum cutoff value identified in the ROC analysis (0.6)"
    dbscsnv_ada_score = models.FloatField(null=True, blank=True)
    dbscsnv_rf_score = models.FloatField(null=True, blank=True)
    dbsnp_rs_id = models.TextField(null=True, blank=True)
    cosmic_id = models.TextField(null=True, blank=True)  # COSV - Genomic Mutation ID
    cosmic_legacy_id = models.TextField(null=True, blank=True)  # COSM
    cosmic_count = models.IntegerField(null=True, blank=True)
    pubmed = models.TextField(null=True, blank=True)
    # Mastermind Cited Variants Reference. @see https://www.genomenon.com/cvr/
    mastermind_count_1_cdna = models.IntegerField(null=True, blank=True)
    mastermind_count_2_cdna_prot = models.IntegerField(null=True, blank=True)
    mastermind_count_3_aa_change = models.IntegerField(null=True, blank=True)
    mastermind_mmid3 = models.TextField(null=True, blank=True)  # gene:key for mastermind_count_3_aa_change

    # dbNSFP pathogenicity rank scores (variant only)
    cadd_raw_rankscore = models.FloatField(null=True, blank=True)
    revel_rankscore = models.FloatField(null=True, blank=True)
    bayesdel_noaf_rankscore = models.FloatField(null=True, blank=True)
    clinpred_rankscore = models.FloatField(null=True, blank=True)
    vest4_rankscore = models.FloatField(null=True, blank=True)
    metalr_rankscore = models.FloatField(null=True, blank=True)

    # Not all builds have all phylop/phastcons
    phylop_30_way_mammalian = models.FloatField(null=True, blank=True)
    phylop_46_way_mammalian = models.FloatField(null=True, blank=True)
    phylop_100_way_vertebrate = models.FloatField(null=True, blank=True)
    phastcons_30_way_mammalian = models.FloatField(null=True, blank=True)
    phastcons_46_way_mammalian = models.FloatField(null=True, blank=True)
    phastcons_100_way_vertebrate = models.FloatField(null=True, blank=True)

    # SpliceAI - @see https://pubmed.ncbi.nlm.nih.gov/30661751/
    # This has so many columns, perhaps a highest score and summary col of "AG: -5, AL: 35, DG: -5, DL: 17"
    # 0.2 = high recal, 0.5 = recommended, 0.8 = high precision
    spliceai_pred_dp_ag = models.IntegerField(null=True, blank=True)
    spliceai_pred_dp_al = models.IntegerField(null=True, blank=True)
    spliceai_pred_dp_dg = models.IntegerField(null=True, blank=True)
    spliceai_pred_dp_dl = models.IntegerField(null=True, blank=True)
    spliceai_pred_ds_ag = models.FloatField(null=True, blank=True)
    spliceai_pred_ds_al = models.FloatField(null=True, blank=True)
    spliceai_pred_ds_dg = models.FloatField(null=True, blank=True)
    spliceai_pred_ds_dl = models.FloatField(null=True, blank=True)
    spliceai_gene_symbol = models.TextField(null=True, blank=True)
    repeat_masker = models.TextField(null=True, blank=True)
    overlapping_symbols = models.TextField(null=True, blank=True)
    # Summary of most_damaging fields for faster DamageNode queries
    predictions_num_pathogenic = models.IntegerField(default=0)
    predictions_num_benign = models.IntegerField(default=0)

    somatic = models.BooleanField(null=True, blank=True)
    variant_class = models.CharField(max_length=2, choices=VariantClass.choices, null=True, blank=True)
    vep_skipped_reason = models.CharField(max_length=1, choices=VEPSkippedReason.choices, null=True, blank=True)

    GNOMAD_FIELDS = {
        GnomADPopulation.AFRICAN_AFRICAN_AMERICAN: 'gnomad_afr_af',
        GnomADPopulation.LATINO: 'gnomad_amr_af',
        GnomADPopulation.ASHKENAZI_JEWISH: 'gnomad_asj_af',
        GnomADPopulation.EAST_ASIAN: 'gnomad_eas_af',
        GnomADPopulation.FINNISH: 'gnomad_fin_af',
        GnomADPopulation.NON_FINNISH_EUROPEAN: 'gnomad_nfe_af',
        GnomADPopulation.OTHER: 'gnomad_oth_af',
        GnomADPopulation.SOUTH_ASIAN: 'gnomad_sas_af',
    }

    MASTERMIND_FIELDS = {
        "mastermind_count_1_cdna": "cDNA",
        "mastermind_count_2_cdna_prot": "cDNA/Prot",
        "mastermind_count_3_aa_change": "AA change",
        "mastermind_mmid3": "MMID3",
    }

    SPLICEAI_DS_DP = {
        "AG": ("spliceai_pred_ds_ag", "spliceai_pred_dp_ag"),
        "AL": ("spliceai_pred_ds_al", "spliceai_pred_dp_al"),
        "DG": ("spliceai_pred_ds_dg", "spliceai_pred_dp_dg"),
        "DL": ("spliceai_pred_ds_dl", "spliceai_pred_dp_dl"),
    }

    # List of filters to describe variants that can be annotated
    VARIANT_ANNOTATION_Q = [
        Variant.get_no_reference_q(),
        ~Q(alt__seq__in=['.', '*', "<DEL>"]),  # Exclude non-standard variants
    ]

    @lazy
    def has_extended_gnomad_fields(self):
        """ I grabbed a few new fields but haven't patched back to GRCh37 yet
            TODO: remove this and if statements in variant_details.html once issue #231 is completed """
        extended_fields = ["gnomad2_liftover_af", "gnomad_ac", "gnomad_an", "gnomad_popmax_ac",
                           "gnomad_popmax_an", "gnomad_popmax_hom_alt"]
        return any(getattr(self, f) is not None for f in extended_fields)

    @property
    def gnomad_url(self):
        GNOMAD2 = "gnomad_r2_1"
        GNOMAD3 = "gnomad_r3"

        gnomad_dataset = None
        if self.version.gnomad.startswith("3"):
            if self.gnomad_af is not None:
                gnomad_dataset = GNOMAD3
            elif self.gnomad2_liftover_af is not None:
                gnomad_dataset = GNOMAD2
        elif self.gnomad_af and self.version.gnomad.startswith("2.1"):
            gnomad_dataset = GNOMAD2

        url = None
        if gnomad_dataset:
            v = self.variant
            gnomad_variant = f"{v.locus.chrom}-{v.locus.position}-{v.locus.ref}-{v.alt}"
            url = f"http://gnomad.broadinstitute.org/variant/{gnomad_variant}?dataset={gnomad_dataset}"
        return url

    @property
    def mastermind_url(self):
        return self.get_mastermind_url(self.mastermind_mmid3)

    @staticmethod
    def get_mastermind_url(mastermind_mmid3):
        if mastermind_mmid3:
            return f"https://mastermind.genomenon.com/detail?mutation={mastermind_mmid3}"
        return None

    def has_spliceai(self):
        return any((self.spliceai_pred_ds_ag, self.spliceai_pred_ds_al,
                    self.spliceai_pred_ds_dg, self.spliceai_pred_ds_dl))

    @staticmethod
    def get_gnomad_population_field(population):
        return VariantAnnotation.GNOMAD_FIELDS.get(population)

    @lazy
    def transcript_annotation(self) -> List['VariantTranscriptAnnotation']:
        return self.variant.varianttranscriptannotation_set.filter(version=self.version)

    def get_search_terms(self):
        extra_terms = []
        if self.dbsnp_rs_id:
            extra_terms.append(self.dbsnp_rs_id)
        return get_variant_search_terms(self.transcript_annotation, extra_terms=extra_terms)

    def get_pubmed_search_terms(self):
        return get_variant_pubmed_search_terms(self.transcript_annotation)

    def __str__(self):
        return f"{self.variant}: {self.version}"

    class Meta:
        unique_together = ("version", "variant")


class VariantTranscriptAnnotation(AbstractVariantAnnotation):
    """ This is ALL transcript predictions for a variant (>=1 per variant) """

    class Meta:
        unique_together = ("version", "variant", "transcript_version")

    @staticmethod
    def get_overlapping_genes_q(variant_annotation_version: VariantAnnotationVersion, gene_ids_qs) -> Q:
        """ Returns Q object where any transcript for a variant matches gene in genes
            Variants can overlap with multiple genes, and the VariantAnnotation (ie "pick" or representative annotation)
            may not be the one in the gene list. Thus we have to check everything that was in transcript annotation too
        """

        vto = VariantGeneOverlap.objects.filter(version=variant_annotation_version, gene__in=gene_ids_qs)
        # Use pk__in so we don't return multiple records per variant
        return Q(pk__in=vto.values_list("variant_id", flat=True))


class VariantGeneOverlap(models.Model):
    """ Created for every gene that overlaps a variant (per annotation version)
        Allows efficient gene list queries (10% larger than variant annotation) while handling multiple transcripts """
    version = models.ForeignKey(VariantAnnotationVersion, on_delete=CASCADE)
    annotation_run = models.ForeignKey(AnnotationRun, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)

    class Meta:
        unique_together = ('version', 'variant', 'annotation_run', 'gene')


class ManualVariantEntryCollection(models.Model):
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=CASCADE)
    created = models.DateTimeField(auto_now_add=True)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    celery_task = models.CharField(max_length=36, null=True)

    @lazy
    def first_entry(self) -> Optional['ManualVariantEntry']:
        """ Often a user will enter a single variant into search and then wait for it to be created/annotated
            We have the "first" one to deal with this """
        return self.manualvariantentry_set.all().order_by("pk").first()

    @lazy
    def first_variant(self) -> Optional['Variant']:
        variant = None
        if mve := self.first_entry:
            if cve := mve.unique_created_variants.filter(variant__isnull=False).order_by("variant_id").first():
                variant = cve.variant
        return variant

    def first_variant_annotation_run(self) -> Optional['AnnotationRun']:
        annotation_run = None
        if mve := self.first_entry:
            if variant := self.first_variant:
                ar_qs = AnnotationRun.objects.filter(annotation_range_lock__version__genome_build=mve.genome_build,
                                                     annotation_range_lock__min_variant__gte=variant.pk,
                                                     annotation_range_lock__max_variant__lte=variant.pk)
                annotation_run = ar_qs.first()
        return annotation_run

    @staticmethod
    def get_for_user(user, mvec_id):
        mvec = ManualVariantEntryCollection.objects.get(pk=mvec_id)
        if not (user.is_staff or user == mvec.user):
            msg = f"User {user} doesn't have permission to access mvec_id={mvec_id}"
            raise PermissionDenied(msg)
        return mvec

    def __str__(self):
        status = self.get_import_status_display()
        if self.import_status == ImportStatus.IMPORTING:
            status += f", task: {self.celery_task}"

        count = self.manualvariantentry_set.count()
        return f"{self.pk}: {count} manual entries ({status})"


class ManualVariantEntry(models.Model):
    manual_variant_entry_collection = models.ForeignKey(ManualVariantEntryCollection, on_delete=CASCADE)
    line_number = models.IntegerField()
    entry_text = models.TextField()
    entry_type = models.CharField(max_length=1, choices=ManualVariantEntryType.choices, default=ManualVariantEntryType.UNKNOWN)
    error_message = models.TextField(blank=True, null=True)  # Set if any error...

    @property
    def unique_created_variants(self):
        return self.createdmanualvariant_set.distinct('variant_id').all()

    @property
    def genome_build(self) -> GenomeBuild:
        return self.manual_variant_entry_collection.genome_build

    @staticmethod
    def get_entry_type(line: str) -> ManualVariantEntryType:
        HGVS_MINIMAL_PATTERN = re.compile(r"[^:].+:[cgp]\..*\d+.*")
        MATCHERS = {
            DBSNP_PATTERN: ManualVariantEntryType.DBSNP,
            HGVS_MINIMAL_PATTERN: ManualVariantEntryType.HGVS,
            VARIANT_PATTERN: ManualVariantEntryType.VARIANT,
        }
        for k, v in MATCHERS.items():
            if k.match(line):
                return v
        return ManualVariantEntryType.UNKNOWN


class CreatedManualVariant(models.Model):
    manual_variant_entry = models.ForeignKey(ManualVariantEntry, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)


class InvalidAnnotationVersionError(Exception):
    pass


class AnnotationVersionManager(models.Manager):

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.select_related('genome_build', 'variant_annotation_version', 'gene_annotation_version',
                                 'clinvar_version', 'human_protein_atlas_version')


# Use this as a way to keep all the versions together.
# We re-use this, if nobody has referenced it - as you may eg update
# variant/gene/clinvar annotations every 6 months etc and no point having so many sub-versions
class AnnotationVersion(models.Model):
    SUB_ANNOTATIONS = ['variant_annotation_version', 'gene_annotation_version', 'clinvar_version', 'human_protein_atlas_version']
    objects = AnnotationVersionManager()  # Always select_related
    annotation_date = models.DateTimeField(auto_now=True)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, null=True, on_delete=PROTECT)
    gene_annotation_version = models.ForeignKey(GeneAnnotationVersion, null=True, on_delete=PROTECT)
    clinvar_version = models.ForeignKey(ClinVarVersion, null=True, on_delete=PROTECT)
    human_protein_atlas_version = models.ForeignKey(HumanProteinAtlasAnnotationVersion, null=True, on_delete=PROTECT)

    def get_sub_annotation_versions(self):
        sub_annotation_versions = []
        for field in self.SUB_ANNOTATIONS:
            annotation_version = getattr(self, field)
            if annotation_version:
                sub_annotation_versions.append(annotation_version)
        return sub_annotation_versions

    def get_partition_names(self):
        partitions_names = {}
        for sub_version in self.get_sub_annotation_versions():
            for base_table_name in sub_version.RECORDS_BASE_TABLE_NAMES:
                partitions_names[base_table_name] = sub_version.get_partition_table(base_table_name=base_table_name)
        return partitions_names

    def sql_partition_transformer(self, sql):
        """ Modifies SQL generated by QuerySet
             @see library.django_utils.django_queryset_sql_transformer.get_queryset_with_transformer_hook  """

        for annotation_version in self.get_sub_annotation_versions():
            sql = annotation_version.sql_partition_transformer(sql)
        return sql

    @lazy
    def is_valid(self) -> bool:
        try:
            self.validate()
            return True
        except InvalidAnnotationVersionError:
            return False

    def validate(self):
        missing_sub_annotations = []
        for field in self.SUB_ANNOTATIONS:
            field_fk = f"{field}_id"  # Avoid fetching related data
            sub_annotation = getattr(self, field_fk)
            if sub_annotation is None:
                missing_sub_annotations.append(field)
        if missing_sub_annotations:
            missing = ", ".join([str(s) for s in missing_sub_annotations])
            raise InvalidAnnotationVersionError(f"AnnotationVersion: {self} missing sub annotations: {missing}")

        if vav_gar_id := self.variant_annotation_version.gene_annotation_release_id:
            gene_gar_id = self.gene_annotation_version.gene_annotation_release_id
            if vav_gar_id != gene_gar_id:
                different_msg = f"GeneAnnotationRelease is different fro Variant {vav_gar_id} vs Gene: {gene_gar_id}"
                raise InvalidAnnotationVersionError(different_msg)

    @staticmethod
    def latest(genome_build: GenomeBuild, validate=True, active=True) -> Optional['AnnotationVersion']:
        av_qs = AnnotationVersion.objects.filter(genome_build=genome_build)
        if active:
            av_qs = av_qs.exclude(variant_annotation_version__active=False)
        av: AnnotationVersion = av_qs.order_by("annotation_date").last()
        if validate:
            if av is None:
                raise AnnotationVersion.DoesNotExist(f"Warning: GenomeBuild {genome_build} has no annotation version!")
            av.validate()
        return av

    @staticmethod
    @transaction.atomic
    def new_sub_version(sub_version_genome_build):
        logging.info("new_sub_version genome_build=%s", sub_version_genome_build)

        def latest(klass):
            return klass.objects.order_by('annotation_date').last()

        def latest_for_build(klass, genome_build: GenomeBuild, genome_build_path: str = "genome_build"):
            qs = klass.objects.filter(**{genome_build_path: genome_build})
            return qs.order_by('annotation_date').last()

        if sub_version_genome_build:
            builds = [sub_version_genome_build]  # Just this one
        else:
            builds = GenomeBuild.builds_with_annotation()

        for genome_build in builds:
            kwargs = {
                "genome_build": genome_build,
                "variant_annotation_version": latest_for_build(VariantAnnotationVersion, genome_build),
                "gene_annotation_version": latest_for_build(GeneAnnotationVersion, genome_build, "gene_annotation_release__genome_build"),
                "clinvar_version": latest_for_build(ClinVarVersion, genome_build),
                "human_protein_atlas_version": latest(HumanProteinAtlasAnnotationVersion)
            }

            # There is a slight race condition here between is_referenced and saving the object, but
            # You've shut the server down to upgrade the annotations, right? Right??
            last_annotation_version = AnnotationVersion.latest(genome_build, validate=False)
            if last_annotation_version is None or object_is_referenced(last_annotation_version):
                av, created = AnnotationVersion.objects.get_or_create(**kwargs)
                if not created:
                    logging.warning(
                        "new_sub_version() found existing AnnotationVersion, probably shouldn't have been called?")
                else:
                    logging.info("Created new AnnotationVersion %s", av.pk)
            else:
                logging.info("Nothing referencing last annotation version, so overwrite")
                av = AnnotationVersion(**kwargs)
                av.pk = last_annotation_version.pk
                av.save()

                logging.info("Annotation (%s) hasn't been referenced by anything, upgrading it.", av)

    def get_variant_annotation(self):
        return VariantAnnotation.objects.filter(version=self.variant_annotation_version)

    def get_gene_annotation(self):
        return GeneAnnotation.objects.filter(version=self.gene_annotation_version)

    def get_clinvar(self):
        return ClinVar.objects.filter(version=self.clinvar_version)

    def get_human_protein_atlas_annotation(self):
        return HumanProteinAtlasAnnotation.objects.filter(version=self.human_protein_atlas_version)

    def __str__(self):
        sub_versions = [f"Variant: {self.variant_annotation_version}",
                        f"Gene: {self.gene_annotation_version}",
                        f"ClinVar: {self.clinvar_version}",
                        f"HPA: {self.human_protein_atlas_version}"]
        sub_versions_str = ", ".join(sub_versions)
        return f"{self.pk} ({self.annotation_date.date()}). {sub_versions_str}"


class CachedWebResource(TimeStampedModel):
    """ These are annotations that can populate themselves via the web

        * Can be manually updated via button kicking off annotation.views.load_cached_web_resource
        * List of names each system uses is in settings
        * Creation of one with name triggers signals which apps can configure to listen to / run populators
    """
    name = models.TextField(primary_key=True)
    description = models.TextField(blank=True)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)

    @staticmethod
    def named_handler_factory(name, celery_task_class):
        def handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
            created = kwargs.get("created")

            if created and instance.name == name:
                task_signature = celery_task_class.si(instance.pk)
                task_signature.apply_async()

        return handler


class GeneSymbolCitation(models.Model):
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    citation = models.ForeignKey(Citation, on_delete=CASCADE)

    class Meta:
        unique_together = ('gene_symbol', 'citation')


class GeneSymbolPubMedCount(TimeStampedModel):
    gene_symbol = models.OneToOneField(GeneSymbol, on_delete=CASCADE)
    count = models.IntegerField()

    @staticmethod
    def get_for_gene_symbol(gene_symbol_id):
        max_days = settings.ANNOTATION_PUBMED_GENE_SYMBOL_COUNT_CACHE_DAYS
        latest_time = timezone.now() - timedelta(days=max_days)
        try:
            return GeneSymbolPubMedCount.objects.get(gene_symbol_id=gene_symbol_id,
                                                     modified__gte=latest_time)
        except GeneSymbolPubMedCount.DoesNotExist:
            pass

        h = Entrez.esearch(db="pubmed", term=gene_symbol_id, rettype='count')
        record = Entrez.read(h)
        count = record["Count"]
        return GeneSymbolPubMedCount.objects.get_or_create(gene_symbol_id=gene_symbol_id,
                                                           defaults={"count": count})[0]

    def __str__(self):
        return f"{self.gene_symbol_id}: {self.count}"


class MutationalSignatureInfo(models.Model):
    signature_id = models.IntegerField(primary_key=True)
    cancer_types = models.TextField()
    proposed_aetiology = models.TextField()
    additional_mutational_features = models.TextField()
    comments = models.TextField()
