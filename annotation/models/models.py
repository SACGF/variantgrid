import json
import logging
import os
import re
from collections import defaultdict
from datetime import datetime, timedelta
from typing import List

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
from django.utils.timesince import timesince
from django.utils.timezone import localtime
from django_extensions.db.models import TimeStampedModel
from lazy import lazy

from annotation.external_search_terms import get_variant_search_terms, get_variant_pubmed_search_terms
from annotation.models.damage_enums import Polyphen2Prediction, FATHMMPrediction, MutationTasterPrediction, \
    SIFTPrediction, PathogenicityImpact, MutationAssessorPrediction
from annotation.models.models_enums import HumanProteinAtlasAbundance, AnnotationStatus, CitationSource, \
    TranscriptStatus, GenomicStrand, ClinGenClassification, VariantClass, ColumnAnnotationCategory, VEPPlugin, \
    VEPCustom, ClinVarReviewStatus, VEPSkippedReason
from annotation.models.models_mim_hpo import MIMMorbid, HPOSynonym, MIMMorbidAlias
from genes.models import GeneSymbol, Gene, TranscriptVersion, Transcript, GeneAnnotationRelease
from genes.models_enums import AnnotationConsortium
from library.django_utils import object_is_referenced
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.utils import invert_dict
from patients.models_enums import GnomADPopulation
from snpdb.models import GenomeBuild, Variant, VariantGridColumn, Q, VCF
from snpdb.models.models_enums import ImportStatus


class MIMGene(models.Model):
    mim_morbid = models.ForeignKey(MIMMorbid, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)

    class Meta:
        unique_together = ("mim_morbid", "gene")

    def __str__(self):
        return f"{self.mim_morbid}: {self.gene}"


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
def clinvar_version_pre_delete_handler(sender, instance, **kwargs):
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
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    clinvar_variation_id = models.IntegerField()
    clinvar_allele_id = models.IntegerField()
    clinvar_preferred_disease_name = models.TextField(null=True, blank=True)
    clinvar_disease_database_name = models.TextField(null=True, blank=True)
    clinvar_review_status = models.CharField(max_length=1, choices=ClinVarReviewStatus.CHOICES)
    clinical_significance = models.TextField(null=True, blank=True)  # Multiple values of above
    highest_pathogenicity = models.IntegerField(default=0)  # Highest of clinical_significance
    clinvar_clinical_sources = models.TextField(null=True, blank=True)
    clinvar_origin = models.IntegerField(default=0)
    clinvar_suspect_reason_code = models.IntegerField(default=0)
    drug_response = models.BooleanField(default=False)

    @property
    def stars(self):
        return ClinVarReviewStatus.STARS[self.clinvar_review_status]

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
    pass


class Citation(models.Model):
    citation_source = models.CharField(max_length=1, choices=CitationSource.CHOICES)
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

    @staticmethod
    def citations_from_text(text):
        """ returns a list of (unsaved) Citation objects from text """
        citation_source_codes = dict({k.lower(): v for (k, v) in CitationSource.CODES.items()})
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


class EnsemblGeneAnnotationVersion(SubVersionPartition):
    RECORDS_BASE_TABLE_NAMES = ["annotation_ensemblgeneannotation"]
    filename = models.TextField()
    md5_hash = models.CharField(max_length=32)
    ensembl_version = models.IntegerField()
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)


@receiver(pre_delete, sender=EnsemblGeneAnnotationVersion)
def ensembl_gene_annotation_version_pre_delete_handler(sender, instance, **kwargs):
    instance.delete_related_objects()


class EnsemblGeneAnnotation(models.Model):
    """ Populated from SACGF Gene Level Annotation (Frank Feng) """
    version = models.ForeignKey(EnsemblGeneAnnotationVersion, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    hgnc_symbol = models.TextField(null=True)
    external_gene_name = models.TextField(null=True)
    hgnc_symbol_lower = models.TextField(null=True)
    hgnc_name = models.TextField(null=True)
    synonyms = models.TextField(null=True)
    previous_symbols = models.TextField(null=True)
    hgnc_chromosome = models.TextField(null=True)
    gene_family_tag = models.TextField(null=True)
    gene_family_description = models.TextField(null=True)
    hgnc_id = models.TextField(null=True)
    entrez_gene_id = models.IntegerField(null=True)
    uniprot_id = models.TextField(null=True)
    ucsc_id = models.TextField(null=True)
    omim_id = models.TextField(null=True)  # Can have multiple (comma sep)
    enzyme_ids = models.TextField(null=True)
    ccds_ids = models.TextField(null=True)
    rgd_id = models.TextField(null=True)
    mgi_id = models.TextField(null=True)
    rvis_percentile = models.TextField(null=True)
    refseq_gene_summary = models.TextField(null=True)
    function_from_uniprotkb = models.TextField(null=True)
    pathway_from_uniprotkb = models.TextField(null=True)
    tissue_specificity_from_uniprotkb = models.TextField(null=True)
    phenotypes_from_ensembl = models.TextField(null=True)
    omim_phenotypes = models.TextField(null=True)
    gene_biotype = models.TextField(null=True)
    status = models.CharField(max_length=1, choices=TranscriptStatus.CHOICES, null=True)
    chromosome_name = models.TextField(null=True)
    start_position = models.IntegerField(null=True)
    end_position = models.IntegerField(null=True)
    band = models.TextField(null=True)
    strand = models.CharField(max_length=1, choices=GenomicStrand.CHOICES, null=True)
    percentage_gc_content = models.FloatField(null=True)
    transcript_count = models.IntegerField(null=True)
    in_cancer_gene_census = models.BooleanField(null=True)

    class Meta:
        unique_together = ("version", "gene")

    @staticmethod
    def latest_qs():
        latest_version = EnsemblGeneAnnotationVersion.objects.order_by("pk").last()
        return EnsemblGeneAnnotation.objects.filter(version=latest_version)

    @staticmethod
    def get_for_symbol(genome_build: GenomeBuild, gene_symbol: GeneSymbol):
        annotation_version = AnnotationVersion.latest(genome_build)
        ega_qs = annotation_version.get_ensembl_gene_annotation()
        ENSEMBL_ANNOTATION = [
            Q(hgnc_symbol=gene_symbol),
            Q(gene__in=gene_symbol.get_genes().filter(annotation_consortium=AnnotationConsortium.ENSEMBL)),
        ]
        gene_annotation = None
        for q in ENSEMBL_ANNOTATION:
            if gene_annotation := ega_qs.filter(q).first():
                break
        return gene_annotation


class HumanProteinAtlasAnnotationVersion(SubVersionPartition):
    RECORDS_BASE_TABLE_NAMES = ["annotation_humanproteinatlasannotation"]
    filename = models.TextField()
    md5_hash = models.CharField(max_length=32)
    hpa_version = models.FloatField()


@receiver(pre_delete, sender=HumanProteinAtlasAnnotationVersion)
def human_protein_annotation_version_pre_delete_handler(sender, instance, **kwargs):
    instance.delete_related_objects()


class HumanProteinAtlasTissueSample(models.Model):
    name = models.TextField()

    def __str__(self):
        return self.name


class HumanProteinAtlasAnnotation(models.Model):
    version = models.ForeignKey(HumanProteinAtlasAnnotationVersion, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)  # Always Ensembl
    tissue_sample = models.ForeignKey(HumanProteinAtlasTissueSample, on_delete=CASCADE)
    abundance = models.CharField(max_length=1, choices=HumanProteinAtlasAbundance.CHOICES)


class ColumnVEPField(models.Model):
    """ For VariantAnnotation/Transcript columns derived from VEP fields """
    column = models.TextField(unique=True)
    variant_grid_column = models.OneToOneField(VariantGridColumn, null=True, on_delete=SET_NULL)
    category = models.CharField(max_length=1, choices=ColumnAnnotationCategory.CHOICES)
    source_field = models.TextField(null=True)  # @see use vep_info_field
    source_field_processing_description = models.TextField(null=True)
    vep_plugin = models.CharField(max_length=1, choices=VEPPlugin.CHOICES, null=True)
    vep_custom = models.CharField(max_length=1, choices=VEPCustom.CHOICES, null=True)
    source_field_has_custom_prefix = models.BooleanField(default=False)

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
    def get_source_fields(**columnvepfield_kwargs):
        qs = ColumnVEPField.objects.filter(**columnvepfield_kwargs).distinct("source_field")
        return list(qs.values_list("source_field", flat=True).order_by("source_field"))


class VariantAnnotationVersion(SubVersionPartition):
    REPRESENTATIVE_TRANSCRIPT_ANNOTATION = "annotation_variantannotation"
    TRANSCRIPT_ANNOTATION = "annotation_varianttranscriptannotation"
    VARIANT_GENE_OVERLAP = "annotation_variantgeneoverlap"
    RECORDS_BASE_TABLE_NAMES = [REPRESENTATIVE_TRANSCRIPT_ANNOTATION, TRANSCRIPT_ANNOTATION, VARIANT_GENE_OVERLAP]

    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    annotation_consortium = models.CharField(max_length=1, choices=AnnotationConsortium.CHOICES)
    # GeneAnnotationRelease - imported GTF we can use to get gene/transcript versions that match VEP
    gene_annotation_release = models.ForeignKey(GeneAnnotationRelease, null=True, on_delete=CASCADE)
    last_checked_date = models.DateTimeField(null=True)

    vep = models.IntegerField()
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

    def get_annotation_status_info_for_variant(self, variant):
        """ returns (annotation_status_display, time_since_last_modified) """
        qs = self.annotationrangelock_set.filter(min_variant__pk__gte=variant.pk,
                                                 max_variant__pk__lte=variant.pk)
        annotation_range_lock = qs.order_by("pk").last()
        annotation_status_display = "No annotation"
        time_since_last_modified = "n/a"
        if annotation_range_lock:
            try:
                annotation_run = annotation_range_lock.annotationrun
                annotation_status_display = annotation_run.get_status_display()
                time_since_last_modified = timesince(annotation_run.modified)
            except AnnotationRun.DoesNotExist:
                pass
        return annotation_status_display, time_since_last_modified

    @staticmethod
    def latest(genome_build):
        return VariantAnnotationVersion.objects.filter(genome_build=genome_build).order_by("annotation_date").last()

    def get_any_annotation_version(self):
        """ Often you don't care what annotation version you use, only that variant annotation version is this one """
        return self.annotationversion_set.last()

    def __str__(self):
        super_str = super().__str__()
        return f"{super_str} VEP: {self.vep}/{self.get_annotation_consortium_display()}/{self.genome_build.name}"


@receiver(pre_delete, sender=VariantAnnotationVersion)
def variant_annotation_version_pre_delete_handler(sender, instance, **kwargs):
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
    status = models.CharField(max_length=1, choices=AnnotationStatus.CHOICES, default=AnnotationStatus.CREATED)
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

    def delete(self, **kwargs):
        self.delete_related_objects()
        super().delete(**kwargs)

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
    transcript = models.ForeignKey(Transcript, null=True, on_delete=CASCADE)
    # Linked from HGVS transcript (to get version) @see BulkVEPVCFAnnotationInserter.get_transcript_version_id
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, on_delete=CASCADE)

    # VEP Fields
    # The best way to see how these map to VEP fields is via the annotation details page
    amino_acids = models.TextField(null=True, blank=True)
    cadd_phred = models.FloatField(null=True, blank=True)
    cadd_raw = models.FloatField(null=True, blank=True)
    canonical = models.BooleanField(null=True, blank=True)
    codons = models.TextField(null=True, blank=True)
    consequence = models.TextField(null=True, blank=True)
    dbsnp_rs_id = models.TextField(null=True, blank=True)
    distance = models.IntegerField(null=True, blank=True)
    domains = models.TextField(null=True, blank=True)
    ensembl_protein = models.TextField(null=True, blank=True)
    exon = models.TextField(null=True, blank=True)
    fathmm_pred = models.TextField(null=True, blank=True)
    fathmm_pred_most_damaging = models.CharField(max_length=1, choices=FATHMMPrediction.CHOICES, null=True, blank=True)
    flags = models.TextField(null=True, blank=True)
    gene_text = models.TextField(null=True, blank=True)
    gerp_pp_rs = models.FloatField(null=True, blank=True)
    grantham = models.IntegerField(null=True, blank=True)
    hgnc_id = models.IntegerField(null=True, blank=True)
    hgvs_c = models.TextField(null=True, blank=True)
    hgvs_p = models.TextField(null=True, blank=True)
    impact = models.CharField(max_length=1, choices=PathogenicityImpact.CHOICES, null=True, blank=True)
    interpro_domain = models.TextField(null=True, blank=True)
    intron = models.TextField(null=True, blank=True)
    loftool = models.FloatField(null=True, blank=True)
    maxentscan_alt = models.FloatField(null=True, blank=True)
    maxentscan_diff = models.FloatField(null=True, blank=True)
    maxentscan_ref = models.FloatField(null=True, blank=True)
    maxentscan_percent_diff_ref = models.FloatField(null=True, blank=True)
    mutation_assessor_pred = models.TextField(null=True, blank=True)
    mutation_assessor_pred_most_damaging = models.CharField(max_length=1, choices=MutationAssessorPrediction.CHOICES, null=True, blank=True)
    mutation_taster_pred = models.TextField(null=True, blank=True)
    mutation_taster_pred_most_damaging = models.CharField(max_length=1, choices=MutationTasterPrediction.CHOICES, null=True, blank=True)
    # Summary of most_damaging fields for faster DamageNode queries
    predictions_num_pathogenic = models.IntegerField(default=0)
    predictions_num_benign = models.IntegerField(default=0)
    # Not all builds have all phylop/phastcons
    phylop_30_way_mammalian = models.FloatField(null=True, blank=True)
    phylop_46_way_mammalian = models.FloatField(null=True, blank=True)
    phylop_100_way_vertebrate = models.FloatField(null=True, blank=True)
    phastcons_30_way_mammalian = models.FloatField(null=True, blank=True)
    phastcons_46_way_mammalian = models.FloatField(null=True, blank=True)
    phastcons_100_way_vertebrate = models.FloatField(null=True, blank=True)
    polyphen2_hvar_pred = models.TextField(null=True, blank=True)
    polyphen2_hvar_pred_most_damaging = models.CharField(max_length=1, choices=Polyphen2Prediction.CHOICES, null=True, blank=True)
    pubmed = models.TextField(null=True, blank=True)
    # protein_position = text as it can be eg indel: "22-23" or splicing: "?-10" or "10-?"
    protein_position = models.TextField(null=True, blank=True)
    revel_score = models.FloatField(null=True, blank=True)
    sift = models.CharField(max_length=1, choices=SIFTPrediction.CHOICES, null=True, blank=True)
    splice_region = models.TextField(null=True, blank=True)
    swissprot = models.TextField(null=True, blank=True)
    symbol = models.TextField(null=True, blank=True)
    symbol_source = models.TextField(null=True, blank=True)
    trembl = models.TextField(null=True, blank=True)
    uniparc = models.TextField(null=True, blank=True)

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

    @lazy
    def flagged_pathogenicity(self):
        PATHOGENICITY_FIELDS = {
            "fathmm_pred_most_damaging": FATHMMPrediction,
            "mutation_assessor_pred_most_damaging": MutationAssessorPrediction,
            "mutation_taster_pred_most_damaging": MutationTasterPrediction,
            "polyphen2_hvar_pred_most_damaging": Polyphen2Prediction,
            "sift": SIFTPrediction,
        }

        fp = {}
        for f, klass in PATHOGENICITY_FIELDS.items():
            level = getattr(self, f)
            fp[f] = klass.is_level_flagged(level)
        return fp

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
    topmed_af = models.FloatField(null=True, blank=True)
    gnomad_filtered = models.BooleanField(null=True, blank=True)
    gnomad_popmax = models.CharField(max_length=3, choices=GnomADPopulation.CHOICES, null=True, blank=True)

    # From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267638/
    # "optimum cutoff value identified in the ROC analysis (0.6)"
    DBSCSNV_CUTOFF = 0.6
    dbscsnv_ada_score = models.FloatField(null=True, blank=True)
    dbscsnv_rf_score = models.FloatField(null=True, blank=True)
    cosmic_id = models.TextField(null=True, blank=True)  # COSV - Genomic Mutation ID
    cosmic_legacy_id = models.TextField(null=True, blank=True)  # COSM
    cosmic_count = models.IntegerField(null=True, blank=True)
    # Mastermind Cited Variants Reference. @see https://www.genomenon.com/cvr/
    mastermind_count_1_cdna = models.IntegerField(null=True, blank=True)
    mastermind_count_2_cdna_prot = models.IntegerField(null=True, blank=True)
    mastermind_count_3_aa_change = models.IntegerField(null=True, blank=True)
    mastermind_mmid3 = models.TextField(null=True, blank=True)  # gene:key for mastermind_count_3_aa_change
    # SpliceAI - @see https://pubmed.ncbi.nlm.nih.gov/30661751/
    # This has so many columns, perhaps a highest score and summary col of "AG: -5, AL: 35, DG: -5, DL: 17"
    SPLICEAI_CUTOFF = 0.5  # 0.2 = high recal, 0.5 = recommended, 0.8 = high precision
    spliceai_pred_dp_ag = models.IntegerField(null=True, blank=True)
    spliceai_pred_dp_al = models.IntegerField(null=True, blank=True)
    spliceai_pred_dp_dg = models.IntegerField(null=True, blank=True)
    spliceai_pred_dp_dl = models.IntegerField(null=True, blank=True)
    spliceai_pred_ds_ag = models.FloatField(null=True, blank=True)
    spliceai_pred_ds_al = models.FloatField(null=True, blank=True)
    spliceai_pred_ds_dg = models.FloatField(null=True, blank=True)
    spliceai_pred_ds_dl = models.FloatField(null=True, blank=True)
    repeat_masker = models.TextField(null=True, blank=True)
    overlapping_symbols = models.TextField(null=True, blank=True)

    somatic = models.BooleanField(null=True, blank=True)
    variant_class = models.CharField(max_length=2, choices=VariantClass.CHOICES, null=True, blank=True)
    vep_skipped_reason = models.CharField(max_length=1, choices=VEPSkippedReason.CHOICES, null=True, blank=True)

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

    @property
    def gnomad_variant(self):
        """ If you have gnomAD frequency - returns variant formatted as per 1-169519049-T-C """
        if self.gnomad_af is not None:
            v = self.variant
            return f"{v.locus.chrom}-{v.locus.position}-{v.locus.ref}-{v.alt}"
        return None

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
        return get_variant_search_terms(self.transcript_annotation)

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
    def get_overlapping_genes_q(gene_ids_qs) -> Q:
        """ Returns Q object where any transcript for a variant matches gene in genes
            Variants can overlap with multiple genes, and the VariantAnnotation (ie "pick" or representative annotation)
            may not be the one in the gene list. Thus we have to check everything that was in transcript annotation too
        """

        vto = VariantGeneOverlap.objects.filter(gene__in=gene_ids_qs)
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
    import_status = models.CharField(max_length=1, choices=ImportStatus.CHOICES, default=ImportStatus.CREATED)
    celery_task = models.CharField(max_length=36, null=True)

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
    error_message = models.TextField(blank=True, null=True)  # Set if any error...
    @property
    def unique_created_variants(self):
        return self.createdmanualvariant_set.distinct('variant_id').all()


class CreatedManualVariant(models.Model):
    manual_variant_entry = models.ForeignKey(ManualVariantEntry, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)


class InvalidAnnotationVersionError(Exception):
    pass


# Use this as a way to keep all the versions together.
# We re-use this, if nobody has referenced it - as you may eg update
# variant/gene/clinvar annotations every 6 months etc and no point having so many sub-versions
class AnnotationVersion(models.Model):
    SUB_ANNOTATIONS = ['variant_annotation_version', 'ensembl_gene_annotation_version', 'clinvar_version',
                       'human_protein_atlas_version']

    annotation_date = models.DateTimeField(auto_now=True)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, null=True, on_delete=PROTECT)
    ensembl_gene_annotation_version = models.ForeignKey(EnsemblGeneAnnotationVersion, null=True, on_delete=PROTECT)
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

    def validate(self):
        missing_sub_annotations = []
        for field in self.SUB_ANNOTATIONS:
            sub_annotation = getattr(self, field)
            if sub_annotation is None:
                missing_sub_annotations.append(field)
        if missing_sub_annotations:
            missing = ", ".join([str(s) for s in missing_sub_annotations])
            raise InvalidAnnotationVersionError(f"AnnotationVersion: {self} missing sub annotations: {missing}")

    @staticmethod
    def latest(genome_build: GenomeBuild, validate=True):
        av: AnnotationVersion = AnnotationVersion.objects.filter(genome_build=genome_build).order_by("annotation_date").last()
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

        def latest_for_build(klass, genome_build):
            qs = klass.objects.filter(genome_build=genome_build)
            return qs.order_by('annotation_date').last()

        if sub_version_genome_build:
            builds = [sub_version_genome_build]  # Just this one
        else:
            builds = GenomeBuild.builds_with_annotation()

        for genome_build in builds:
            kwargs = {
                "genome_build": genome_build,
                "variant_annotation_version": latest_for_build(VariantAnnotationVersion, genome_build),
                "ensembl_gene_annotation_version": latest_for_build(EnsemblGeneAnnotationVersion, genome_build),
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

    def get_ensembl_gene_annotation(self):
        return EnsemblGeneAnnotation.objects.filter(version=self.ensembl_gene_annotation_version)

    def get_clinvar(self):
        return ClinVar.objects.filter(version=self.clinvar_version)

    def get_human_protein_atlas_annotation(self):
        return HumanProteinAtlasAnnotation.objects.filter(version=self.human_protein_atlas_version)

    @property
    def long_description(self):
        sub_versions = [f"Variant: {self.variant_annotation_version}",
                        f"Gene: {self.ensembl_gene_annotation_version}",
                        f"ClinVar: {self.clinvar_version}",
                        f"HPA: {self.human_protein_atlas_version}"]
        sub_versions_str = ", ".join(sub_versions)
        return f"{self}. {sub_versions_str}"

    def __str__(self):
        return f"{self.pk} ({self.annotation_date.date()})"


class MonarchDiseaseOntology(models.Model):
    # PK = Mondo ID
    name = models.TextField()
    definition = models.TextField(null=True)

    # TODO: Mappings to all the other ontologies
    # phenotype_mim = models.ForeignKey(MIMMorbid, on_delete=CASCADE)

    @staticmethod
    def mondo_id_as_int(mondo_text) -> int:
        if isinstance(mondo_text, int):
            return mondo_text
        """ MONDO:0005045 -> 0005045 """
        return int(mondo_text.split(":")[1])

    @property
    def id_str(self) -> str:
        return MonarchDiseaseOntology.mondo_int_as_id(self.id)

    @staticmethod
    def mondo_int_as_id(mondo_id: int) -> str:
        num_part = str(mondo_id).rjust(7, '0')
        return f"MONDO:{num_part}"


class MonarchDiseaseOntologyGeneRelationship(models.Model):
    mondo = models.ForeignKey(MonarchDiseaseOntology, on_delete=CASCADE)
    relationship = models.TextField()
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)

    class Meta:
        unique_together = ('mondo', 'relationship', 'gene_symbol')


class MonarchDiseaseOntologyMIMMorbid(models.Model):
    mondo = models.ForeignKey(MonarchDiseaseOntology, on_delete=CASCADE)
    relationship = models.TextField()
    omim_id = models.IntegerField()  # could also link this to MIMMorbid records, but that gets out of date

    class Meta:
        unique_together = ('mondo', 'relationship', 'omim_id')


class CachedWebResource(TimeStampedModel):
    """ These are annotations that can populate themselves via the web

        * Can be manually updated via button kicking off annotation.views.load_cached_web_resource
        * List of names each system uses is in settings
        * Creation of one with name triggers signals which apps can configure to listen to / run populators
    """
    name = models.TextField(primary_key=True)
    description = models.TextField(blank=True)
    import_status = models.CharField(max_length=1, choices=ImportStatus.CHOICES, default=ImportStatus.CREATED)

    @staticmethod
    def named_handler_factory(name, celery_task_class):
        def handler(sender, instance, **kwargs):
            created = kwargs.get("created")

            if created and instance.name == name:
                task_signature = celery_task_class.si(instance.pk)
                task_signature.apply_async()

        return handler


class GeneDiseaseCurator(models.Model):
    name = models.TextField(primary_key=True)
    css_class = models.TextField(blank=True)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    cached_web_resource = models.ForeignKey(CachedWebResource, null=True, on_delete=CASCADE)  # If from web


class DiseaseValidity(models.Model):
    gene_disease_curator = models.ForeignKey(GeneDiseaseCurator, on_delete=CASCADE)
    text_phenotype = models.TextField(blank=True)
    sop = models.TextField(blank=True)  # ClinGen Gene Clinical Validity SOP
    mondo = models.ForeignKey(MonarchDiseaseOntology, null=True, on_delete=SET_NULL)
    hpo_synonym = models.ForeignKey(HPOSynonym, null=True, on_delete=SET_NULL)
    mim_morbid_alias = models.ForeignKey(MIMMorbidAlias, null=True, on_delete=SET_NULL)
    classification = models.CharField(max_length=1, choices=ClinGenClassification.CHOICES)
    date = models.DateField(null=True)
    validity_summary_url = models.TextField()


# Because there can be multiple Ensembl Genes for a symbol, link those via a FK
class GeneDiseaseValidity(models.Model):
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    disease_validity = models.ForeignKey(DiseaseValidity, on_delete=CASCADE)


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


class GeneDiseaseValidityEvidence(models.Model):
    gene_disease = models.ForeignKey(GeneDiseaseValidity, on_delete=CASCADE)
    citation = models.ForeignKey(Citation, null=True, on_delete=CASCADE)


class MutationalSignatureInfo(models.Model):
    signature_id = models.IntegerField(primary_key=True)
    cancer_types = models.TextField()
    proposed_aetiology = models.TextField()
    additional_mutational_features = models.TextField()
    comments = models.TextField()
