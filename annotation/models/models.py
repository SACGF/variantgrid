import logging
import operator
import os
import re
from collections import defaultdict
from datetime import datetime
from functools import cached_property, reduce
from typing import Optional, Callable, Iterable

from Bio.Data.IUPACData import protein_letters_1to3
from django.conf import settings
from django.contrib.auth.models import User
from django.contrib.postgres.fields import ArrayField
from django.core.exceptions import PermissionDenied
from django.db import models, transaction, connection
from django.db.models import QuerySet, Subquery, OuterRef
from django.db.models.deletion import PROTECT, CASCADE, SET_NULL
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.urls import reverse
from django_extensions.db.models import TimeStampedModel
from psqlextra.models import PostgresPartitionedModel
from psqlextra.types import PostgresPartitioningMethod

from annotation.external_search_terms import get_variant_search_terms, get_variant_pubmed_search_terms
from annotation.models.damage_enums import Polyphen2Prediction, FATHMMPrediction, MutationTasterPrediction, \
    SIFTPrediction, PathogenicityImpact, MutationAssessorPrediction, ALoFTPrediction
from annotation.models.models_citations import Citation, CitationFetchRequest, CitationFetchResponse
from annotation.models.models_enums import AnnotationStatus, \
    ColumnAnnotationCategory, VEPPlugin, VEPCustom, ClinVarReviewStatus, VEPSkippedReason, \
    ManualVariantEntryType, HumanProteinAtlasAbundance, EssentialGeneCRISPR, EssentialGeneCRISPR2, \
    EssentialGeneGeneTrap, VariantAnnotationPipelineType
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from classification.enums import AlleleOriginBucket
from genes.models import GeneSymbol, Gene, TranscriptVersion, Transcript, GeneAnnotationRelease
from genes.models_enums import AnnotationConsortium
from library.django_utils import object_is_referenced
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.genomics import parse_gnomad_coord
from library.genomics.vcf_enums import VariantClass
from library.utils import invert_dict, name_from_filename, first, all_equal
from ontology.models import OntologyVersion
from patients.models_enums import GnomADPopulation
from snpdb.models import GenomeBuild, Variant, VariantGridColumn, Q, VCF, DBSNP_PATTERN, VARIANT_PATTERN, \
    HGVS_UNCLEANED_PATTERN, Allele, VARIANT_SYMBOLIC_PATTERN
from snpdb.models.models_enums import ImportStatus


class SubVersionPartition(RelatedModelsPartitionModel):
    RECORDS_FK_FIELD_TO_THIS_MODEL = "version_id"
    PARTITION_LABEL_TEXT = "version"
    created = models.DateTimeField(auto_now_add=True)  # Inserted into DB
    annotation_date = models.DateTimeField(auto_now_add=True)  # Date of annotation (what we sort by)

    class Meta:
        abstract = True

    def save(self, *args, **kwargs):
        created = not self.pk
        super().save(*args, **kwargs)
        if created:
            genome_build = getattr(self, "genome_build", None)
            AnnotationVersion.new_sub_version(genome_build)

    def __str__(self):
        date_str = self.annotation_date.strftime("%d %B %Y")
        return f"v{self.pk}. ({date_str})"


class ClinVarVersion(SubVersionPartition):
    RECORDS_BASE_TABLE_NAMES = ["annotation_clinvar"]
    filename = models.TextField()
    sha256_hash = models.TextField()
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)

    @staticmethod
    def get_annotation_date_from_filename(filename) -> datetime:
        # file looks like: clinvar_20160302.vcf.gz
        base_name = os.path.basename(filename)
        # Allow a bit of lee-way as may want to call it e.g. clinvar_20190708_grch37.vcf.gz
        CLINVAR_PATTERN_37_38 = r"^clinvar_(\d{8}).*\.vcf"
        if m := re.match(CLINVAR_PATTERN_37_38, base_name):
            date_time = m.group(1)
            return datetime.strptime(date_time, "%Y%m%d")
        CLINVAR_PATTERN_T2T = r"^Homo_sapiens-GCA_009914755.4-(20\d{2})_(10)-clinvar.vcf.gz"
        if m := re.match(CLINVAR_PATTERN_T2T, base_name):
            year = m.group(1)
            month = m.group(0)
            return datetime(int(year), int(month), 1)  # Just go with 1st of month

        patterns = ", ".join([CLINVAR_PATTERN_37_38, CLINVAR_PATTERN_T2T])
        msg = f"File name '{base_name}' didn't match {patterns=}"
        raise ValueError(msg)

    def __str__(self):
        date_str = self.annotation_date.strftime("%d %B %Y")
        return f"v{self.pk}. {self.genome_build} ({date_str})"

@receiver(pre_delete, sender=ClinVarVersion)
def clinvar_version_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    instance.delete_related_objects()


class ClinVar(models.Model):

    class Meta:
        verbose_name = "ClinVar annotation"
        verbose_name_plural = "ClinVar annotations"

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

    ALLELE_ORIGIN_BUCKETS = {
        0: AlleleOriginBucket.UNKNOWN,
        1: AlleleOriginBucket.GERMLINE,
        2: AlleleOriginBucket.SOMATIC
    }

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
    preferred_disease_name = models.TextField(null=True, blank=True)
    disease_database_name = models.TextField(null=True, blank=True)
    review_status = models.CharField(max_length=1, null=True, choices=ClinVarReviewStatus.choices)
    clinical_significance = models.TextField(null=True, blank=True)

    # If clinical_significance = 'Conflicting_interpretations_of_pathogenicity'
    conflicting_clinical_significance = models.TextField(null=True, blank=True)
    highest_pathogenicity = models.IntegerField(default=0)  # Highest of clinical_significance
    clinical_sources = models.TextField(null=True, blank=True)
    origin = models.IntegerField(default=0)
    suspect_reason_code = models.IntegerField(default=0)
    drug_response = models.BooleanField(default=False)

    # ONCREVSTAT
    oncogenic_review_status = models.CharField(max_length=1, null=True, choices=ClinVarReviewStatus.choices)
    # ONCINCL
    oncogenic_classification = models.TextField(null=True, blank=True)
    # ONCCONF
    oncogenic_conflicting_classification = models.TextField(null=True, blank=True)

    # ONCDN
    oncogenic_preferred_disease_name = models.TextField(null=True, blank=True)
    # ONCDISDB
    oncogenic_disease_database_name = models.TextField(null=True, blank=True)

    # SCIREVSTAT
    somatic_review_status = models.CharField(max_length=1, null=True, choices=ClinVarReviewStatus.choices)
    # SCI
    somatic_clinical_significance = models.TextField(null=True, blank=True)

    # SCIDN
    somatic_preferred_disease_name = models.TextField(null=True, blank=True)
    # SCIDISDB
    somatic_disease_database_name = models.TextField(null=True, blank=True)

    @staticmethod
    def _database_terms(clinvar_db_value) -> list[str]:
        if db_name_text := clinvar_db_value:
            def fix_name(name: str):
                name = name.strip()
                if name.startswith("MONDO:MONDO:"):
                    name = name.replace("MONDO:MONDO:", "MONDO:")
                elif name.startswith("Orphanet:ORPHA"):
                    name = "ORPHA:" + name[len("Orphanet:ORPHA"):]
                elif name.startswith("Human_Phenotype_Ontology"):
                    name = name[25:]
                return name

            db_names = list(sorted(fix_name(db_name) for db_name in re.split("[|,]", db_name_text)))
            return db_names
        return []

    @property
    def clinvar_disease_database_terms(self) -> list[str]:
        """
        Deprecated, use germline_disease_database_terms
        """
        return self.germline_disease_database_terms

    @cached_property
    def germline_disease_database_terms(self) -> list[str]:
        return ClinVar._database_terms(self.disease_database_name)

    @cached_property
    def somatic_disease_database_terms(self) -> list[str]:
        return ClinVar._database_terms(self.somatic_disease_database_name)

    @cached_property
    def oncogenic_disease_database_terms(self) -> list[str]:
        return ClinVar._database_terms(self.oncogenic_disease_database_name)

    # repeat the above for oncogenic and somatic

    @property
    def clinical_sources_list(self) -> list[str]:
        if clinical_sources := self.clinical_sources:
            return [name.strip() for name in clinical_sources.split("|")]
        return []

    @property
    def stars(self):
        """
        deprecated - use .germline_stars
        """
        return self.germline_stars

    @staticmethod
    def _stars_for(value: Optional[str]) -> int:
        if not value:
            return 0
        return ClinVarReviewStatus(value).stars()

    @property
    def germline_stars(self) -> int:
        return ClinVar._stars_for(self.review_status)

    @property
    def somatic_stars(self) -> int:
        return ClinVar._stars_for(self.somatic_review_status)

    @property
    def oncogenic_stars(self) -> int:
        return ClinVar._stars_for(self.oncogenic_review_status)

    @property
    def is_expert_panel_or_greater(self):
        return max(self.germline_stars, self.somatic_stars, self.oncogenic_stars) >= CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE

    def get_origin_display(self):
        """
        FIXME use allele_origins clinvar_origin is a flag
        it might contain multiple values
        """
        return ClinVar.ALLELE_ORIGIN.get(self.origin)

    @property
    def get_allele_origins_display(self):
        return ", ".join(self.allele_origins)

    @cached_property
    def allele_origins(self) -> list[str]:
        origins: list[str] = []
        for flag, label in ClinVar.ALLELE_ORIGIN.items():
            if flag:
                if self.origin & flag:
                    origins.append(label)
        if not origins:
            origins.append("unknown")
        return origins

    @property
    def allele_origin_bucket(self) -> AlleleOriginBucket:
        is_germline = "germline" in self.allele_origins
        is_somatic = "somatic" in self.allele_origins
        if is_germline and is_somatic:
            return AlleleOriginBucket.UNKNOWN  # technically both, but unknown gets treated like both
        elif is_germline:
            return AlleleOriginBucket.GERMLINE
        elif is_somatic:
            return AlleleOriginBucket.SOMATIC
        else:
            return AlleleOriginBucket.UNKNOWN

    def get_suspect_reason_code_display(self):
        return ClinVar.SUSPECT_REASON_CODES.get(self.suspect_reason_code)

    def get_loaded_citations(self) -> CitationFetchResponse:
        cvc_qs = ClinVarCitation.objects.filter(clinvar_variation_id=self.clinvar_variation_id,
                                                clinvar_allele_id=self.clinvar_allele_id)
        return CitationFetchRequest.fetch_all_now(Citation.objects.filter(clinvarcitation__in=cvc_qs))

    @property
    def citation_ids(self) -> list[str]:
        return sorted(set(ClinVarCitation.objects.filter(clinvar_variation_id=self.clinvar_variation_id,
                                       clinvar_allele_id=self.clinvar_allele_id).values_list('citation_id', flat=True)))

    def json_summary(self) -> dict:
        return {
            "clinvar_variation_id": self.clinvar_variation_id,
            "ClinSig": self.clinical_significance,
            "stars": self.germline_stars,
            "disease": self.preferred_disease_name,
            "origin": self.get_origin_display(),
        }

    def short_summary(self) -> str:
        info = self.json_summary()
        return " ".join(f"{k}:{v}" for k,v in info.items())

    def __str__(self):
        return f"ClinVar: variant: {self.variant}, path: {self.highest_pathogenicity}"


class ClinVarRecordCollection(TimeStampedModel):
    """
    Stores data about when we've retrieved individual ClinVar records for a clinvar variation id.
    Importantly, when we did it, and what was the minimum number of stars on a record that we kept.
    Let's us know if we can re-use the cached ClinVarRecords or if we should retrieve them fresh from ClinVar.
    """

    class Meta:
        verbose_name = "ClinVar record collection"
        ordering = ["-max_stars", "-pk"]

    clinvar_variation_id = models.IntegerField(primary_key=True)
    allele = models.ForeignKey(Allele, null=True, on_delete=SET_NULL)
    urls = ArrayField(base_field=models.TextField(), blank=True, null=True, size=None)
    last_loaded = models.DateTimeField(blank=True, null=True)
    parser_version = models.IntegerField(blank=True, null=True)

    max_stars = models.IntegerField(blank=True, null=True)
    expert_panel = models.OneToOneField('ClinVarRecord', on_delete=SET_NULL, null=True, blank=True)

    def records_with_min_stars(self, min_stars: int) -> list['ClinVarRecord']:
        return list(sorted(self.clinvarrecord_set.filter(stars__gte=min_stars), reverse=True))

    def update_with_records_and_save(self, records: list['ClinVarRecord']):
        records = list(sorted(records, reverse=True))
        self.clinvarrecord_set.all().delete()
        for record in records:
            record.clinvar_record_collection = self
        ClinVarRecord.objects.bulk_create(records)
        self.expert_panel = None
        self.max_stars = None
        if best_record := first(records):
            self.max_stars = best_record.stars
            if best_record.is_expert_panel_or_greater:
                self.expert_panel = best_record
        self.save()


class ClinVarRecord(TimeStampedModel):
    """
    Represents a single record within ClinVar.
    Important to note this has been retrieved from ClinVar, and not our submission to ClinVar.
    """

    class Meta:
        verbose_name = "ClinVar record"
        ordering = ["-stars", "-date_last_evaluated"]

    clinvar_record_collection = models.ForeignKey(ClinVarRecordCollection, on_delete=CASCADE)
    record_id = models.TextField(primary_key=True)  # SCV
    stars = models.IntegerField()
    org_id = models.TextField()
    genome_build = models.TextField(null=True, blank=True)
    review_status = models.TextField()
    submitter = models.TextField()

    submitter_date = models.DateField(null=True, blank=True)
    date_last_evaluated = models.DateField(null=True, blank=True)
    date_clinvar_created = models.DateField(null=True, blank=True)
    date_clinvar_updated = models.DateField(null=True, blank=True)

    c_hgvs = models.TextField(null=True, blank=True)
    variant_coordinate = models.TextField(null=True, blank=True)
    condition = models.TextField(null=True, blank=True)
    clinical_significance = models.TextField(null=True, blank=True)
    somatic_clinical_significance = models.TextField(null=True, blank=True)
    gene_symbol = models.TextField(null=True, blank=True)
    interpretation_summary = models.TextField(null=True, blank=True)
    assertion_method = models.TextField(null=True, blank=True)
    allele_origin = models.TextField(null=True, blank=True)
    allele_origin_bucket = models.TextField(choices=AlleleOriginBucket.choices, null=True, blank=True)

    @property
    def conditions(self) -> list[str]:
        if condition := self.condition:
            if ":" in condition and ";" in condition:
                parts = [p.strip() for p in condition.split(";")]
                return parts
            else:
                return [condition]
        return []

    def mark_invalid(self):
        setattr(self, '_invalid', True)

    def __bool__(self) -> bool:
        if hasattr(self, '_invalid'):
            return not getattr(self, '_invalid')
        return True

    def __lt__(self, other):
        def sort_key(record: ClinVarRecord):
            return record.stars, record.date_last_evaluated or record.submitter_date
        return sort_key(self) < sort_key(other)

    def __str__(self):
        date_last_evaluated_str = ""
        if date_last_evaluated := self.date_last_evaluated or self.submitter_date:
            date_last_evaluated_str = date_last_evaluated.strftime('%Y-%m-%d')

        return f"{self.record_id} {self.get_allele_origin_bucket_display()}, {self.stars} stars, class={self.clinical_significance} somatic clin sig={self.somatic_clinical_significance}, condition=\"{self.condition}\" {date_last_evaluated_str}"

    @property
    def is_expert_panel_or_greater(self):
        return self.stars >= CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE


class ClinVarCitationsCollection(models.Model):
    cached_web_resource = models.ForeignKey('CachedWebResource', null=True, on_delete=CASCADE)


class ClinVarCitation(models.Model):
    clinvar_citations_collection = models.ForeignKey(ClinVarCitationsCollection, on_delete=CASCADE)
    clinvar_variation_id = models.IntegerField()
    clinvar_allele_id = models.IntegerField()
    citation = models.ForeignKey(Citation, on_delete=CASCADE)


class DBNSFPGeneAnnotationVersion(TimeStampedModel):
    """ @see https://sites.google.com/site/jpopgen/dbNSFP
        This isn't updated every release, so can have same hash across diff versions """
    version = models.TextField(primary_key=True)
    sha256_hash = models.TextField()

    class Meta:
        unique_together = ('version', 'sha256_hash')

    def save(self, *args, **kwargs):
        created = not self.pk
        super().save(*args, **kwargs)
        if created:
            logging.info("Creating new DBNSFPGeneAnnotation partition")
            version = self.pk
            connection.schema_editor().add_list_partition(
                model=DBNSFPGeneAnnotation,
                name=f"version_{version}",
                values=[version],
            )

    @staticmethod
    def latest() -> Optional['DBNSFPGeneAnnotationVersion']:
        return DBNSFPGeneAnnotationVersion.objects.order_by("created").last()


class DBNSFPGeneAnnotation(PostgresPartitionedModel, TimeStampedModel):
    """ @see https://sites.google.com/site/jpopgen/dbNSFP """
    version = models.ForeignKey(DBNSFPGeneAnnotationVersion, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE)
    refseq_transcript = models.ForeignKey(Transcript, related_name="refseq_dbnsfp_gene",
                                          null=True, blank=True, on_delete=CASCADE)
    ensembl_transcript = models.ForeignKey(Transcript, related_name="ensembl_dbnsfp_gene",
                                           null=True, blank=True, on_delete=CASCADE)
    gene_damage_index_score = models.FloatField(null=True)
    gene_damage_index_phred = models.FloatField(null=True)
    phi = models.FloatField(null=True)
    ghis = models.FloatField(null=True)
    prec = models.FloatField(null=True)
    hipred_score = models.FloatField(null=True)
    gnomad_pli = models.FloatField(null=True, blank=True)
    gnomad_prec = models.FloatField(null=True, blank=True)
    gnomad_pnull = models.FloatField(null=True, blank=True)
    loftool = models.FloatField(null=True, blank=True)
    gene_indispensability_score = models.FloatField(null=True, blank=True)
    hipred_prediction = models.BooleanField(null=True)
    gene_indispensability_pred = models.BooleanField(null=True)
    pathway_biocarta_full = models.TextField(null=True, blank=True)
    pathway_consensus_pathdb = models.TextField(null=True, blank=True)
    pathway_kegg_id = models.TextField(null=True, blank=True)
    pathway_kegg_full = models.TextField(null=True, blank=True)
    gwas_trait_association = models.TextField(null=True, blank=True)
    go_biological_process = models.TextField(null=True, blank=True)
    go_cellular_component = models.TextField(null=True, blank=True)
    go_molecular_function = models.TextField(null=True, blank=True)
    interactions_biogrid = models.TextField(null=True, blank=True)
    interactions_consensus_pathdb = models.TextField(null=True, blank=True)
    expression_egenetics = models.TextField(null=True, blank=True)
    expression_gnf_atlas = models.TextField(null=True, blank=True)
    essential_gene_crispr = models.CharField(max_length=1, null=True, blank=True,
                                             choices=EssentialGeneCRISPR.choices)
    essential_gene_crispr2 = models.CharField(max_length=1, null=True, blank=True,
                                              choices=EssentialGeneCRISPR2.choices)
    essential_gene_gene_trap = models.CharField(max_length=1, null=True, blank=True,
                                                choices=EssentialGeneGeneTrap.choices)

    def get_gene_indispensability_pred_display(self) -> str:
        if self.gene_indispensability_pred is None:
            return "n/a"
        elif self.gene_indispensability_pred:
            return "Essential"
        return "LoF tolerant"

    class PartitioningMeta:
        method = PostgresPartitioningMethod.LIST
        key = ["version_id"]


class GeneAnnotationVersion(SubVersionPartition):
    RECORDS_BASE_TABLE_NAMES = ["annotation_geneannotation"]
    gene_annotation_release = models.ForeignKey(GeneAnnotationRelease, on_delete=CASCADE)
    ontology_version = models.ForeignKey(OntologyVersion, null=True, on_delete=PROTECT)
    dbnsfp_gene_version = models.ForeignKey(DBNSFPGeneAnnotationVersion, null=True, on_delete=CASCADE)
    gnomad_import_date = models.DateTimeField()

    @property
    def genome_build(self):
        return self.gene_annotation_release.genome_build

    def __str__(self):
        return f"GAR: {self.gene_annotation_release}, Ont: {self.ontology_version}, DBnsFP: {self.dbnsfp_gene_version}"


@receiver(pre_delete, sender=GeneAnnotationVersion)
def gene_annotation_version_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    instance.delete_related_objects()


class GeneAnnotation(models.Model):
    """ This is generated against genes found via GeneAnnotationRelease
        so that data matches up in analyses """
    version = models.ForeignKey(GeneAnnotationVersion, on_delete=CASCADE)
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    dbnsfp_gene = models.ForeignKey(DBNSFPGeneAnnotation, null=True, on_delete=SET_NULL)
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
    sha256_hash = models.TextField()
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
    column = models.TextField(unique=True)  # TODO: Do we actually need this?
    variant_grid_column = models.ForeignKey(VariantGridColumn, on_delete=CASCADE)  # VariantAnnotation.column dest
    genome_build = models.ForeignKey(GenomeBuild, null=True, on_delete=CASCADE)  # null = all builds
    pipeline_type = models.CharField(max_length=1, choices=VariantAnnotationPipelineType.choices,
                                     null=True, blank=True)  # Null = all pipeline types
    category = models.CharField(max_length=1, choices=ColumnAnnotationCategory.choices)
    source_field = models.TextField(null=True)  # @see use vep_info_field
    source_field_processing_description = models.TextField(null=True)
    vep_plugin = models.CharField(max_length=1, choices=VEPPlugin.choices, null=True)
    vep_custom = models.CharField(max_length=1, choices=VEPCustom.choices, null=True)
    source_field_has_custom_prefix = models.BooleanField(default=False)
    # We can use these min/max versions to turn on/off columns over time
    min_columns_version = models.IntegerField(null=True)
    max_columns_version = models.IntegerField(null=True)
    min_vep_version = models.IntegerField(null=True)
    max_vep_version = models.IntegerField(null=True)
    summary_stats = models.TextField(blank=True, null=True)  # Only used VEP 111 and on...

    def __str__(self) -> str:
        return self.column

    @property
    def vep_info_field(self):
        """ For VCFs, be sure to set source_field_has_custom_prefix=True
            Annotating with a VCF (short name = TopMed) brings in the DB column as "TopMed" and
            prefixes INFO fields, e.g. 'TOPMED' => TopMed_TOPMED.
            We need to adjust for this in BulkVEPVCFAnnotationInserter """

        vif = self.source_field
        if self.vep_custom and self.source_field_has_custom_prefix:
            if vif:
                vif = self.get_vep_custom_display() + "_" + vif
            else:
                vif = self.get_vep_custom_display()  # Just the prefix - used eg to return ID in VEP custom VCFs
        return vif

    @cached_property
    def columns_version_description(self) -> str:
        limits = []
        if self.min_columns_version:
            limits.append(f"column version >= {self.min_columns_version}")
        if self.max_columns_version:
            limits.append(f"column version <= {self.max_columns_version}")
        return " and ".join(limits)

    @staticmethod
    def get_columns_version_q(columns_version: int) -> Q:
        q_min = Q(min_columns_version__isnull=True) | Q(min_columns_version__lte=columns_version)
        q_max = Q(max_columns_version__isnull=True) | Q(max_columns_version__gte=columns_version)
        return q_min & q_max

    @staticmethod
    def get_vep_version_q(vep_version: int) -> Q:
        q_min = Q(min_vep_version__isnull=True) | Q(min_vep_version__lte=vep_version)
        q_max = Q(max_vep_version__isnull=True) | Q(max_vep_version__gte=vep_version)
        return q_min & q_max

    @staticmethod
    def get_pipeline_type_q(pipeline_type) -> Q:
        return Q(pipeline_type__isnull=True) | Q(pipeline_type=pipeline_type)

    @staticmethod
    def get_genome_build_q(genome_build) -> Q:
        return Q(genome_build=genome_build) | Q(genome_build__isnull=True)

    @staticmethod
    def get_q(genome_build: GenomeBuild, vep_version, columns_version, pipeline_type) -> Q:
        filters = [
            ColumnVEPField.get_genome_build_q(genome_build),
            ColumnVEPField.get_columns_version_q(columns_version),
            ColumnVEPField.get_vep_version_q(vep_version),
            ColumnVEPField.get_pipeline_type_q(pipeline_type),
        ]
        return reduce(operator.and_, filters)

    @staticmethod
    def filter(genome_build: GenomeBuild, vep_version: int, columns_version: int, pipeline_type):
        q = ColumnVEPField.get_q(genome_build, vep_version, columns_version, pipeline_type)
        return ColumnVEPField.objects.filter(q)

    @staticmethod
    def filter_for_build(genome_build: GenomeBuild):
        """ genome_build = NULL (no build) or matches provided build """
        return ColumnVEPField.objects.filter(ColumnVEPField.get_genome_build_q(genome_build))

    @staticmethod
    def get(genome_build: GenomeBuild, *columnvepfield_args, **columnvepfield_kwargs) -> QuerySet['ColumnVEPField']:
        qs = ColumnVEPField.filter_for_build(genome_build)
        qs = qs.filter(*columnvepfield_args, **columnvepfield_kwargs)
        return qs.distinct("source_field").order_by("source_field")

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
    genebuild = models.TextField()
    gnomad = models.TextField(blank=True, null=True)  # 37/38 only
    refseq = models.TextField(blank=True)
    regbuild = models.TextField(blank=True, null=True)  # 37/38 only
    sift = models.TextField(blank=True, null=True)  # 37/38 only
    dbnsfp = models.TextField(blank=True, null=True)  # 37/38 only
    distance = models.IntegerField(default=5000)  # VEP --distance parameter

    @property
    def gnomad_major_version(self) -> int:
        return int(self.gnomad.split(".", maxsplit=1)[0])

    @staticmethod
    def latest(genome_build, active=True) -> 'VariantAnnotationVersion':
        qs = VariantAnnotationVersion.objects.filter(genome_build=genome_build, active=active)
        return qs.order_by("annotation_date").last()

    @staticmethod
    def latest_for_all_builds(active=True) -> QuerySet:
        # Default to latest
        latest_qs = VariantAnnotationVersion.objects.filter(
            genome_build=OuterRef('genome_build'),
            active=active
        ).order_by('-annotation_date')

        return VariantAnnotationVersion.objects.filter(
            pk__in=Subquery(latest_qs.values('pk')[:1]),
            genome_build__in=GenomeBuild.builds_with_annotation()
        )

    def get_any_annotation_version(self):
        """ Often you don't care what annotation version you use, only that variant annotation version is this one """
        return self.annotationversion_set.last()

    def get_pathogenic_prediction_funcs(self) -> dict[str, Callable]:
        if self.columns_version == 1:
            return {
                'sift': lambda d: d in SIFTPrediction.get_damage_or_greater_levels(),
                'fathmm_pred_most_damaging': lambda d: d in FATHMMPrediction.get_damage_or_greater_levels(),
                'mutation_assessor_pred_most_damaging': lambda d: d in MutationAssessorPrediction.get_damage_or_greater_levels(),
                'mutation_taster_pred_most_damaging': lambda d: d in MutationTasterPrediction.get_damage_or_greater_levels(),
                'polyphen2_hvar_pred_most_damaging': lambda d: d in Polyphen2Prediction.get_damage_or_greater_levels(),
            }
        elif self.columns_version in (2, 3):
            pathogenic_rankscore = settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE
            pathogenic_prediction_columns = ['bayesdel_noaf_rankscore', 'cadd_raw_rankscore', 'clinpred_rankscore',
                                             'revel_rankscore', 'metalr_rankscore', 'vest4_rankscore']
            if self.columns_version == 3:
                pathogenic_prediction_columns.append("alphamissense_rankscore")

            return {c: lambda d: float(d) >= pathogenic_rankscore for c in pathogenic_prediction_columns}

        raise ValueError(f"Don't know fields for {self.columns_version=}")

    @cached_property
    def damage_predictions_description(self) -> str:
        pathogenic_prediction = list(self.get_pathogenic_prediction_funcs())
        columns = ", ".join(pathogenic_prediction)
        description = ""
        if self.columns_version == 1:
            description = f"Count of {columns} at the most damaging level."
        elif self.columns_version >= 2:
            description = f"Count of {columns} that exceed {settings.ANNOTATION_MIN_PATHOGENIC_RANKSCORE}"
        return description

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

    @cached_property
    def _gene_annotation_release_and_gff_url(self) -> tuple[Optional[str], Optional[str]]:
        release = None
        gff_url = None
        # See issue: https://github.com/Ensembl/ensembl-vep/issues/833 - perhaps this can be done more easily now
        if self.annotation_consortium == AnnotationConsortium.ENSEMBL:
            if self.genome_build.name == "GRCh37":
                # This has remained unchanged (last checked v108)
                release = "87"
                gff_url = "ftp://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz"
            elif self.genome_build.name == "GRCh38":
                release = str(self.vep)
                gff_url = f"ftp://ftp.ensembl.org/pub/release-{self.vep}/gff3/homo_sapiens/Homo_sapiens.GRCh38.{self.vep}.gff3.gz"
        else:
            if self.genome_build.name == "GRCh37":
                # This is the last GRCh37 release (checked VEP v108)
                if self.refseq == '2020-10-26 17:03:42 - GCF_000001405.25_GRCh37.p13_genomic.gff':
                    release = "105.20201022"
                    gff_url = f"http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/{release}/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz"
            elif self.genome_build.name == "GRCh38":
                if m := re.match(r"(109.20\d{6}) - GCF_000001405.39_GRCh38.p13_genomic.gff", self.refseq):
                    release = m.group(1)
                    gff_url = f"http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/{release}/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz"
                else:
                    (release, gff_filename) = self.refseq.split(" - ")
                    if not gff_filename.endswith(".gz"):
                        gff_filename += ".gz"
                    # This is good for VEP v108 (will need to keep on top of this)
                    patch_version = 14
                    gff_url = f"https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/{release}/GCF_000001405.40_GRCh38.p{patch_version}/{gff_filename}"
        return release, gff_url

    @property
    def gene_annotation_release_gff_url(self):
        return self._gene_annotation_release_and_gff_url[1]

    @cached_property
    def cdot_gene_release_filename(self) -> str:
        """ returns blank if unknown """
        name_components = []
        _release, gff_url = self._gene_annotation_release_and_gff_url
        if gff_url:
            name_components = [name_from_filename(gff_url)]
            # These ones got renamed as the filename wasn't unique
            if self.annotation_consortium == AnnotationConsortium.REFSEQ:
                if self.genome_build.name == "GRCh37":
                    if m := re.match(r"http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/(105.20\d{6})/GCF_000001405.25_GRCh37.p13/(GCF_000001405.25_GRCh37.p13_genomic).gff.gz", gff_url):
                        name_components = [m.group(2), m.group(1), "gff"]
                elif self.genome_build.name == "GRCh38":
                    # GCF_000001405.39_GRCh38.p13_genomic.109.20210514.gff.json.gz
                    if m := re.match(r"http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/(109.20\d{6})/GCF_000001405.39_GRCh38.p13/(GCF_000001405.39_GRCh38.p13_genomic).gff.gz", gff_url):
                        name_components = [m.group(2), m.group(1), "gff"]
            name_components.append("json.gz")

        return ".".join(name_components)

    @property
    def suggested_gene_annotation_release_name(self) -> str:
        release, _ = self._gene_annotation_release_and_gff_url
        if release:
            return f"{self.get_annotation_consortium_display()}_{release}"
        return "TOOD_RELEASE_NAME"

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
    dump_count = models.IntegerField(null=True)
    annotated_count = models.IntegerField(null=True)
    celery_task_logs = models.JSONField(null=False, default=dict)  # Key=task_id, so we keep logs from multiple runs

    class Meta:
        unique_together = ('annotation_range_lock', 'pipeline_type')

    def get_absolute_url(self):
        return reverse('view_annotation_run', kwargs={'annotation_run_id': self.pk})

    @staticmethod
    def count_not_successful_runs_for_version(version: VariantAnnotationVersion) -> int:
        qs = AnnotationRun.objects.filter(annotation_range_lock__version=version)
        return qs.exclude(status=AnnotationStatus.FINISHED).count()

    @staticmethod
    def get_for_variant(variant: Variant, genome_build) -> Optional['AnnotationRun']:
        if variant.is_symbolic:
            pipeline_type = VariantAnnotationPipelineType.STRUCTURAL_VARIANT
        else:
            pipeline_type = VariantAnnotationPipelineType.STANDARD
        # For newly created variants, there will only be one per build for latest annotation version
        ar: Optional[AnnotationRun]
        ar = AnnotationRun.objects.filter(annotation_range_lock__version__genome_build=genome_build,
                                          annotation_range_lock__min_variant__lte=variant.pk,
                                          annotation_range_lock__max_variant__gte=variant.pk,
                                          pipeline_type=pipeline_type).first()
        return ar

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

    def get_dump_filename(self) -> str:
        PIPELINE_TYPE = {
            VariantAnnotationPipelineType.STANDARD: "standard",
            VariantAnnotationPipelineType.STRUCTURAL_VARIANT: "structural_variant",
        }
        type_desc = PIPELINE_TYPE.get(self.pipeline_type, str(self.pipeline_type))
        vcf_base_name = f"dump_{self.pk}_{type_desc}.vcf"
        return os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, vcf_base_name)

    def delete(self, using=None, keep_parents=False):
        self.delete_related_objects()
        super().delete(using=using, keep_parents=keep_parents)

    def set_task_log(self, key, value):
        assert self.task_id is not None
        task_log = self.celery_task_logs.get(self.task_id, {})
        task_log[key] = value

    def __str__(self):
        return f"AnnotationRun: {self.pk}/{self.get_pipeline_type_display()}: ({self.get_status_display()})"


class AbstractVariantAnnotation(models.Model):
    """ Common fields between VariantAnnotation and VariantTranscriptAnnotation
        These fields are PER-TRANSCRIPT """
    SV_HGVS_TOO_LONG_MESSAGE = "HGVS not calculated due to length"
    SV_HGVS_ERROR_MESSAGE = "Error creating HGVS"

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
    # TODO: This doesn't need to be nullable (default=False) - but will be slow. Change with next schema change
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
    gerp_pp_rs = models.FloatField(null=True, blank=True)  # Genomic Evolutionary Rate Profiling (GERP++)
    grantham = models.IntegerField(null=True, blank=True)
    hgvs_c = models.TextField(null=True, blank=True)
    hgvs_p = models.TextField(null=True, blank=True)
    impact = models.CharField(max_length=1, choices=PathogenicityImpact.CHOICES, null=True, blank=True)
    interpro_domain = models.TextField(null=True, blank=True)
    intron = models.TextField(null=True, blank=True)
    maxentscan_alt = models.FloatField(null=True, blank=True)
    maxentscan_diff = models.FloatField(null=True, blank=True)
    maxentscan_ref = models.FloatField(null=True, blank=True)
    maxentscan_percent_diff_ref = models.FloatField(null=True, blank=True)
    mutation_assessor_pred_most_damaging = models.CharField(max_length=1, choices=MutationAssessorPrediction.CHOICES, null=True, blank=True)
    mutation_taster_pred_most_damaging = models.CharField(max_length=1, choices=MutationTasterPrediction.CHOICES, null=True, blank=True)
    polyphen2_hvar_pred_most_damaging = models.CharField(max_length=1, choices=Polyphen2Prediction.CHOICES, null=True, blank=True)
    # protein_position = text as it can be e.g. indel: "22-23" or splicing: "?-10" or "10-?"
    protein_position = models.TextField(null=True, blank=True)
    revel_score = models.FloatField(null=True, blank=True)
    sift = models.CharField(max_length=1, choices=SIFTPrediction.CHOICES, null=True, blank=True)
    splice_region = models.TextField(null=True, blank=True)
    symbol = models.TextField(null=True, blank=True)

    mavedb_score = models.FloatField(null=True, blank=True)
    mavedb_urn = models.TextField(null=True, blank=True)

    class Meta:
        abstract = True

    @property
    def has_hgvs_c(self) -> bool:
        return self.hgvs_c and self.hgvs_c != self.SV_HGVS_TOO_LONG_MESSAGE

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
        """ p.Ala2351Thr -> p.A2351T """

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

    def get_hgvs_c_with_symbol(self) -> str:
        hgvs_c = self.hgvs_c
        if self.has_hgvs_c and self.symbol:
            from genes.hgvs import HGVSMatcher
            hgvs_matcher = HGVSMatcher(self.version.genome_build)
            hgvs_variant = hgvs_matcher.create_hgvs_variant(hgvs_c)
            hgvs_variant.gene = self.symbol
            hgvs_c = str(hgvs_variant)
        return hgvs_c


class VariantAnnotation(AbstractVariantAnnotation):
    """ This is the "representative transcript" chosen (1 per variant/annotation version) """
    GENE_COLUMN = "variantannotation__gene"

    # Only need this once per variant
    hgvs_g = models.TextField(null=True, blank=True)

    # Population frequency
    af_1kg = models.FloatField(null=True, blank=True)
    af_uk10k = models.FloatField(null=True, blank=True)
    topmed_af = models.FloatField(null=True, blank=True)
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
    gnomad_mid_af = models.FloatField(null=True, blank=True)  # Middle East - added in gnomADv4
    gnomad_nfe_af = models.FloatField(null=True, blank=True)
    gnomad_oth_af = models.FloatField(null=True, blank=True)
    gnomad_sas_af = models.FloatField(null=True, blank=True)
    # filtering allele frequencies (new in gnomADv4)
    gnomad_faf95 = models.FloatField(null=True, blank=True)
    gnomad_faf99 = models.FloatField(null=True, blank=True)
    gnomad_fafmax_faf95_max = models.FloatField(null=True, blank=True)
    gnomad_fafmax_faf99_max = models.FloatField(null=True, blank=True)
    gnomad_xy_af = models.FloatField(null=True, blank=True)
    gnomad_xy_ac = models.IntegerField(null=True, blank=True)
    gnomad_xy_an = models.IntegerField(null=True, blank=True)
    gnomad_hemi_count = models.IntegerField(null=True, blank=True)  # This is set from gnomad_xy_ac if gnomad_non_par
    gnomad_popmax_af = models.FloatField(null=True, blank=True)
    gnomad_popmax_ac = models.IntegerField(null=True, blank=True)
    gnomad_popmax_an = models.IntegerField(null=True, blank=True)
    gnomad_popmax_hom_alt = models.IntegerField(null=True, blank=True)
    gnomad_filtered = models.BooleanField(null=True, blank=True)
    gnomad_non_par = models.BooleanField(null=True, blank=True)  # Not pseudoautosomal regions
    gnomad_popmax = models.CharField(max_length=3, choices=GnomADPopulation.choices, null=True, blank=True)

    # These are populated for SVs using VEP plugin StructuralVariantOverlap
    # They are text as they can have multiple entries joined via '&'
    gnomad_sv_overlap_af = models.TextField(null=True, blank=True)
    gnomad_sv_overlap_percent = models.TextField(null=True, blank=True)
    gnomad_sv_overlap_name = models.TextField(null=True, blank=True)
    gnomad_sv_overlap_coords = models.TextField(null=True, blank=True)
    # Can't use filters unfortunately due to VEP custom bug, @see https://github.com/Ensembl/ensembl-vep/issues/1646
    # gnomad_sv_overlap_filters = models.TextField(null=True, blank=True)

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
    alphamissense_rankscore = models.FloatField(null=True, blank=True)

    # ALoFT (from dbNSFP)
    aloft_prob_tolerant = models.FloatField(null=True, blank=True)
    aloft_prob_recessive = models.FloatField(null=True, blank=True)
    aloft_prob_dominant = models.FloatField(null=True, blank=True)
    aloft_pred = models.CharField(max_length=1, choices=ALoFTPrediction.choices, null=True, blank=True)
    aloft_high_confidence = models.BooleanField(null=True, blank=True)
    aloft_ensembl_transcript = models.TextField(null=True, blank=True)  # Transcript of most damaging prediction chosen

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
        GnomADPopulation.MIDDLE_EASTERN: 'gnomad_mid_af',
        GnomADPopulation.NON_FINNISH_EUROPEAN: 'gnomad_nfe_af',
        GnomADPopulation.OTHER: 'gnomad_oth_af',
        GnomADPopulation.SOUTH_ASIAN: 'gnomad_sas_af',
    }

    GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS = [
        'gnomad_sv_overlap_af',
        'gnomad_sv_overlap_name',
        'gnomad_sv_overlap_percent',
        'gnomad_sv_overlap_coords',
    ]

    ALOFT_FIELDS = {
        "aloft_pred": "Pred",
        "aloft_high_confidence": "High Confidence",
        "aloft_prob_tolerant": "Tol",
        "aloft_prob_recessive": "Rec",
        "aloft_prob_dominant": "Dom",
        "aloft_ensembl_transcript": "Transcript",
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

    CONSERVATION_SCORES = {
        "gerp_pp_rs": {
            # UCSC says RS scores range from a maximum of 6.18 down to a below-zero minimum, which we cap at -12.36
            "min": -12.36,
            "max": 6.18,
        },
        # BigWig stats obtained via kent-335 bigWigInfo
        "phylop_30_way_mammalian": {
            "mean": 0.097833,
            "min": -20.0,
            "max": 1.312,
            "std": 0.727453,
        },
        "phylop_46_way_mammalian": {
            "mean": 0.035934,
            "min": -13.796,
            "max": 2.941,
            "std": 0.779426,
        },
        "phylop_100_way_vertebrate": {
            "mean": 0.093059,
            "min": -20.0,
            "max": 10.003,
            "std": 1.036944,
        },
        "phastcons_30_way_mammalian": {
            "mean": 0.128025,
            "min": 0.0,
            "max": 1.0,
            "std": 0.247422,
        },
        "phastcons_46_way_mammalian": {
            "mean": 0.088576,
            "min": 0.0,
            "max": 1.0,
            "std": 0.210242,
        },
        "phastcons_100_way_vertebrate": {
            "mean": 0.101765,
            "min": 0.0,
            "max": 1.0,
            "std": 0.237072,
        }
    }

    # List of filters to describe variants that can be annotated
    VARIANT_ANNOTATION_Q = [
        Variant.get_no_reference_q(),
        ~Q(alt__seq__in=['.', '*']),  # Exclude non-standard variants
    ]

    @cached_property
    def is_standard_annotation(self) -> bool:
        return self.annotation_run.pipeline_type == VariantAnnotationPipelineType.STANDARD

    @property
    def has_pathogenicity(self) -> bool:
        return self.is_standard_annotation

    @property
    def has_conservation(self) -> bool:
        """ Thanks to summary stats we can now do this in VEP112 """
        return self.is_standard_annotation or self.version.vep >= 112

    @property
    def has_splicing(self) -> bool:
        return self.is_standard_annotation

    @property
    def has_gnomad(self) -> bool:
        return bool(self.gnomad_af or self.gnomad2_liftover_af)

    @property
    def has_non_gnomad_population_frequency(self) -> bool:
        return self.is_standard_annotation

    @property
    def has_cosmic(self) -> bool:
        return self.is_standard_annotation

    @property
    def has_mavedb(self) -> bool:
        return self.is_standard_annotation and self.version.columns_version >= 3

    @property
    def has_gnomad_faf(self) -> bool:
        return self.is_standard_annotation and self.gnomad4_or_later

    @cached_property
    def has_extended_gnomad_fields(self):
        """ I grabbed a few new fields but haven't patched back to GRCh37 yet
            TODO: remove this and if statements in variant_details.html once issue #231 is completed """
        extended_fields = ["gnomad2_liftover_af", "gnomad_ac", "gnomad_an", "gnomad_popmax_ac",
                           "gnomad_popmax_an", "gnomad_popmax_hom_alt"]
        return any(getattr(self, f) is not None for f in extended_fields)

    @property
    def gnomad4_or_later(self) -> bool:
        return self.version.gnomad_major_version >= 4

    @property
    def has_hemi(self):
        return self.gnomad4_or_later and self.variant.locus.contig.name == 'X'

    @property
    def has_mid(self):
        return self.gnomad4_or_later

    @property
    def gnomad_url(self):
        GNOMAD2 = "gnomad_r2_1"
        GNOMAD3 = "gnomad_r3"
        GNOMAD4 = "gnomad_r4"

        gnomad_dataset = None
        if self.gnomad_af and self.version.gnomad.startswith("2.1"):
            gnomad_dataset = GNOMAD2
        elif self.version.gnomad.startswith("3"):
            if self.gnomad_af is not None:
                gnomad_dataset = GNOMAD3
        elif self.gnomad_af and self.version.gnomad.startswith("4"):
            gnomad_dataset = GNOMAD4

        url = None
        if gnomad_dataset:
            v = self.variant
            gnomad_variant = f"{v.locus.chrom}-{v.locus.position}-{v.locus.ref}-{v.alt}"
            url = self.get_gnomad_url(gnomad_variant, gnomad_dataset)
        return url

    @staticmethod
    def get_gnomad_url(gnomad_variant: str, dataset: str) -> str:
        return f"http://gnomad.broadinstitute.org/variant/{gnomad_variant}?dataset={dataset}"

    @staticmethod
    def mavedb_urn_to_urls(mavedb_urn: str) -> dict[str, str]:
        experiment_sets = set()
        EXPERIMENT_URL_PREFIX = "https://www.mavedb.org/#/experiment-sets/"
        if mavedb_urn:
            for urn in mavedb_urn.split("&"):
                es = urn.rsplit("-", 2)[0]
                experiment_sets.add(es)
        return {urn: EXPERIMENT_URL_PREFIX + urn for urn in sorted(experiment_sets)}

    def get_mave_urls(self) -> dict[str, str]:
        """ dict of label: url
            MAVE team intend on creating a variant landing page, so we'll replace this with that when ready
        """
        return self.mavedb_urn_to_urls(self.mavedb_urn)

    @property
    def mmid3_mastermind_urls(self) -> dict:
        return self.get_mmid3_mastermind_urls(self.mastermind_mmid3)

    @staticmethod
    def get_mmid3_mastermind_urls(mastermind_mmid3) -> dict:
        mmid3_mastermind_urls = {}
        if mastermind_mmid3:
            for mmid3 in mastermind_mmid3.split("&"):
                mmid3_mastermind_urls[mmid3] = f"https://mastermind.genomenon.com/detail?mutation={mmid3}"
        return mmid3_mastermind_urls

    def has_spliceai(self):
        return any((self.spliceai_pred_ds_ag, self.spliceai_pred_ds_al,
                    self.spliceai_pred_ds_dg, self.spliceai_pred_ds_dl))

    def highest_spliceai(self) -> int|None:
        values = []
        for (ds, _) in self.SPLICEAI_DS_DP.values():
            if v := getattr(self, ds):
                values.append(v)
        if values:
            return max(values)
        return None

    @staticmethod
    def get_gnomad_population_field(population):
        return VariantAnnotation.GNOMAD_FIELDS.get(population)

    @cached_property
    def transcript_annotation(self) -> list['VariantTranscriptAnnotation']:
        return self.variant.varianttranscriptannotation_set.filter(version=self.version)

    def get_search_terms(self):
        extra_terms = []
        if self.dbsnp_rs_id:
            extra_terms.append(self.dbsnp_rs_id)
        return get_variant_search_terms(self.transcript_annotation, extra_terms=extra_terms)

    def get_pubmed_search_terms(self):
        return get_variant_pubmed_search_terms(self.transcript_annotation)

    @staticmethod
    def vep_multi_fields_to_list_of_dicts(data: dict, fields: Iterable[str]) -> list[dict]:
        list_values = {}
        for field in fields:
            if value := data.get(field):
                list_values[field] = value.split("&")

        # Ensure they are all same length
        first_field_len = len(next(iter(list_values.values())))
        field_lengths = {k: len(v) for k, v in list_values.items()}

        if not all_equal(field_lengths.values()):
            logging.error("All split VEP multi-values must be equal length")
            for f, l in field_lengths.items():
                logging.error("%s: %d", f, l)
            raise ValueError(
                f"All split multi-values must be equal length: {field_lengths}")

        records: list[dict] = []
        for i in range(first_field_len):
            values_dict = {}
            for field in fields:
                if v := list_values.get(field):
                    values_dict[field] = v[i]
            records.append(values_dict)
        return records

    @property
    def gnomad_sv_overlap(self) -> list[dict]:
        sv_overlap_list = []
        for sv_overlap in self.vep_multi_fields_to_list_of_dicts(self.__dict__,
                                                      VariantAnnotation.GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS):
            sv_overlap['gnomad_sv_overlap_af'] = float(sv_overlap['gnomad_sv_overlap_af'])

            dataset = None
            remove_prefix = None
            if self.version.gnomad.startswith("2.1"):
                dataset = "gnomad_sv_r2_1"
                remove_prefix = "gnomAD-SV_v2.1_"
            elif self.version.gnomad.startswith("4"):
                dataset = "gnomad_sv_r4"
                remove_prefix = "gnomAD-SV_v3_"

            gnomad_variant = sv_overlap["gnomad_sv_overlap_name"]
            if remove_prefix:
                gnomad_variant = gnomad_variant.replace(remove_prefix, "")

            sv_overlap['gnomad_sv_overlap_url'] = self.get_gnomad_url(gnomad_variant, dataset)
            if coords := sv_overlap.get("gnomad_sv_overlap_coords"):
                _chrom, start, end = parse_gnomad_coord(coords)
                sv_overlap['gnomad_sv_overlap_length'] = end - start

            sv_overlap_list.append(sv_overlap)
        return sv_overlap_list

    @staticmethod
    def get_hgvs_g(variant: Variant) -> Optional[str]:
        """ This is very slow - only use on a page with just 1 variant, not in a loop """
        # They should all be the same so any will be ok - but take latest in case we fixed HGVS
        qs = variant.variantannotation_set.filter(hgvs_g__isnull=False).order_by("-version")
        data = qs.values_list("hgvs_g", flat=True)[:1]
        hgvs_g = None
        if data:
            hgvs_g = data[0]
        if hgvs_g is None:
            # Reference variants have no annotation - so we'll have to fall back to generating it
            from genes.hgvs import HGVSMatcher
            matcher = HGVSMatcher(variant.any_genome_build)
            hgvs_g = matcher.variant_to_g_hgvs(variant)

        return hgvs_g

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

    @cached_property
    def first_entry(self) -> Optional['ManualVariantEntry']:
        """ Often a user will enter a single variant into search and then wait for it to be created/annotated
            We have the "first" one to deal with this """
        return self.manualvariantentry_set.all().order_by("pk").first()

    @cached_property
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
                annotation_run = AnnotationRun.get_for_variant(variant, mve.genome_build)
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
    warning_message = models.TextField(blank=True, null=True)
    error_message = models.TextField(blank=True, null=True)  # Set if any error...

    @property
    def unique_created_variants(self):
        return self.createdmanualvariant_set.distinct('variant_id').all()

    @property
    def genome_build(self) -> GenomeBuild:
        return self.manual_variant_entry_collection.genome_build

    @staticmethod
    def get_entry_type(line: str) -> ManualVariantEntryType:
        MATCHERS = {
            DBSNP_PATTERN: ManualVariantEntryType.DBSNP,
            HGVS_UNCLEANED_PATTERN: ManualVariantEntryType.HGVS,
            VARIANT_PATTERN: ManualVariantEntryType.VARIANT,
            VARIANT_SYMBOLIC_PATTERN: ManualVariantEntryType.VARIANT,
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
# We re-use this, if nobody has referenced it - as you may e.g. update
# variant/gene/clinvar annotations every 6 months etc and no point having so many sub-versions
class AnnotationVersion(models.Model):
    objects = AnnotationVersionManager()  # Always select_related
    annotation_date = models.DateTimeField(auto_now=True)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    variant_annotation_version = models.ForeignKey(VariantAnnotationVersion, null=True, on_delete=PROTECT)
    gene_annotation_version = models.ForeignKey(GeneAnnotationVersion, null=True, on_delete=PROTECT)
    clinvar_version = models.ForeignKey(ClinVarVersion, null=True, on_delete=PROTECT)
    human_protein_atlas_version = models.ForeignKey(HumanProteinAtlasAnnotationVersion, null=True, on_delete=PROTECT)
    ontology_version = models.ForeignKey(OntologyVersion, null=True, on_delete=PROTECT)

    @property
    def sub_annotations_inheritance_partitioning(self) -> list[str]:
        """ Old style Postgres partitions (inherits) """
        _sub_list = ['variant_annotation_version', 'clinvar_version', 'human_protein_atlas_version']
        if settings.ANNOTATION_GENE_ANNOTATION_VERSION_ENABLED:
            _sub_list.append('gene_annotation_version')
        return _sub_list

    @property
    def sub_annotations(self) -> list[str]:
        _sub_declarative_partitioning_list = ['ontology_version']
        return self.sub_annotations_inheritance_partitioning + _sub_declarative_partitioning_list

    def _get_inheritance_partitioning_sub_annotation_versions(self):
        sub_annotation_versions = []
        for field in self.sub_annotations_inheritance_partitioning:
            if annotation_version := getattr(self, field):
                sub_annotation_versions.append(annotation_version)
        return sub_annotation_versions

    def get_partition_names(self):
        partitions_names = {}
        for sub_version in self._get_inheritance_partitioning_sub_annotation_versions():
            for base_table_name in sub_version.RECORDS_BASE_TABLE_NAMES:
                partitions_names[base_table_name] = sub_version.get_partition_table(base_table_name=base_table_name)
        return partitions_names

    def sql_partition_transformer(self, sql):
        """ Modifies SQL generated by QuerySet
             @see library.django_utils.django_queryset_sql_transformer.get_queryset_with_transformer_hook  """

        for annotation_version in self._get_inheritance_partitioning_sub_annotation_versions():
            sql = annotation_version.sql_partition_transformer(sql)
        return sql

    @cached_property
    def is_valid(self) -> bool:
        try:
            self.validate()
            return True
        except InvalidAnnotationVersionError:
            return False

    def validate(self):
        missing_sub_annotations = []
        for field in self.sub_annotations:
            field_fk = f"{field}_id"  # Avoid fetching related data
            sub_annotation = getattr(self, field_fk)
            if sub_annotation is None:
                missing_sub_annotations.append(field)
        if missing_sub_annotations:
            missing = ", ".join([str(s) for s in missing_sub_annotations])
            raise InvalidAnnotationVersionError(f"AnnotationVersion: {self} missing sub annotations: {missing}")

        if self.gene_annotation_version:
            if vav_gar_id := self.variant_annotation_version.gene_annotation_release_id:
                gene_gar_id = self.gene_annotation_version.gene_annotation_release_id
                if vav_gar_id != gene_gar_id:
                    different_msg = f"Inconsistent GeneAnnotationRelease. Variant {vav_gar_id} vs Gene: {gene_gar_id}"
                    raise InvalidAnnotationVersionError(different_msg)

            ov_id = self.ontology_version_id
            gav_ov_id = self.gene_annotation_version.ontology_version_id
            if (ov_id and gav_ov_id) and (ov_id != gav_ov_id):
                msg = f"OntologyVersion {ov_id} != GeneAnnotationVersion OntologyVersion {gav_ov_id}"
                raise InvalidAnnotationVersionError(msg)

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

        def latest(klass, order_by: str):
            return klass.objects.order_by(order_by).last()

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
                "human_protein_atlas_version": latest(HumanProteinAtlasAnnotationVersion, 'annotation_date'),
                "ontology_version": latest(OntologyVersion, 'pk'),
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

    def __lt__(self, other):
        return self.pk < other.pk

    def __str__(self):
        sub_versions = [f"Variant: {self.variant_annotation_version}",
                        f"Gene: {self.gene_annotation_version}",
                        f"ClinVar: {self.clinvar_version}",
                        f"HPA: {self.human_protein_atlas_version}",
                        f"Ontology: {self.ontology_version}"]
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
    citation = models.ForeignKey(Citation, null=True, on_delete=CASCADE)

    class Meta:
        unique_together = ('gene_symbol', 'citation')


class GenePubMedCount(models.Model):
    """ From https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz """
    gene = models.ForeignKey(Gene, on_delete=CASCADE)
    count = models.IntegerField()
    cached_web_resource = models.ForeignKey(CachedWebResource, on_delete=CASCADE)

    @staticmethod
    def get_count_and_note_for_gene_symbol(gene_symbol_id) -> tuple[int, str]:
        gene_symbol = GeneSymbol.objects.get(pk=gene_symbol_id)
        genes_qs = gene_symbol.get_genes().filter(annotation_consortium=AnnotationConsortium.REFSEQ)
        notes = []
        count = 0
        if gene_counts := list(GenePubMedCount.objects.filter(gene__in=genes_qs)):
            for gc in gene_counts:
                count = max(count, gc.count)
                notes.append(f"{gc.gene}={gc.count}")
            modified = gene_counts[0].retrieved_date
        else:
            if first_gpmc := GenePubMedCount.objects.first():
                count = 0
                notes.append(f"No results (could be no RefSeq gene for {gene_symbol=})")
                modified = first_gpmc.retrieved_date
            else:
                raise ValueError("No GenePubMedCount entries")

        retrieved = modified.date().isoformat()
        notes.insert(0, f"NCBI curated gene2pubmed retrieved {retrieved}. Homo Sapien only. ")
        return count, " ".join(notes)

    @staticmethod
    def get_for_gene_symbol(gene_symbol_id):
        # TODO
        pass

    @property
    def retrieved_date(self):
        return self.cached_web_resource.created


class MutationalSignatureInfo(models.Model):
    signature_id = models.IntegerField(primary_key=True)
    cancer_types = models.TextField()
    proposed_aetiology = models.TextField()
    additional_mutational_features = models.TextField()
    comments = models.TextField()
