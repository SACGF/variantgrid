"""
ClinVar: ClinVarVersion (one ClinVar VCF per build, partitioned) and its per-variant ClinVar rows, plus the
ClinVarRecordCollection / ClinVarRecord cache of individual submissions fetched from ClinVar's XML
(annotation/clinvar_fetch_request.py fills it) and ClinVar's citations.
"""

import os
import re
from datetime import datetime
from functools import cached_property
from typing import Optional

import django
from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models.deletion import CASCADE, DO_NOTHING, PROTECT, SET_NULL
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from model_utils.models import TimeStampedModel

from annotation.enums import ClinVarOncogenicity, ClinVarReviewStatus
from annotation.models.models_citations import Citation, CitationFetchRequest, CitationFetchResponse
from annotation.models.models_version import SubVersionPartition
from annotation.utils.clinvar_constants import CLINVAR_REVIEW_EXPERT_PANEL_STARS_VALUE
from library.django_utils.data_archive_mixin import DataArchiveMixin
from library.utils import first
from snpdb.enums import AlleleOriginBucket, SomaticClinicalSignificance
from snpdb.models import (
    Allele,
    GenomeBuild,
    Variant,
)
from snpdb.models.models_cached_web_resource import CachedWebResource


class ClinVarVersion(DataArchiveMixin, SubVersionPartition):
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
        CLINVAR_PATTERN_T2T = r"^Homo_sapiens-GCA_009914755.4-(20\d{2})_(\d+)-clinvar.vcf.gz"
        if m := re.match(CLINVAR_PATTERN_T2T, base_name):
            year = m.group(1)
            month = m.group(2)
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

    # DO_NOTHING as records are deleted by dropping the partition (a CASCADE locks every other version's)
    version = models.ForeignKey(ClinVarVersion, on_delete=DO_NOTHING)
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

    highest_oncogenicity = models.IntegerField(default=0)  # Highest of oncogenic_classification

    # ONCDN
    oncogenic_preferred_disease_name = models.TextField(null=True, blank=True)
    # ONCDISDB
    oncogenic_disease_database_name = models.TextField(null=True, blank=True)

    # SCIREVSTAT
    somatic_review_status = models.CharField(max_length=1, null=True, choices=ClinVarReviewStatus.choices)
    # SCI
    somatic_clinical_significance = models.TextField(null=True, blank=True)
    # Derived from SCI
    somatic_tier = models.CharField(max_length=20, null=True, blank=True,
                                    choices=SomaticClinicalSignificance.CHOICES)

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
                elif name.startswith("Human_Phenotype_Ontology:"):
                    name = name.removeprefix("Human_Phenotype_Ontology:")
                return name

            db_names = sorted(fix_name(db_name) for db_name in re.split("[|,]", db_name_text))
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
    def preferred_disease_names(self) -> list[str]:
        """ CLNDN packs every condition into one pipe joined run of underscored words - split it
            back out so it reads as a list and wraps """
        return ClinVar._pipe_names(self.preferred_disease_name)

    @staticmethod
    def _pipe_names(value: Optional[str]) -> list[str]:
        if not value:
            return []
        return [name.lstrip("_").replace("_", " ").strip() for name in value.split("|")]

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
    def somatic_tier_label(self) -> Optional[str]:
        return SomaticClinicalSignificance.LABELS.get(self.somatic_tier)

    def get_highest_oncogenicity_display(self) -> Optional[str]:
        if self.highest_oncogenicity:
            return ClinVarOncogenicity(self.highest_oncogenicity).label
        return None

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

    def __str__(self):
        return f"ClinVar: variant: {self.variant}, path: {self.highest_pathogenicity}"


clinvar_record_collection_refreshed = django.dispatch.Signal()
"""
Takers ClinVarRecordCollection as sender, and the instance of a ClinVarRecordCollection that has just been refreshed
as the 2nd parameter
"""


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

    @classmethod
    def set_allele_for_variants(cls, allele):
        """Proactively link ClinVarRecordCollections to allele for all its associated variants.
        Called when viewing the allele/variant page, so the allele is set before the ClinVar AJAX fires."""
        clinvar_variation_ids = list(
            ClinVar.objects.filter(variant__variantallele__allele=allele)
            .values_list('clinvar_variation_id', flat=True).distinct()
        )
        if clinvar_variation_ids:
            cls.objects.filter(
                clinvar_variation_id__in=clinvar_variation_ids,
                allele__isnull=True
            ).update(allele=allele)

    def records_with_min_stars(self, min_stars: int) -> list['ClinVarRecord']:
        return sorted(self.clinvarrecord_set.filter(stars__gte=min_stars), reverse=True)

    def update_with_records_and_save(self, records: list['ClinVarRecord']):
        records = sorted(records, reverse=True)
        self.clinvarrecord_set.exclude(record_id__in={r.record_id for r in records}).delete()
        for record in records:
            record.clinvar_record_collection = self

        # since the primary key of a ClinVarRecord is the record_id, we should just be able to upsert
        # If these newly created "ClinVarRecord"s have a record_id that matches a record already in the DB
        # it should just update anyway
        all_field_names = [field.name for field in ClinVarRecord._meta.concrete_fields if field.name != "record_id"]

        ClinVarRecord.objects.bulk_create(records, update_conflicts=True, unique_fields=["record_id"], update_fields=all_field_names)
        self.expert_panel = None
        self.max_stars = None
        if best_record := first(records):
            self.max_stars = best_record.stars
            if best_record.is_expert_panel_or_greater:
                self.expert_panel = best_record
        self.save()
        clinvar_record_collection_refreshed.send(ClinVarRecordCollection, instance=self)


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
        self._invalid = True

    def __bool__(self) -> bool:
        if hasattr(self, '_invalid'):
            return not self._invalid
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
    cached_web_resource = models.ForeignKey(CachedWebResource, null=True, on_delete=CASCADE)


class ClinVarCitation(models.Model):
    clinvar_citations_collection = models.ForeignKey(ClinVarCitationsCollection, on_delete=CASCADE)
    clinvar_variation_id = models.IntegerField()
    clinvar_allele_id = models.IntegerField()
    citation = models.ForeignKey(Citation, on_delete=CASCADE)
