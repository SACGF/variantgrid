import logging
from functools import cached_property
from typing import Union

import celery
from django.contrib.auth.models import User, Group
from django.contrib.postgres.fields.array import ArrayField
from django.core.exceptions import PermissionDenied
from django.db import models
from django.db.models import QuerySet
from django.db.models.aggregates import Max
from django.db.models.deletion import CASCADE, DO_NOTHING, PROTECT
from django.db.models.expressions import F, Value
from django.db.models.query_utils import Q, FilteredRelation
from django.db.models.signals import pre_delete, post_delete
from django.dispatch.dispatcher import receiver
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from guardian.shortcuts import get_objects_for_user

from library.django_utils import SortByPKMixin
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.django_utils.django_postgres import PostgresRealField
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from library.guardian_utils import DjangoPermission
from library.preview_request import PreviewModelMixin, PreviewKeyValue
from library.utils import invert_dict
from patients.models_enums import Zygosity
from snpdb.models.models_enums import ImportStatus, CohortGenotypeCollectionType
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import Variant
from snpdb.models.models_vcf import VCF, Sample


class Cohort(GuardianPermissionsAutoInitialSaveMixin, PreviewModelMixin, SortByPKMixin, TimeStampedModel):
    """ Cohort - a collection of samples

        We pack data from all of the samples (zygosity, allele_depth, read_depth, genotype_quality, phred_likelihood) into 1 row
        for quicker group level queries without joins to all of the sample table partitions.

        * We create cohort automatically for a VCF
        * You can also create a cohort from samples from different VCF files
            - This is done in SQL in snpdb.tasks.cohort_genotype_tasks.cohort_count_task

        Data is stored in CohortGenotype rows, which are partitioned by CohortGenotypeCollection """
    name = models.TextField()
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    version = models.IntegerField(null=False, default=0)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    vcf = models.OneToOneField(VCF, null=True, on_delete=CASCADE)  # Will be NULL for custom Cohorts from multi-vcfs
    # Deal with parent_cohort delete in snpdb.signals.signal_handlers.pre_delete_cohort
    parent_cohort = models.ForeignKey("self", null=True, related_name="sub_cohort_set", on_delete=DO_NOTHING)
    sample_count = models.IntegerField(null=True)

    def can_view(self, user_or_group: Union[User, Group]) -> bool:
        """ Also uses VCF permission """
        if self.vcf and self.vcf.can_view(user_or_group):
            return True
        return super().can_view(user_or_group)

    def can_write(self, user_or_group: Union[User, Group]) -> bool:
        """ Also uses VCF permission """
        if self.vcf and self.vcf.can_write(user_or_group):
            return True
        return super().can_write(user_or_group)

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-people-arrows"

    @classmethod
    def preview_if_url_visible(cls) -> str:
        return "patients"

    @property
    def preview(self) -> 'PreviewData':
        extras = []
        if samples := PreviewKeyValue.count(Sample, self.sample_count):
            extras.append(samples)
        return self.preview_with(
            identifier=self.name,
            summary_extra=extras
        )

    @property
    def has_genotype(self):
        if self.vcf:
            return self.vcf.has_genotype
        return True  # Created cohorts must contain genotype

    def increment_version(self):
        # Check if any samples not in parent cohort (can no longer be a sub cohort)
        if self.parent_cohort:
            my_samples = self.get_samples()
            parent_samples = self.parent_cohort.get_samples()
            if my_samples.exclude(pk__in=parent_samples).exists():
                self.parent_cohort = None  # No longer a sub cohort

        self.version += 1
        self.sample_count = self.cohortsample_set.count()
        self.save()

        if self.vcf is None:
            # For custom built Cohorts - need to ensure cohort_genotype_packed_field_index entries are correct
            if highest_packed_index := self._get_cohort_sample_field_max("cohort_genotype_packed_field_index"):
                # packed fields are unique per cohort so need to move them sufficiently high to be able to rearrange
                out_of_range = F("cohort_genotype_packed_field_index") + Value(highest_packed_index + 1)
                self.cohortsample_set.all().update(cohort_genotype_packed_field_index=out_of_range)

                # Cohort genotype packed fields are in sample order
                cohort_samples = []
                for i, cs in enumerate(self.cohortsample_set.order_by("sample")):
                    cs.cohort_genotype_packed_field_index = i
                    cohort_samples.append(cs)
                if cohort_samples:
                    # Do a bulk update so it doesn't trigger a save then call this func again
                    CohortSample.objects.bulk_update(cohort_samples, ["cohort_genotype_packed_field_index"])

        self.delete_old_counts()

    def delete_old_counts(self):
        qs = CohortGenotypeCollection.objects.filter(cohort=self, cohort_version__lt=self.version)
        num_to_delete = qs.update(marked_for_deletion=True)
        if num_to_delete:
            task = delete_old_cohort_genotypes_task.si()
            task.apply_async()

    def _get_cohort_sample_field_max(self, field):
        qs = self.cohortsample_set.all()
        data = qs.aggregate(highest_value=Max(field))
        return data["highest_value"]

    def _get_next_sort_order(self):
        existing_max = self._get_cohort_sample_field_max("sort_order")
        if existing_max is not None:
            return existing_max + 1
        return 0

    def add_sample(self, sample_id):
        # TODO: Check for existing
        try:
            ss = self.cohortsample_set.get(cohort=self, sample_id=sample_id)
        except CohortSample.DoesNotExist:
            i = self._get_next_sort_order()
            ss = CohortSample.objects.create(cohort=self,
                                             sample_id=sample_id,
                                             cohort_genotype_packed_field_index=i,
                                             sort_order=i)
            # Will call increment_version() to bump cohort
        return ss

    def get_cohort_samples(self):
        return self.cohortsample_set.all().select_related("sample", "sample__vcf").order_by("sort_order")

    def get_samples(self) -> QuerySet[Sample]:
        """ In sample.pk order (consistent regardless of cohort samples order)
            This is different than get_samples_qs so that we have the option of caching it
        """
        return self.get_samples_qs()

    def get_samples_qs(self) -> QuerySet[Sample]:
        return Sample.objects.filter(cohortsample__cohort=self).order_by("pk")

    def get_sample_ids(self):
        return self.get_samples_qs().values_list("pk", flat=True)

    def is_sub_cohort(self):
        return self.parent_cohort is not None

    def get_base_cohort(self):
        """ Underlying cohort (ie deals with sub cohorts) """
        if self.is_sub_cohort():
            # recursion in case we ever decide to allow sub-sub cohorts
            cohort = self.parent_cohort.get_base_cohort()
        else:
            cohort = self
        return cohort

    @cached_property
    def cohort_genotype_collection(self):
        """ This is used to get the alias - which is either UNCOMMON if doing a rare filter, or BOTH common/uncommon
            if not doing rare pop filter """
        cohort = self.get_base_cohort()
        return CohortGenotypeCollection.objects.get(cohort=cohort,
                                                    cohort_version=cohort.version,
                                                    collection_type=CohortGenotypeCollectionType.UNCOMMON)

    def get_vcf(self):
        return self.get_base_cohort().vcf

    def contains_all_samples(self, sample_ids):
        sample_set = set(sample_ids)
        cohort_sample_set = set(self.get_sample_ids())
        return not sample_set - cohort_sample_set

    def get_absolute_url(self):
        return reverse('view_cohort', kwargs={"cohort_id": self.pk})

    @classmethod
    def get_listing_url(cls):
        return reverse('cohorts')

    def create_sub_cohort(self, user, sample_list):
        sub_cohort_name = f"{self.name} sub cohort"
        sub_cohort = Cohort.objects.create(name=sub_cohort_name,
                                           user=user,
                                           parent_cohort=self,
                                           import_status=self.import_status,
                                           genome_build=self.genome_build)

        sample_index = {cs.sample: cs.cohort_genotype_packed_field_index for cs in self.cohortsample_set.all()}
        for i, sample in enumerate(sample_list):
            field_index = sample_index[sample]
            CohortSample.objects.create(cohort=sub_cohort, sample=sample,
                                        cohort_genotype_packed_field_index=field_index, sort_order=i)
        return sub_cohort

    @classmethod
    def get_for_user(cls, user, pk, write=False):
        cohort = get_object_or_404(Cohort, pk=pk)
        if write:
            dp = DjangoPermission.WRITE
        else:
            dp = DjangoPermission.READ

        if cohort.vcf:
            perm = DjangoPermission.perm(cohort.vcf, dp)
            if user.has_perm(perm, cohort.vcf):
                return cohort
        else:
            perm = DjangoPermission.perm(cohort, dp)
            if user.has_perm(perm, cohort):
                return cohort

        raise PermissionDenied(f"You do not have permissions to access cohort id {cohort.pk}")

    @classmethod
    def filter_for_user(cls, user, queryset=None, group_data=True, success_status_only=True, **kwargs):
        if queryset is None:
            cohort_qs = get_objects_for_user(user, 'snpdb.view_cohort', accept_global_perms=False)

            # Also ones we have access to vcfs
            vcfs_qs = VCF.filter_for_user(user, group_data=group_data)
            vcf_cohorts = Cohort.objects.filter(pk__in=vcfs_qs.values_list("cohort__pk", flat=True))
            queryset = cohort_qs | vcf_cohorts

        if success_status_only:
            queryset = queryset.filter(import_status=ImportStatus.SUCCESS)
        return queryset

    @staticmethod
    def get_cohort_containing_all_samples(sample_ids, extra_cohort_filter_kwargs=None,
                                          extra_q=None):
        """ Returns cohort if exists, or None """

        cohort_filter_kwargs = {'import_status': ImportStatus.SUCCESS,
                                'cohortsample__sample_id__in': sample_ids}
        if extra_cohort_filter_kwargs:
            cohort_filter_kwargs.update(extra_cohort_filter_kwargs)
        q = Q(**cohort_filter_kwargs)
        if extra_q:
            q &= extra_q

        for cohort in Cohort.objects.filter(q):
            if cohort.contains_all_samples(sample_ids):
                return cohort

        return None

    def __str__(self):
        return f"{self.name} ({self.sample_count} samples)"


@receiver(pre_delete, sender=Cohort)
def pre_delete_cohort(sender, instance, **kwargs):  # pylint: disable=unused-argument
    # Make sub cohorts not depend on this one (they'll have to be regenerated)
    for sub_cohort in instance.sub_cohort_set.all():
        # logging.info("Updating sub cohorts: %s", sub_cohort)
        sub_cohort.parent_cohort = None
        sub_cohort.import_status = ImportStatus.CREATED
        sub_cohort.save()


class CohortSample(models.Model):
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    cohort_genotype_packed_field_index = models.IntegerField()
    sort_order = models.IntegerField()

    class Meta:
        unique_together = ('cohort', 'cohort_genotype_packed_field_index')

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)
        self.cohort.increment_version()

    def delete(self, *args, **kwargs):
        super().delete(*args, **kwargs)
        self.cohort.increment_version()

    @staticmethod
    def get_for_user(user, cohort_sample_id):
        cohort_sample = get_object_or_404(CohortSample, pk=cohort_sample_id)
        Cohort.get_for_user(user, cohort_sample.cohort.pk)  # Check cohort for permissions
        return cohort_sample

    @staticmethod
    def filter_for_user(user):
        cohort_ids = Cohort.filter_for_user(user).values_list("pk", flat=True)
        return CohortSample.objects.filter(cohort__in=cohort_ids)

    @property
    def name(self):
        return self.sample.name

    @property
    def cohort_genotype_sql_index(self):
        return self.cohort_genotype_packed_field_index + 1

    def __str__(self):
        return self.name


class CohortGenotypeTaskVersion(TimeStampedModel):
    """ Used to track what version of snpdb.tasks.cohort_genotype_tasks.cohort_genotype_task was run """

    name = models.TextField(unique=True)

    def __str__(self):
        return self.name


class CohortGenotypeCommonFilterVersion(TimeStampedModel):
    """
        We split a VCF up into 2 partitions based on how common the variants are (also classification)

        @see https://github.com/SACGF/variantgrid/wiki/Cohort-Genotype-Common-Filters

        For utilities on this method, see "common_variants.py" """
    gnomad_version = models.TextField()
    gnomad_af_min = models.FloatField()
    # This value is from classification.enums.classification_enums.ClinicalSignificance
    # but don't want to bring dependency in from classification
    clinical_significance_max = models.CharField(max_length=1, null=True)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)

    def __str__(self):
        description = f"gnomAD: {self.gnomad_version} AF>{self.gnomad_af_min}"
        if self.clinical_significance_max:
            description += f" and classification <= {self.clinical_significance_max}"
        return description


class CommonVariantClassified(TimeStampedModel):
    """ We store this so we know we've handled a variant being classified """
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    common_filter = models.ForeignKey(CohortGenotypeCommonFilterVersion, on_delete=CASCADE)

    class Meta:
        unique_together = ('variant', 'common_filter')


class CohortGenotypeCollection(RelatedModelsPartitionModel):
    """ A collection of Multiple-genotype records for a set of variants, for fast multi-sample zygosity queries.
        Storing genotypes individually per sample requires lots of joins when doing trios/cohorts

        We pack together genotype info into a single row (CohortGenotype) similar to Gemini:
        @see https://gemini.readthedocs.io/en/latest/content/database_schema.html#genotype-information

        These records are created via:
            * VCF import (in which case cohort.vcf will be set)
            * cohort_genotype_task (joining samples from diff VCFs together into a cohort)

        A VCF can be split into 2 CGCs (2 partitions) based on common / uncommon
        see CohortGenotypeCommonFilterVersion above
    """

    RECORDS_BASE_TABLE_NAMES = ["snpdb_cohortgenotype"]
    RECORDS_FK_FIELD_TO_THIS_MODEL = "collection_id"
    PARTITION_LABEL_TEXT = "collection"
    DEFAULT_COMBINED_COUNT_COLUMN = "vc_zygosity_count"

    name = models.TextField(null=True)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    cohort_version = models.IntegerField()
    num_samples = models.IntegerField()  # Number of elements in CohortGenotype arrays (samples may be deleted)
    # task fields only populated when created via cohort_genotype_task
    celery_task = models.CharField(max_length=36, null=True)
    task_version = models.ForeignKey(CohortGenotypeTaskVersion, null=True, on_delete=CASCADE)
    marked_for_deletion = models.BooleanField(null=False, default=False)
    collection_type = models.CharField(max_length=1, choices=CohortGenotypeCollectionType.choices,
                                       default=CohortGenotypeCollectionType.UNCOMMON)
    # common_collection will be set on the 'uncommon' (interesting/rare) CGC
    common_collection = models.OneToOneField('self', null=True, related_name="uncommon", on_delete=CASCADE)
    # common filter will be set on the 'common' CGC
    common_filter = models.ForeignKey(CohortGenotypeCommonFilterVersion, null=True, on_delete=PROTECT)

    class Meta:
        unique_together = ('cohort', 'cohort_version', 'collection_type')

    def __str__(self) -> str:
        parts = [f"CohortGenotypeCollection: {self.cohort}"]
        if self.cohort_version != self.cohort.version:
            parts.append(f"(v.{self.cohort.version} != {self.cohort_version})")
        if self.common_filter:
            parts.append(str(self.common_filter))
        return " ".join(parts)

    def get_common_filter_info(self) -> str:
        default_or_rare = self.cohortgenotype_set.count()
        if cc := self.common_collection:
            common = self.common_collection.cohortgenotype_set.count()
            total = default_or_rare + common
            if total:
                common_percent = 100 * common / (default_or_rare + common)
                common_percent = f"{common_percent:.1f}%"
            else:
                common_percent = "n/a"

            common_filter_info = [
                str(cc.common_filter),
                f"default/rare: {default_or_rare:,}",
                f"common: {common:,} ({common_percent})",
            ]
        else:
            common_filter_info = [
                "No filter applied",
                f"total: {default_or_rare:,}",
            ]
        return " ".join(common_filter_info)

    @property
    def cohortgenotype_alias(self):
        return f"cohortgenotype_{self.pk}"

    def get_packed_column_alias(self, column) -> str:
        return f"{self.cohortgenotype_alias}_packed_{column}"

    def get_partition_table(self, base_table_name=None):
        if base_table_name == "common":
            partition_table = self.common_collection.get_partition_table()
        else:
            partition_table = super().get_partition_table()
        return partition_table

    def sql_partition_transformer(self, sql):
        # No need to do anything as we're joining to Variant
        return sql

    def get_annotation_kwargs(self, **kwargs) -> dict:
        """ For Variant.objects.annotate """
        annotation_kwargs = {}
        already_set = self.cohortgenotype_alias in kwargs.get("existing_annotation_kwargs", set())
        if not (already_set and kwargs.get("override") is False):
            collections = [self.pk]
            if kwargs.get("common_variants", True):
                collections.append(self.common_collection_id)
            cgc_condition = Q(cohortgenotype__collection__in=collections)
            annotation_kwargs[self.cohortgenotype_alias] = FilteredRelation('cohortgenotype', condition=cgc_condition)
        return annotation_kwargs

    def get_sample_zygosity_regex(self, sample_zygosities: dict, sample_require_zygosity: dict):
        """ sample_zygosities_dict = {sample_id : zygosities_set} """

        zygosity_sample_matches = ["."] * self.num_samples

        for cs in self.cohort.get_cohort_samples():
            sample_zyg = sample_zygosities.get(cs.sample)
            require_zygosity = sample_require_zygosity.get(cs.sample, True)
            regex_match = Zygosity.get_regex_match(sample_zyg, require_zygosity=require_zygosity)
            zygosity_sample_matches[cs.cohort_genotype_packed_field_index] = regex_match

        regex_string = ''.join(zygosity_sample_matches)
        return regex_string

    @cached_property
    def get_packed_index_by_sample_id(self):
        return dict(self.cohort.cohortsample_set.values_list("sample_id", "cohort_genotype_packed_field_index"))

    def get_sql_index_for_sample_id(self, sample_id):
        return self.get_packed_index_by_sample_id[sample_id] + 1  # 1 based

    def get_array_index_for_sample_id(self, sample_id):
        # For https://docs.djangoproject.com/en/2.1/ref/contrib/postgres/fields/#index-transforms
        return self.get_packed_index_by_sample_id[sample_id]  # 0 based

    def get_zygosity_q(self, sample_zygosities: dict, sample_require_zygosity: dict = None, exclude=False) -> Q:
        """ sample_zygosities = {sample : zygosities_set}
            sample_require_zygosity = {sample : True/False} - defaults to True
            exclude - invert query (not equals)
        """
        if all(not v for v in sample_zygosities.values()):
            # nothing selected
            q_none = Q(pk__isnull=True)
            return q_none

        if sample_require_zygosity is None:
            sample_require_zygosity = {}

        regex_string = self.get_sample_zygosity_regex(sample_zygosities, sample_require_zygosity)
        if exclude:
            # Inverting a query via ~Q() leads to extremely slow queries so inverting regex.
            # Use ! (negative lookahead) only works where regex is anchored at start of line
            regex_string = f"^((?!{regex_string}))"

        # If regex string is all "." (ie everything) then can optimise away
        non_wildcard = regex_string.replace(".", "")
        if not non_wildcard:
            # Show everything in cohort
            q = Q(**{f"{self.cohortgenotype_alias}__isnull": False})
        else:
            q = Q(**{f"{self.cohortgenotype_alias}__samples_zygosity__regex": regex_string})
        return q

    @staticmethod
    def annotate_all_counts(qs, column_name=DEFAULT_COMBINED_COUNT_COLUMN):
        """ returns annotated_qs, column_name """
        both = F("cohortgenotype__het_count") + F("cohortgenotype__hom_count")
        return qs.annotate(**{column_name: both}), column_name


@receiver(pre_delete, sender=CohortGenotypeCollection)
def cohort_genotype_collection_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    try:
        instance.delete_related_objects()
    except:
        pass

@receiver(post_delete, sender=CohortGenotypeCollection)
def cohort_genotype_collection_post_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    # Handled in post_delete as common_collection has on_delete=CASCADE which would delete instance...
    try:
        if instance.common_collection:
            instance.common_collection.delete()
    except CohortGenotypeCollection.DoesNotExist:
        pass


class SampleGenotype:
    """ CohortGenotype for a sample """
    def __init__(self, cohort_genotype: 'CohortGenotype', sample: Sample, sample_index: int):
        self._cohort_genotype = cohort_genotype
        self.variant = cohort_genotype.variant
        self.sample = sample
        self.zygosity = self._get_sample_value("samples_zygosity", sample_index)
        self.allele_depth = self._get_sample_value("samples_allele_depth", sample_index)
        allele_frequency = self._get_sample_value("samples_allele_frequency", sample_index)
        if sample.vcf.allele_frequency_percent:
            allele_frequency = VCF.convert_from_percent_to_unit(allele_frequency)
        self.allele_frequency = allele_frequency
        self.read_depth = self._get_sample_value("samples_read_depth", sample_index)
        self.genotype_quality = self._get_sample_value("samples_genotype_quality", sample_index)
        self.phred_likelihood = self._get_sample_value("samples_phred_likelihood", sample_index)

    def _get_sample_value(self, cg_field_name, sample_index):
        (_, empty_value) = CohortGenotype.COLUMN_IS_ARRAY_EMPTY_VALUE[cg_field_name]
        array = getattr(self._cohort_genotype, cg_field_name)
        if array is None:
            value = empty_value
        else:
            value = array[sample_index]
        return value

    def get_vcf_filters(self) -> list:
        """ List of VCF Filters (string) """
        vcf_filters = ["."]  # May be missing (if not VCFLocusFilter)
        filter_lookup = invert_dict(self.sample.vcf.get_filter_dict())
        if self._cohort_genotype.filters:
            vcf_filters = [filter_lookup[f] for f in set(self._cohort_genotype.filters)]
        else:
            vcf_filters = ["PASS"]
        return vcf_filters

    def __str__(self):
        return f"{self.variant} ({self.sample})"


class CohortGenotype(models.Model):
    """ Genotype information for multiple samples in a single database row. """
    MISSING_NUMBER_VALUE = -1
    MISSING_FT_VALUE = "NULL"  # This can't be '.' as that could be a code....

    COLUMN_IS_ARRAY_EMPTY_VALUE = {
        "samples_zygosity": (False, "."),
        "samples_allele_depth": (True, MISSING_NUMBER_VALUE),
        "samples_allele_frequency": (True, MISSING_NUMBER_VALUE),
        "samples_read_depth": (True, MISSING_NUMBER_VALUE),
        "samples_genotype_quality": (True, MISSING_NUMBER_VALUE),
        "samples_phred_likelihood": (True, MISSING_NUMBER_VALUE),
        "samples_filters": (True, MISSING_FT_VALUE),
    }

    collection = models.ForeignKey(CohortGenotypeCollection, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    ref_count = models.IntegerField(default=0)
    het_count = models.IntegerField(default=0)
    hom_count = models.IntegerField(default=0)  # hom_alt
    unk_count = models.IntegerField(default=0)  # Unknown (ie ./.)
    filters = models.TextField(null=True)
    # samples_ fields below packed in sorted order of CohortSample.cohort_genotype_packed_field_index
    samples_zygosity = models.TextField()  # 1 character per sample of Zygosity
    samples_allele_depth = ArrayField(models.IntegerField(), null=True)
    samples_allele_frequency = ArrayField(PostgresRealField(), null=True)
    samples_read_depth = ArrayField(models.IntegerField(), null=True)
    samples_genotype_quality = ArrayField(models.IntegerField(), null=True)
    samples_phred_likelihood = ArrayField(models.IntegerField(), null=True)
    # Same codes as filters, at sample level. None for no entry, PASS = ""
    samples_filters = ArrayField(models.TextField(), null=True)
    # This stores array of format fields (array of dicts - 1 per sample)
    format = models.JSONField(null=False, blank=True, default=list)
    # This stores any remaining INFO fields from VCF record not used in fields above
    info = models.JSONField(null=False, blank=True, default=dict)

    def get_sample_genotype(self, sample: Sample) -> SampleGenotype:
        sample_index = self.collection.get_array_index_for_sample_id(sample.pk)
        return SampleGenotype(self, sample, sample_index)

    def get_sample_genotypes(self) -> list[SampleGenotype]:
        sample_genotypes = []
        for i, sample in enumerate(self.collection.cohort.get_samples()):
            sample_genotypes.append(SampleGenotype(self, sample, i))
        return sample_genotypes

    class Meta:
        unique_together = ("collection", "variant")


class Trio(GuardianPermissionsAutoInitialSaveMixin, SortByPKMixin, TimeStampedModel):
    """ A simple pedigree used frequently for Mendellian disease (TrioNode in analysis)
        and karyomapping """
    name = models.TextField(blank=True)
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)
    cohort = models.ForeignKey(Cohort, on_delete=CASCADE)
    mother = models.ForeignKey(CohortSample, related_name='trio_mother', on_delete=CASCADE)
    mother_affected = models.BooleanField(default=False)
    father = models.ForeignKey(CohortSample, related_name='trio_father', on_delete=CASCADE)
    father_affected = models.BooleanField(default=False)
    proband = models.ForeignKey(CohortSample, related_name='trio_proband', on_delete=CASCADE)

    @classmethod
    def get_permission_class(cls):
        return Cohort

    def get_permission_object(self):
        # Trio permissions based on cohort
        return self.cohort

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(cohort__in=queryset)

    @property
    def genome_build(self):
        return self.cohort.genome_build

    def get_cohort_samples(self):
        return [self.mother, self.father, self.proband]

    def get_samples(self):
        return Sample.objects.filter(cohortsample__in=self.get_cohort_samples()).order_by("pk")

    def get_absolute_url(self):
        return reverse('view_trio', kwargs={"pk": self.pk})

    def get_listing_url(self):
        return reverse('trios')

    @property
    def mother_details(self):
        affected = "affected" if self.mother_affected else "unaffected"
        return f"{self.mother} ({affected})"

    @property
    def father_details(self):
        affected = "affected" if self.father_affected else "unaffected"
        return f"{self.father} ({affected})"

    def __str__(self):
        return self.name or f"Trio {self.pk}"


# This has to be in this file so we don't end up with circular references
@celery.shared_task(ignore_result=False)
def delete_old_cohort_genotypes_task():
    logging.debug("About to delete old CohortGenotypeCollection")

    for cc in CohortGenotypeCollection.objects.filter(marked_for_deletion=True):
        cc.delete()
    logging.debug("Finished deleting old cohort counts")
