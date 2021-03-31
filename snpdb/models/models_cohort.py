from typing import List, Dict

import celery
from django.contrib.postgres.fields.array import ArrayField
from django.core.exceptions import PermissionDenied
from django.db import models
from django.db.models.aggregates import Max
from django.db.models.deletion import CASCADE, DO_NOTHING
from django.db.models.expressions import F
from django.db.models.query_utils import Q, FilteredRelation
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from guardian.shortcuts import get_objects_for_user
from lazy import lazy
import logging

from library.django_utils import SortByPKMixin
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.django_utils.django_postgres import PostgresRealField
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.guardian_utils import DjangoPermission, assign_permission_to_user_and_groups
from library.utils import invert_dict
from patients.models_enums import Zygosity
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import Variant
from snpdb.models.models_vcf import VCF, Sample
from snpdb.models.models_enums import ImportStatus


class Cohort(SortByPKMixin, TimeStampedModel):
    """ Cohort - a collection of samples

        We pack data from all of the samples (zygosity, allele_depth, read_depth, genotype_quality, phred_likelihood) into 1 row
        for quicker group level queries without joins to all of the sample table partitions.

        * We create cohort automatically for a VCF
        * You can also create a cohort from samples from different VCF files
            - This is done in SQL in snpdb.tasks.cohort_genotype_tasks.cohort_count_task

        Data is stored in CohortGenotype rows, which are partitioned by CohortGenotypeCollection """
    name = models.TextField()
    version = models.IntegerField(null=False, default=0)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    vcf = models.OneToOneField(VCF, null=True, on_delete=CASCADE)
    # Deal with parent_cohort delete in snpdb.signals.signal_handlers.pre_delete_cohort
    parent_cohort = models.ForeignKey("self", null=True, related_name="sub_cohort_set", on_delete=DO_NOTHING)
    sample_count = models.IntegerField(null=True)

    @property
    def has_genotype(self):
        if self.vcf:
            return self.vcf.has_genotype
        return True  # Created cohorts must contain genotype

    def can_write(self, user):
        write_perm = DjangoPermission.perm(self, DjangoPermission.WRITE)
        return user.has_perm(write_perm, self)

    def increment_version(self):
        self.version += 1
        self.sample_count = self.cohortsample_set.all().count()
        self.save()

        self.delete_old_counts()

    def delete_old_counts(self):
        qs = CohortGenotypeCollection.objects.filter(cohort=self, cohort_version__lt=self.version)
        num_to_delete = qs.update(marked_for_deletion=True)
        if num_to_delete:
            task = delete_old_cohort_genotypes_task.si()
            task.apply_async()

    def get_next_sort_order(self):
        qs = self.cohortsample_set.all()
        data = qs.aggregate(Max("sort_order"))
        existing_max = data["sort_order__max"]
        if existing_max is not None:
            return existing_max + 1
        return 0

    def add_sample(self, sample_id):
        # TODO: Check for existing
        try:
            ss = self.cohortsample_set.get(cohort=self, sample_id=sample_id)
        except CohortSample.DoesNotExist:
            i = self.get_next_sort_order()
            ss = CohortSample.objects.create(cohort=self,
                                             sample_id=sample_id,
                                             cohort_genotype_packed_field_index=i,
                                             sort_order=i)
            # Will call increment_version() to bump cohort
        return ss

    def get_cohort_samples(self):
        return self.cohortsample_set.all().select_related("sample").order_by("sort_order")

    def get_samples(self):
        """ In sample.pk order (consistent regardless of cohort samples order) """
        return Sample.objects.filter(cohortsample__cohort=self).order_by("pk")

    def get_sample_ids(self):
        return self.get_samples().values_list("pk", flat=True)

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

    @lazy
    def cohort_genotype_collection(self):
        cohort = self.get_base_cohort()
        return CohortGenotypeCollection.objects.get(cohort=cohort, cohort_version=cohort.version)

    def get_vcf(self):
        return self.get_base_cohort().vcf

    def contains_all_samples(self, sample_ids):
        sample_set = set(sample_ids)
        cohort_sample_set = set(self.get_sample_ids())
        return not sample_set - cohort_sample_set

    def get_absolute_url(self):
        return reverse('view_cohort', kwargs={"cohort_id": self.pk})

    def get_sample_column_order_by(self, sample, column) -> str:
        """ Used to sort analysis grid """
        cgc = self.cohort_genotype_collection
        i = cgc.get_sql_index_for_sample_id(sample.pk)
        source = f"{cgc.get_partition_table()}.{column}"
        is_array = column != "samples_zygosity"
        if is_array:
            order_by = f"{source}[{i}]"
        else:
            order_by = f"substring({source}, {i}, 1)"
        return order_by

    @classmethod
    def get_listing_url(cls):
        return reverse('cohorts')

    def create_sub_cohort(self, user, sample_list):
        sub_cohort_name = f"{self.name} sub cohort"
        sub_cohort = Cohort.objects.create(name=sub_cohort_name,
                                           parent_cohort=self,
                                           import_status=self.import_status,
                                           genome_build=self.genome_build)
        assign_permission_to_user_and_groups(user, sub_cohort)

        sample_index = {cs.sample: cs.cohort_genotype_packed_field_index for cs in self.cohortsample_set.all()}
        for i, sample in enumerate(sample_list):
            field_index = sample_index[sample]
            CohortSample.objects.create(cohort=sub_cohort, sample=sample,
                                        cohort_genotype_packed_field_index=field_index, sort_order=i)
        return sub_cohort

    @staticmethod
    def get_for_user(user, cohort_id):
        cohort = get_object_or_404(Cohort, pk=cohort_id)
        if cohort.vcf:
            if user.has_perm('view_vcf', cohort.vcf):
                return cohort

        elif user.has_perm('view_cohort', cohort):
            return cohort

        raise PermissionDenied(f"You do not have permissions to access cohort id {cohort.pk}")

    @staticmethod
    def filter_for_user(user, group_data=True, success_status_only=True):
        cohort_qs = get_objects_for_user(user, 'snpdb.view_cohort', accept_global_perms=False)

        # Also ones we have access to vcfs
        vcfs_qs = VCF.filter_for_user(user, group_data=group_data)
        vcf_cohorts = Cohort.objects.filter(pk__in=vcfs_qs.values_list("cohort__pk", flat=True))
        cohort_qs |= vcf_cohorts
        if success_status_only:
            cohort_qs = cohort_qs.filter(import_status=ImportStatus.SUCCESS)
        return cohort_qs

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
def pre_delete_cohort(sender, instance, *args, **kwargs):
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

    def save(self, **kwargs):
        super().save(**kwargs)
        self.cohort.increment_version()

    def delete(self, **kwargs):
        super().delete(**kwargs)
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


class CohortGenotypeCollection(RelatedModelsPartitionModel):
    """ A collection of Multiple-genotype records for a set of variants, used to perform fast multi-sample zygosity queries.

        Storing genotypes as ObservedVariant records per sample requires lots of joins when doing multi-sample zygosity queries

        We pack together genotype info into a single row (CohortGenotype) similar to Gemini:
        @see https://gemini.readthedocs.io/en/latest/content/database_schema.html#genotype-information

        These records are created via:
            * VCF import (in which case cohort.vcf will be set)
            * cohort_genotype_task
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

    @property
    def cohortgenotype_alias(self):
        return f"cohortgenotype_{self.pk}"

    def get_annotation_kwargs(self) -> Dict:
        """ For Variant.objects.annotate """
        cgc_condition = Q(cohortgenotype__collection=self)
        return {self.cohortgenotype_alias: FilteredRelation('cohortgenotype', condition=cgc_condition)}

    def get_sample_zygosity_regex(self, sample_zygosities: dict, sample_require_zygosity: dict):
        """ sample_zygosities_dict = {sample_id : zygosities_set} """

        zygosity_sample_matches = []

        for sample in self.cohort.get_samples():
            sample_zyg = sample_zygosities.get(sample)
            require_zygosity = sample_require_zygosity.get(sample, True)
            regex_match = Zygosity.get_regex_match(sample_zyg, require_zygosity=require_zygosity)
            zygosity_sample_matches.append(regex_match)

        regex_string = ''.join(zygosity_sample_matches)
        return regex_string

    @lazy
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

        if sample_require_zygosity is None:
            sample_require_zygosity = {}

        regex_string = self.get_sample_zygosity_regex(sample_zygosities, sample_require_zygosity)
        if exclude:
            # Inverting a query via ~Q() leads to extremely slow queries so inverting regex.
            # Use ! (negative lookahead) only works where regex is anchored at start of line
            regex_string = "^((?!%s))" % regex_string

        # If regex string is all "." (ie everything) then can optimise away
        non_wildcard = regex_string.replace(".", "")
        if not non_wildcard:
            q = Q(pk__isnull=False)
        else:
            q = Q(**{f"{self.cohortgenotype_alias}__samples_zygosity__regex": regex_string})
        return q

    @staticmethod
    def annotate_all_counts(qs, column_name=DEFAULT_COMBINED_COUNT_COLUMN):
        """ returns annotated_qs, column_name """
        both = F("cohortgenotype__het_count") + F("cohortgenotype__hom_count")
        return qs.annotate(**{column_name: both}), column_name


@receiver(pre_delete, sender=CohortGenotypeCollection)
def cohort_genotype_collection_pre_delete_handler(sender, instance, **kwargs):
    try:
        instance.delete_related_objects()
    except:
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

    COLUMN_IS_ARRAY_EMPTY_VALUE = {
        "samples_zygosity": (False, "."),
        "samples_allele_depth": (True, -1),
        "samples_allele_frequency": (True, -1),
        "samples_read_depth": (True, -1),
        "samples_genotype_quality": (True, -1),
        "samples_phred_likelihood": (True, -1),
        "samples_filters": (True, "NULL"),
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

    def get_sample_genotype(self, sample: Sample) -> SampleGenotype:
        sample_index = self.collection.get_array_index_for_sample_id(sample.pk)
        return SampleGenotype(self, sample, sample_index)

    def get_sample_genotypes(self) -> List[SampleGenotype]:
        sample_genotypes = []
        for i, sample in enumerate(self.cohort.get_samples()):
            sample_genotypes.append(SampleGenotype(self, sample, i))
        return sample_genotypes


class Trio(GuardianPermissionsMixin, SortByPKMixin, models.Model):
    """ A simple pedigree used frequently for Mendellian disease (TrioNode in analysis)
        and karyomapping """
    name = models.TextField(blank=True)
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
@celery.task(ignore_result=False)
def delete_old_cohort_genotypes_task():
    logging.debug("About to delete old CohortGenotypeCollection")

    for cc in CohortGenotypeCollection.objects.filter(marked_for_deletion=True):
        cc.delete()
    logging.debug("Finished deleting old cohort counts")
