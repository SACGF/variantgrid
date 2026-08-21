import logging
import operator
import re
from collections import namedtuple
from functools import cached_property, reduce
from typing import Optional, Union

from django.conf import settings
from django.contrib.auth.models import Group, User
from django.contrib.postgres.fields.array import ArrayField
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied, ValidationError
from django.db import models
from django.db.models import Field, Lookup
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.functions import Substr
from django.db.models.query_utils import Q
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from guardian.shortcuts import get_objects_for_user

from library.django_utils import SortByPKMixin
from library.django_utils.data_archive_mixin import DataArchiveMixin
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.genomics.vcf_enums import VariantClass
from library.guardian_utils import DjangoPermission
from library.log_utils import log_traceback
from library.preview_request import PreviewKeyValue, PreviewModelMixin
from patients.models import ExtractionMatchMixin, FakeData, Patient, Specimen
from snpdb.models.models import LabProject
from snpdb.models.models_enums import (
    ImportStatus,
    ProcessingStatus,
    SampleFileType,
    VariantsType,
    VCFInfoTypes,
)
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_genomic_interval import GenomicIntervalsCollection
from snpdb.models.models_variant import AlleleSource, Variant, VariantCollection


@Field.register_lookup
class NotEqual(Lookup):
    """ From https://docs.djangoproject.com/en/4.0/howto/custom-lookups/#a-lookup-example """
    lookup_name = 'ne'

    def as_sql(self, compiler, connection):
        lhs, lhs_params = self.process_lhs(compiler, connection)
        rhs, rhs_params = self.process_rhs(compiler, connection)
        params = lhs_params + rhs_params
        return f'{lhs} <> {rhs}', params


class Project(models.Model):
    """ A way to group VCFs together """
    name = models.TextField(primary_key=True)
    description = models.TextField(null=True, blank=True)

    def __str__(self):
        name = self.name
        if self.description:
            name += f" ({self.description})"
        return name


class VCF(GuardianPermissionsMixin, DataArchiveMixin, PreviewModelMixin):
    name = models.TextField(null=True)
    date = models.DateTimeField()
    # genome_build will be set if imported successfully
    genome_build = models.ForeignKey(GenomeBuild, null=True, blank=True, on_delete=CASCADE)
    project = models.ForeignKey(Project, null=True, blank=True, on_delete=SET_NULL)
    user = models.ForeignKey(User, on_delete=CASCADE)
    genotype_samples = models.IntegerField()
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    fake_data = models.ForeignKey(FakeData, null=True, blank=True, on_delete=CASCADE)
    header = models.TextField(null=True)
    source = models.TextField(blank=True)
    # Most callers put allele depths in AD e.g. AD=[10,12] but some can split into separate ref/alt fields
    allele_depth_field = models.TextField(null=True)
    # If AF is provided, we use it, otherwise if it is null we calculate it ourselves (post normalization w/bcftools)
    # which can sometimes cause issues with splitting multi-alts
    allele_frequency_field = models.TextField(null=True)
    ref_depth_field = models.TextField(null=True)
    alt_depth_field = models.TextField(null=True)
    read_depth_field = models.TextField(null=True)
    genotype_field = models.TextField(null=True)
    genotype_quality_field = models.TextField(null=True)
    phred_likelihood_field = models.TextField(null=True)
    sample_filters_field = models.TextField(null=True)
    allele_frequency_percent = models.BooleanField(default=False)  # Legacy data used AF as percent
    # We don't want some VCFs to add to variant zygosity count (see VCFSourceSettings)
    variant_zygosity_count = models.BooleanField(default=True)
    # Set when an async archive (snpdb.tasks.vcf_archive_tasks) is queued/running. Cleared
    # in both terminal states: on success data_archived_date is stamped, on failure
    # data_archive_error is stamped (so a failed archive stays visible instead of the page
    # silently reverting to looking un-archived).
    data_archive_started_date = models.DateTimeField(null=True, blank=True)
    # Set if the async archive task raised; holds the error message. Cleared when a fresh
    # archive is queued (mark_vcf_archive_started) or once successfully archived.
    data_archive_error = models.TextField(null=True, blank=True)

    class Meta:
        verbose_name = 'VCF'
        verbose_name_plural = 'VCFs'

    @property
    def data_archive_in_progress(self) -> bool:
        """ True between queueing the archive task and it reaching a terminal state. """
        return self.data_archive_started_date is not None and not self.data_archived

    @property
    def data_archive_failed(self) -> bool:
        """ True when the last archive attempt errored and the VCF is still un-archived. """
        return bool(self.data_archive_error) and not self.data_archived

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-regular fa-file-lines"

    @classmethod
    def preview_if_url_visible(cls) -> bool:
        return 'data'

    @cached_property
    def has_filters(self):
        return self.vcffilter_set.exists()

    def get_filter_dict(self):
        filter_dict = dict(self.vcffilter_set.all().values_list("filter_id", "filter_code"))
        # Don't need this for CyVCF as it returns PASS as None but need for sample filters
        filter_dict["PASS"] = ""
        return filter_dict

    @staticmethod
    def convert_from_percent_to_unit(percent):
        from snpdb.models import CohortGenotype  # Circular import

        if percent != CohortGenotype.MISSING_NUMBER_VALUE:
            percent /= 100.0
        return percent

    @staticmethod
    def filter_for_user(user, group_data=True, has_write_permission=False, include_archived=True):
        if has_write_permission:
            perm = DjangoPermission.perm(VCF, DjangoPermission.WRITE)
        else:
            perm = DjangoPermission.perm(VCF, DjangoPermission.READ)

        if group_data:
            queryset = get_objects_for_user(user, perm, klass=VCF, accept_global_perms=True)
        else:
            queryset = VCF.objects.filter(user=user)

        queryset = queryset.exclude(import_status__in=ImportStatus.DELETION_STATES)
        if not include_archived:
            queryset = queryset.filter(data_archived_date__isnull=True)
        return queryset

    @staticmethod
    def get_for_user(user, vcf_id):
        vcf = get_object_or_404(VCF, pk=vcf_id)
        read_perm = DjangoPermission.perm(VCF, DjangoPermission.READ)
        if not user.has_perm(read_perm, vcf):
            raise PermissionDenied()
        return vcf

    def save(self, *args, **kwargs):
        super().save(*args, **kwargs)
        try:
            if self.cohort.name != self.name:
                self.cohort.name = self.name
                self.cohort.save()
        except ObjectDoesNotExist:
            pass

    def __str__(self):
        return f"{self.name}"

    def get_absolute_url(self):
        return reverse('view_vcf', kwargs={"vcf_id": self.pk})

    def get_warnings(self) -> list:
        VCFImportWarning = namedtuple('VCFImportWarning', ['message', 'has_more_details'])
        warnings = []
        if self.import_status == ImportStatus.REQUIRES_USER_INPUT and self.genome_build is None:
            warnings.append(VCFImportWarning("Unable to detect build, please select manually and save VCF.", False))
        return warnings

    @classmethod
    def get_listing_url(cls):
        return reverse('data')

    @property
    def has_genotype(self):
        return self.genotype_samples > 0

    @cached_property
    def samples_by_vcf_name(self) -> dict[str, 'Sample']:
        return {s.vcf_sample_name: s for s in self.sample_set.all()}

    def get_sample_ids(self):
        """ Often we just want this not all the objects """
        return self.sample_set.all().values_list("pk", flat=True)

    def get_variant_qs(self, qs=None):
        cgc = self.cohort.cohort_genotype_collection
        if qs is None:
            qs = Variant.objects.all()
        qs = qs.annotate(**cgc.get_annotation_kwargs())
        return qs.filter(**{f"{cgc.cohortgenotype_alias}__isnull": False})  # Inner join to CohortGenotype

    def delete_internal_data(self, recreate_partitions: bool = True):
        """ Remove internal data but keep VCF and samples for reloading in place.

            recreate_partitions=False is used by archive (snpdb/archive.py) — we want
            the partitions dropped, not replaced with empty ones.
        """

        # Remove VCF filters - some old ones had diff symbol/filter combos that cause errors trying to re-use
        self.vcffilter_set.all().delete()

        for sample in self.sample_set.all():
            sample.delete_internal_data()

        try:
            self.uploadedvcf.uploadedvcfpendingannotation.delete()
        except ObjectDoesNotExist:
            pass

        self.variantzygositycountforvcf_set.all().delete()

        partitions = []
        try:
            partitions += list(self.cohort.cohortgenotypecollection_set.all())
        except ObjectDoesNotExist:
            pass

        if recreate_partitions:
            # Delete and recreate rather than truncate as we may have had a schema change since
            logging.warning("*** Deleting then recreating partitions ***")
            for p in partitions:
                p.delete_related_objects()
                p.create_partition()
        else:
            logging.warning("*** Deleting partitions and CohortGenotypeCollection rows ***")
            for p in partitions:
                p.delete_related_objects()
                p.delete()


@receiver(pre_delete, sender=VCF)
def vcf_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    vcf = instance
    # logging.info("vcf_pre_delete_handler: %s", vcf.pk)
    try:
        cohort = vcf.cohort
        # logging.info("Deleting cohort for vcf %s", vcf.pk)
        cohort.delete()
    except ObjectDoesNotExist:
        pass


class AbstractVCFField(models.Model):
    """ Base class of INFO/FORMAT read from VCF Header """
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)
    identifier = models.TextField()  # ID
    number = models.TextField()  # not int so as to allow values like "." or "A"
    data_type = models.CharField(max_length=1, choices=VCFInfoTypes.choices)
    description = models.TextField()

    class Meta:
        abstract = True
        unique_together = ('vcf', 'identifier')


class VCFInfo(AbstractVCFField):
    COHORT_GENOTYPE_FIELD = "info"

class VCFFormat(AbstractVCFField):
    COHORT_GENOTYPE_FIELD = "format"


# These are VCF Filters (which are per-row - hence per-locus)
# We read them from the VCF header, and assign each filter to a character code.
# This is then stored (in alphabetical order) for each locus with filters.
class VCFFilter(models.Model):
    BASE_RECORDS_TABLE_NAME = "snpdb_vcflocusfilter"
    ASCII_MIN = 32
    ASCII_MAX = 126
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)
    filter_code = models.CharField(max_length=1)  # ASCII printable characters
    filter_id = models.TextField()
    description = models.TextField(null=True)

    class Meta:
        unique_together = (('vcf', 'filter_code'), ('vcf', 'filter_id'))

    @staticmethod
    def get_code_lookup(vcf: VCF) -> dict[str, str]:
        """ Filter codes are assigned per-VCF, so this can only decode that VCF's records """
        return {vf.filter_code: vf.filter_id for vf in vcf.vcffilter_set.all()}

    @staticmethod
    def format_filter_codes(lookup: dict[str, str], filter_string: Optional[str]) -> str:
        """ Empty/null means the record passed all of the VCF's filters """
        if filter_string:
            return ','.join(lookup.get(f, f) for f in filter_string)
        return "PASS"

    @staticmethod
    def get_formatter(vcf: VCF):
        lookup = VCFFilter.get_code_lookup(vcf)

        def filter_string_formatter(row, field):
            return VCFFilter.format_filter_codes(lookup, row[field])

        return filter_string_formatter


class Sample(GuardianPermissionsMixin, SortByPKMixin, PreviewModelMixin, ExtractionMatchMixin, models.Model):
    """ A VCF sample storing genotype information
        Sample data is stored as packed fields in CohortGenotype (via vcf.cohort.cohortgenotypecollection) """
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)
    vcf_sample_name = models.TextField()  # Original - can't be modified
    name = models.TextField()  # Initially set from vcf_sample_name
    no_dna_control = models.BooleanField(default=False)
    research_consent = models.BooleanField(null=True, blank=True)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    variants_type = models.CharField(max_length=1, choices=VariantsType.choices, default=VariantsType.UNKNOWN)

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-microscope"

    @classmethod
    def preview_if_url_visible(cls) -> Optional[str]:
        return 'patients'

    @property
    def preview(self) -> 'PreviewData':
        return self.preview_with(
            identifier=self.name,
            genome_builds={self.vcf.genome_build},
            summary_extra=[PreviewKeyValue("VCF", str(self.vcf))]
        )

    @property
    def specimen(self) -> Optional[Specimen]:
        if self.extraction:
            return self.extraction.specimen
        return None

    @property
    def genome_build(self):
        return self.vcf.genome_build

    @property
    def has_genotype(self):
        return self.vcf.has_genotype

    @property
    def data_archived(self) -> bool:
        return self.vcf.data_archived

    @property
    def is_somatic(self):
        return self.variants_type in VariantsType.SOMATIC_TYPES

    @cached_property
    def enrichment_kit(self):
        try:
            return self.samplefromsequencingsample.sequencing_sample.enrichment_kit
        except ObjectDoesNotExist:
            return None

    def get_minimum_coverage_required(self) -> int:
        if self.enrichment_kit:
            if self.enrichment_kit.min_coverage is not None:
                return self.enrichment_kit.min_coverage
        return settings.SEQAUTO_MIN_COVERAGE

    def can_view(self, user_or_group: Union[User, Group]) -> bool:
        read_perm = DjangoPermission.perm(self, DjangoPermission.READ)
        return self.vcf.can_view(user_or_group) or user_or_group.has_perm(read_perm, self)

    def check_can_view(self, user_or_group: Union[User, Group]):
        if not self.can_view(user_or_group):
            msg = f"You do not have permission to access sample {self.pk} (vcf {self.vcf.pk})"
            raise PermissionDenied(msg)

    def can_write(self, user_or_group: Union[User, Group]) -> bool:
        write_perm = DjangoPermission.perm(self, DjangoPermission.WRITE)
        return self.vcf.can_write(user_or_group) or user_or_group.has_perm(write_perm, self)

    def check_can_write(self, user_or_group: Union[User, Group]):
        if not self.can_write(user_or_group):
            msg = f"You do not have permission to modify sample {self.pk} (vcf {self.vcf.pk})"
            raise PermissionDenied(msg)

    @classmethod
    def allow_group_permission_delete(cls) -> bool:
        return True  # User data; deletable via the group_permissions delete view

    def delete_internal_data(self):
        """ for reloading in place """
        # Only the per-sample (sample IS NOT NULL) CohortGenotype*Stats rows belong to this Sample -
        # aggregate / filter-keyed rows are owned by the CGC and die when it's deleted
        related_objects = [
            self.samplelocuscount_set,
            self.cohortgenotypestats_set,
            self.cohortgenotypevariantannotationstats_set,
            self.cohortgenotypegeneannotationstats_set,
            self.cohortgenotypeclinvarannotationstats_set,
        ]

        for o in related_objects:
            o.all().delete()

    @cached_property
    def cohort_genotype_collection(self):
        return self.vcf.cohort.cohort_genotype_collection

    def get_genotype(self, variant: Variant) -> 'SampleGenotype':
        from snpdb.models import CohortGenotype

        cgc = self.cohort_genotype_collection
        collections = [cgc]
        if cgc.common_collection:
            collections.append(cgc.common_collection)
        try:
            cg = CohortGenotype.objects.filter(collection__in=collections).get(variant=variant)
            sample_genotype = cg.get_sample_genotype(self)
        except CohortGenotype.DoesNotExist:
            sample_genotype = None
        return sample_genotype

    def __str__(self):
        return f"{self.name} ({self.vcf})"

    def get_absolute_url(self):
        return reverse('view_sample', kwargs={"sample_id": self.pk})

    @property
    def zygosity_alias(self):
        return f"sample_{self.pk}"

    def get_annotation_kwargs(self, **kwargs) -> dict:
        """ For annotating Variant queries.

        SUBSTRING+IN beats regex for small numbers of constrained samples (single sample
        here, ~6.5x faster). The cohort path keeps regex because regex wins for wider
        cohorts — see CohortGenotypeCollection.get_zygosity_q and #1494.
        """
        cgc = self.cohort_genotype_collection
        i = cgc.get_sql_index_for_sample_id(self.pk)
        sample_zygosity = Substr(f"{cgc.cohortgenotype_alias}__samples_zygosity", i, length=1)
        return {self.zygosity_alias: sample_zygosity}

    def get_cohort_genotype_alias_and_field(self, field_name) -> tuple[str, str]:
        if not field_name.startswith("samples_"):
            field_name = f"samples_{field_name}"

        # samples_zygosity is a string not array, and should be added to w/get_annotation_kwargs
        if field_name == "samples_zygosity":
            alias = self.zygosity_alias
            field = self.zygosity_alias
        else:
            cgc = self.cohort_genotype_collection
            i = cgc.get_array_index_for_sample_id(self.pk)
            alias = cgc.cohortgenotype_alias
            field = f"{alias}__{field_name}__{i}"
        return alias, field

    def get_variant_qs(self, qs=None):
        """ Returns a Variant queryset inner joined to CohortGenotype with annotation aliases
            CohortGenotypeCollection.cohortgenotype_alias and Sample.zygosity_alias """
        if qs is None:
            qs = Variant.objects.all()
        annotation_kwargs = self.cohort_genotype_collection.get_annotation_kwargs()
        annotation_kwargs.update(self.get_annotation_kwargs())
        qs = qs.annotate(**annotation_kwargs)
        filter_kwargs = {
            f"{self.cohort_genotype_collection.cohortgenotype_alias}__isnull": False,  # Inner join to CohortGenotype
        }
        return qs.filter(**filter_kwargs)

    @classmethod
    def get_listing_url(cls):
        return reverse('data')

    @staticmethod
    def filter_for_user(user, group_data=True, has_write_permission=False, include_archived=True):
        """ May be given permission for whole VCF or just a sample """
        vcfs = VCF.filter_for_user(user, group_data, has_write_permission=has_write_permission,
                                   include_archived=include_archived)
        q_filters = [Q(vcf__in=vcfs)]
        if group_data:
            if has_write_permission:
                perm = DjangoPermission.perm(Sample, DjangoPermission.WRITE)
            else:
                perm = DjangoPermission.perm(Sample, DjangoPermission.READ)
            queryset = get_objects_for_user(user, perm, klass=Sample, accept_global_perms=True)
            sample_q = Q(pk__in=queryset.values_list("pk", flat=True))
            if not include_archived:
                sample_q &= Q(vcf__data_archived_date__isnull=True)
            q_filters.append(sample_q)

        ored_q_filters = reduce(operator.or_, q_filters)
        return Sample.objects.filter(ored_q_filters).exclude(import_status__in=ImportStatus.DELETION_STATES)

    @staticmethod
    def get_for_user(user, sample_id):
        sample = get_object_or_404(Sample.objects.select_related("vcf"), pk=sample_id)
        sample.check_can_view(user)
        return sample

    @cached_property
    def sequencing_run(self):
        try:
            return self.samplefromsequencingsample.sequencing_run
        except Exception:
            return None

    @staticmethod
    def soft_delete_samples_with_deleted_vcfs():
        vcfs_marked_for_deletion = Q(vcf__import_status=ImportStatus.MARKED_FOR_DELETION)
        not_already_deleting = ~Q(import_status=ImportStatus.DELETING)
        sample_mask = vcfs_marked_for_deletion & not_already_deleting
        return Sample.objects.filter(sample_mask).update(import_status=ImportStatus.MARKED_FOR_DELETION)

    def get_bam_files(self) -> list[str]:
        sfp_qs = SampleFilePath.objects.filter(sample=self, file_type=SampleFileType.BAM)
        return list(sfp_qs.values_list("file_path", flat=True))

    def _get_sample_formatter_params(self) -> dict[str, str]:
        params = {
            "sample_id": self.pk,
            "sample": self.name,
        }
        if patient := self.patient:
            params["patient_id"] = patient.pk
            params["patient_code"] = patient.patient_code or ""
            params["patient"] = str(patient)

        if specimen := self.specimen:
            params["specimen_id"] = specimen.pk
            params["specimen"] = str(specimen)
        return params

    @staticmethod
    def _get_sample_formatter_func(sample_label_template, fallback=True):
        """ This is for rendering sample names on analysis grids """
        def _sample_formatter_func(sample):
            if sample_label_template:
                params = sample._get_sample_formatter_params()
                for t in sample_label_template.split("||"):
                    try:
                        return t % params
                    except (ValueError, KeyError):
                        pass
            # In theory this should be valid due to form validator, but just in case
            return sample.name
        return _sample_formatter_func


class SampleFilePath(models.Model):
    """ Objects like a BAM/CRAM that can be associated with a sample """
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    file_type = models.CharField(max_length=1, null=True, blank=True, choices=SampleFileType.choices)
    label = models.TextField(null=True, blank=True)
    file_path = models.TextField()

    def __str__(self):
        if self.label:
            label = f" ({self.label}) "
        else:
            label = ''
        return f"{self.sample}: {label}{self.file_path}"


class VCFAlleleSource(AlleleSource):
    """ Used to link a VCF's variants to alleleles and liftover to other builds """
    vcf = models.ForeignKey(VCF, null=True, on_delete=SET_NULL)

    def get_genome_build(self):
        if self.vcf:
            genome_build = self.vcf.genome_build
        else:
            genome_build = None
        return genome_build

    def get_variant_qs(self):
        if self.vcf and not self.vcf.data_archived:
            qs = self.vcf.get_variant_qs()
        else:
            qs = Variant.objects.none()
        return qs


class SampleStatsCodeVersion(TimeStampedModel):
    """ Track the version and code used to calculate sample stats, in case there are bugs/changes needed """
    name = models.TextField()
    version = models.IntegerField()
    code_git_hash = models.TextField()

    class Meta:
        unique_together = ("name", "version", "code_git_hash")

    def __str__(self):
        return f"{self.name} v{self.version}, git: {self.code_git_hash}, {self.created}"


class VCFLengthStatsCollection(TimeStampedModel):
    vcf = models.OneToOneField(VCF, on_delete=CASCADE)
    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)


class VCFLengthStats(TimeStampedModel):
    collection = models.ForeignKey(VCFLengthStatsCollection, on_delete=CASCADE)
    variant_class = models.CharField(max_length=2, choices=VariantClass.choices, null=True)  # Null = unknown type
    is_log = models.BooleanField(default=False)
    histogram_counts = ArrayField(models.IntegerField())
    histogram_bin_edges = ArrayField(models.FloatField())

    class Meta:
        unique_together = ("collection", "variant_class")


class SampleLocusCount(models.Model):
    """ Count of variants at a locus for a sample """
    sample = models.ForeignKey(Sample, on_delete=CASCADE)
    locus_count = models.IntegerField()
    count = models.IntegerField()

    class Meta:
        unique_together = ("sample", "locus_count")


class SampleLabProject(models.Model):
    sample = models.OneToOneField(Sample, on_delete=CASCADE)
    lab_project = models.ForeignKey(LabProject, on_delete=CASCADE)


class VCFSourceSettings(models.Model):
    """ Modifies VCF based on 'source' header field - applied in upload.vcf.vcf_import.handle_vcf_source """

    # VCF sample-field columns a source may override. configure_vcf_from_header binds these by name, which
    # is wrong for callers that reuse a standard ID for something else (e.g. SpliceGirl's DP is *reference*
    # reads, not total depth) - so a source says what its fields actually mean
    OVERRIDABLE_SAMPLE_FIELDS = frozenset({
        "allele_depth_field",
        "allele_frequency_field",
        "ref_depth_field",
        "alt_depth_field",
        "read_depth_field",
        "genotype_field",
        "genotype_quality_field",
        "phred_likelihood_field",
        "sample_filters_field",
    })

    source_regex = models.TextField()
    sample_variants_type = models.CharField(max_length=1, choices=VariantsType.choices, default=VariantsType.UNKNOWN)
    variant_zygosity_count = models.BooleanField(default=True)
    # The build a source's files are called against - used when nothing else resolves one, ie the header
    # has no contigs/reference to detect from and the submitter declared none (@see resolve_genome_build)
    genome_build = models.ForeignKey(GenomeBuild, null=True, blank=True, on_delete=CASCADE)
    # A key present sets that field, including to null - which is how you clear a by-name default. A JSON
    # blob rather than nullable columns because null can't tell "no override" from "clear this field"
    sample_field_overrides = models.JSONField(default=dict, blank=True)

    @staticmethod
    def get_for_source(source: str) -> list['VCFSourceSettings']:
        """ Every setting whose regex matches the VCF's source, in the order they'll be applied """
        if not source:
            return []
        return [vss for vss in VCFSourceSettings.objects.all() if re.match(vss.source_regex, source)]

    def clean(self):
        super().clean()
        if unknown_fields := set(self.sample_field_overrides or {}) - self.OVERRIDABLE_SAMPLE_FIELDS:
            unknown = ", ".join(sorted(unknown_fields))
            valid = ", ".join(sorted(self.OVERRIDABLE_SAMPLE_FIELDS))
            raise ValidationError({"sample_field_overrides": f"Unknown VCF field(s): {unknown}. Valid: {valid}"})

    def apply_sample_field_overrides(self, vcf):
        """ Applied after the by-name defaults, so overrides land on top """
        for field, value in (self.sample_field_overrides or {}).items():
            if field not in self.OVERRIDABLE_SAMPLE_FIELDS:
                raise ValueError(f"{self}: '{field}' is not an overridable VCF sample field")
            setattr(vcf, field, value)

    def __str__(self):
        return f"{self.source_regex} sample_variants_type={self.get_sample_variants_type_display()}, " \
               f"variant_zygosity_count={self.variant_zygosity_count}"


class VCFBedIntersection(models.Model):
    name = models.TextField()
    status = models.CharField(max_length=1, choices=ProcessingStatus.choices, default=ProcessingStatus.CREATED)
    error_exception = models.TextField(null=True, blank=True)
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)
    genomic_intervals = models.ForeignKey(GenomicIntervalsCollection, on_delete=CASCADE)
    left_padding = models.IntegerField(default=0)
    right_padding = models.IntegerField(default=0)
    variant_collection = models.ForeignKey(VariantCollection, null=True, on_delete=CASCADE)

    @staticmethod
    def get_for_vcf_and_enrichment_kit(vcf, enrichment_kit):
        pbi = None
        try:
            kwargs = {'vcf': vcf,
                      'genomic_intervals': enrichment_kit.genomic_intervals,
                      'left_padding': settings.DEFAULT_ENRICHMENT_KIT_LEFT_PADDING,
                      'right_padding': settings.DEFAULT_ENRICHMENT_KIT_RIGHT_PADDING}

            pbi = VCFBedIntersection.objects.get(**kwargs)
        except VCFBedIntersection.DoesNotExist:
            pass
        except Exception:
            log_traceback()
        return pbi

    @staticmethod
    def get_with_enrichment_kit_for_sample(sample):
        if sample.enrichment_kit:
            pbi = VCFBedIntersection.get_for_vcf_and_enrichment_kit(sample.vcf, sample.enrichment_kit)
            return pbi, sample.enrichment_kit
        return None, None

    def __str__(self):
        num_records = 0
        if self.variant_collection:
            num_records = self.variant_collection.variantcollectionrecord_set.count()

        name = f"{self.name} ({self.get_status_display()}) proj: {self.vcf}, genomic_intervals: {self.genomic_intervals} variant_collection: {self.variant_collection} ({num_records} records)"
        if self.error_exception:
            name += self.error_exception

        return name
