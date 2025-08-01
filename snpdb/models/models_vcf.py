import logging
import operator
from collections import namedtuple
from functools import cached_property, reduce
from typing import Optional, Union

from django.conf import settings
from django.contrib.auth.models import User, Group
from django.contrib.postgres.fields.array import ArrayField
from django.core.exceptions import PermissionDenied, ObjectDoesNotExist
from django.db import models
from django.db.models import Lookup, Field
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.functions import Substr
from django.db.models.query_utils import Q
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from guardian.shortcuts import get_objects_for_user

from library.django_utils import SortByPKMixin
from library.genomics.vcf_enums import VariantClass
from library.guardian_utils import DjangoPermission
from library.log_utils import log_traceback, report_event
from library.preview_request import PreviewModelMixin, PreviewKeyValue
from patients.models import FakeData, Patient, Specimen
from patients.models_enums import Sex
from snpdb.models.models import Tag, LabProject
from snpdb.models.models_enums import ImportStatus, VariantsType, ProcessingStatus, SampleFileType, VCFInfoTypes
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_genomic_interval import GenomicIntervalsCollection
from snpdb.models.models_variant import Variant, VariantCollection, AlleleSource


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


class VCF(models.Model, PreviewModelMixin):
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
    # If AF is provided, we use it, otherwise if it is null we calculate it ourselves (post normalization w/VT)
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

    class Meta:
        verbose_name = 'VCF'
        verbose_name_plural = 'VCFs'

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
    def filter_for_user(user, group_data=True, has_write_permission=False):
        if has_write_permission:
            perm = DjangoPermission.perm(VCF, DjangoPermission.WRITE)
        else:
            perm = DjangoPermission.perm(VCF, DjangoPermission.READ)

        if group_data:
            queryset = get_objects_for_user(user, perm, klass=VCF, accept_global_perms=True)
        else:
            queryset = VCF.objects.filter(user=user)

        return queryset.exclude(import_status__in=ImportStatus.DELETION_STATES)

    @staticmethod
    def get_for_user(user, vcf_id):
        vcf = get_object_or_404(VCF, pk=vcf_id)
        read_perm = DjangoPermission.perm(VCF, DjangoPermission.READ)
        if not user.has_perm(read_perm, vcf):
            raise PermissionDenied()
        return vcf

    def can_view(self, user_or_group: Union[User, Group]) -> bool:
        read_perm = DjangoPermission.perm(VCF, DjangoPermission.READ)
        return user_or_group.has_perm(read_perm, self)

    def can_write(self, user_or_group: Union[User, Group]) -> bool:
        write_perm = DjangoPermission.perm(VCF, DjangoPermission.WRITE)
        return user_or_group.has_perm(write_perm, self)

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

    def delete_internal_data(self):
        """ Remove internal data but keep VCF and samples for reloading in place """

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

        # Delete and recreate rather than truncate as we may have had a schema change since
        logging.warning("*** Deleting then recreating partitions ***")
        for p in partitions:
            p.delete_related_objects()
            p.create_partition()


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
    def get_formatter(vcf: VCF):
        lookup = {vf.filter_code: vf.filter_id for vf in vcf.vcffilter_set.all()}

        def filter_string_formatter(row, field):
            if filter_string := row[field]:
                formatted_filters = []
                for f in filter_string:
                    formatted_filters.append(lookup[f])
                formatted_filters = ','.join(formatted_filters)
            else:
                formatted_filters = "PASS"
            return formatted_filters

        return filter_string_formatter


class VCFTag(models.Model):
    # Is this being used?
    tag = models.ForeignKey(Tag, on_delete=CASCADE)
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)

    def __str__(self):
        return f"{self.tag}:{self.vcf}"


class Sample(SortByPKMixin, PreviewModelMixin, models.Model):
    """ A VCF sample storing genotype information
        Sample data is stored as packed fields in CohortGenotype (via vcf.cohort.cohortgenotypecollection) """
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)
    vcf_sample_name = models.TextField()  # Original - can't be modified
    name = models.TextField()  # Initially set from vcf_sample_name
    no_dna_control = models.BooleanField(default=False)
    patient = models.ForeignKey(Patient, null=True, blank=True, on_delete=SET_NULL)
    # TODO: A sample may have >1 specimens (eg tumor/normal subtraction)
    specimen = models.ForeignKey(Specimen, null=True, blank=True, on_delete=SET_NULL)
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
    def genome_build(self):
        return self.vcf.genome_build

    @property
    def has_genotype(self):
        return self.vcf.has_genotype

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

    def delete_internal_data(self):
        """ for reloading in place """
        try:
            self.samplestats.delete()
        except:
            pass

        try:
            self.samplestatspassingfilter.delete()
        except:
            pass

        related_objects = [
            self.sampleclinvarannotationstats_set,
            self.sampleclinvarannotationstatspassingfilter_set,
            self.samplegeneannotationstats_set,
            self.samplegeneannotationstatspassingfilter_set,
            self.samplelocuscount_set,
            self.samplevariantannotationstats_set,
            self.samplevariantannotationstatspassingfilter_set,
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
        """ For annotating Variant queries """
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
    def filter_for_user(user, group_data=True, has_write_permission=False):
        """ May be given permission for whole VCF or just a sample """
        vcfs = VCF.filter_for_user(user, group_data, has_write_permission=has_write_permission)
        q_filters = [Q(vcf__in=vcfs)]
        if group_data:
            if has_write_permission:
                perm = DjangoPermission.perm(Sample, DjangoPermission.WRITE)
            else:
                perm = DjangoPermission.perm(Sample, DjangoPermission.READ)
            queryset = get_objects_for_user(user, perm, klass=Sample, accept_global_perms=True)
            q_filters.append(Q(pk__in=queryset.values_list("pk", flat=True)))

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
        except:
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
            params["patient_code"] = patient.last_name
            params["patient"] = str(patient)

        if specimen := self.specimen:
            params["specimen_id"] = specimen.pk
            params["specimen"] = str(specimen)
        return params

    @staticmethod
    def _validate_sample_formatter_func(sample_label_template):
        """ Throws error if invalid """
        specimen = Specimen(reference_id='refId', description='description')
        patient = Patient(pk=2, first_name='first_name', last_name='last_name')
        sample = Sample(pk=1, name="sample", patient=patient, specimen=specimen)
        params = sample._get_sample_formatter_params()
        errors = []
        for i, t in enumerate(sample_label_template.split("||")):
            try:
                t % params
            except (ValueError, KeyError) as exception:
                errors.append(f"{i+1}: '{t}: {exception=}'")
        if errors:
            error_msg = '\n'.join(errors)
            raise ValueError(f"Sample formatter function failed: {error_msg}")

    @staticmethod
    def _get_sample_formatter_func(sample_label_template, fallback=True):
        """ This is for rendering sample names on analysis grids """
        def _sample_formatter_func(sample):
            if sample_label_template:
                params = sample._get_sample_formatter_params()
                for t in sample_label_template.split("||"):
                    try:
                        return t % params
                    except (ValueError, KeyError) as e:
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


class SampleTag(models.Model):
    tag = models.ForeignKey(Tag, on_delete=CASCADE)
    sample = models.ForeignKey(Sample, on_delete=CASCADE)

    def __str__(self):
        return f"{self.tag}:{self.sample}"


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
        if self.vcf:
            qs = self.vcf.get_variant_qs()
        else:
            qs = Variant.objects.none()
        return qs

    def liftover_complete(self, genome_build: GenomeBuild):
        report_event('Completed VCF liftover',
                     extra_data={'vcf_id': self.vcf.pk, 'allele_count': self.get_allele_qs().count()})


class SampleStatsCodeVersion(TimeStampedModel):
    """ Track the version and code used to calculate sample stats, in case there are bugs/changes needed """
    name = models.TextField()
    version = models.IntegerField()
    code_git_hash = models.TextField()

    class Meta:
        unique_together = ("name", "version", "code_git_hash")

    def __str__(self):
        return f"{self.name} v{self.version}, git: {self.code_git_hash}, {self.created}"


class AbstractVariantStats(TimeStampedModel):
    """ Base class used for Cohort/Sample stats (note don't have Cohort stats yet)
        @see also annotation.models.models_sample_stats """

    code_version = models.ForeignKey(SampleStatsCodeVersion, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    variant_count = models.IntegerField(default=0)
    snp_count = models.IntegerField(default=0)
    insertions_count = models.IntegerField(default=0)
    deletions_count = models.IntegerField(default=0)
    ref_count = models.IntegerField(default=0)
    het_count = models.IntegerField(default=0)
    hom_count = models.IntegerField(default=0)
    unk_count = models.IntegerField(default=0)
    x_hom_count = models.IntegerField(default=0)
    x_het_count = models.IntegerField(default=0)
    x_unk_count = models.IntegerField(default=0)

    class Meta:
        abstract = True

    @staticmethod
    def percent(a, b):
        percent = float('NaN')
        if b:
            percent = 100.0 * a / b
        return percent

    @property
    def total_count(self):
        return self.variant_count

    @property
    def variant_percent(self):
        return AbstractVariantStats.percent(self.variant_count, self.total_count)

    @property
    def snp_percent(self):
        return AbstractVariantStats.percent(self.snp_count, self.total_count)

    @property
    def insertions_percent(self):
        return AbstractVariantStats.percent(self.insertions_count, self.total_count)

    @property
    def deletions_percent(self):
        return AbstractVariantStats.percent(self.deletions_count, self.total_count)

    def count_for_zygosity(self, zygosity_ref, zygosity_het, zygosity_hom, zygosity_unk, label=None):
        count = 0

        if zygosity_ref:
            count += self.ref_count
        if zygosity_het:
            count += self.het_count
        if zygosity_hom:
            count += self.hom_count
        if zygosity_unk:
            count += self.unk_count

        return count


class AbstractSampleStats(AbstractVariantStats):
    """ @see also annotation.models.models_sample_stats """
    sample = models.OneToOneField(Sample, on_delete=CASCADE)

    class Meta:
        abstract = True

    @classmethod
    def load_version(cls, sample, annotation_version):
        # Ignores annotation_version as not used, but want consistent interface
        return cls.objects.get(sample=sample)

    @property
    def chrx_sex_guess(self):
        """ returns sex by using hom/het ratio <0.2 = female, >0.8=male """

        sex = Sex.UNKNOWN
        if self.x_het_count and self.x_hom_count:
            hom_het_ratio = self.x_hom_count / self.x_het_count
            if hom_het_ratio < 0.2:
                sex = Sex.FEMALE
            elif hom_het_ratio > 0.8:
                sex = Sex.MALE

        return Sex(sex).label


class SampleStats(AbstractSampleStats):
    pass


class SampleStatsPassingFilter(AbstractSampleStats):
    pass


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
    source_regex = models.TextField()
    sample_variants_type = models.CharField(max_length=1, choices=VariantsType.choices, default=VariantsType.UNKNOWN)
    variant_zygosity_count = models.BooleanField(default=True)

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
        except:
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

        params = (self.name, self.get_status_display(), self.vcf, self.genomic_intervals, self.variant_collection, num_records)
        name = "%s (%s) proj: %s, genomic_intervals: %s variant_collection: %s (%d records)" % params
        if self.error_exception:
            name += self.error_exception

        return name
