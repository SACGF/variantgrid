from collections import namedtuple
from typing import Dict, List

from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied, ObjectDoesNotExist
from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django.db.models.functions import Substr
from django.db.models.query_utils import Q
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.shortcuts import get_object_or_404
from django.urls.base import reverse
from functools import reduce
from guardian.shortcuts import get_objects_for_user
from lazy import lazy
import logging
import operator

from library.django_utils import SortByPKMixin
from library.guardian_utils import DjangoPermission
from library.log_utils import log_traceback, report_event
from patients.models import FakeData, Patient, Specimen
from patients.models_enums import Sex
from snpdb.models.models import Tag, LabProject
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_genomic_interval import GenomicIntervalsCollection
from snpdb.models.models_variant import Variant, VariantCollection, AlleleSource
from snpdb.models.models_enums import ImportStatus, VariantsType, ProcessingStatus


class Project(models.Model):
    """ A way to group VCFs together """
    name = models.TextField(primary_key=True)
    description = models.TextField(null=True, blank=True)

    def __str__(self):
        name = self.name
        if self.description:
            name += f" ({self.description})"
        return name


class VCF(models.Model):
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
    # Most callers put allele depths in AD eg AD=[10,12] but some can split into separate ref/alt fields
    allele_depth_field = models.TextField(null=True)
    # If AF is provided, we use it, otherwise if it is null we calculate it ourselves (post normalization w/VT)
    # which can sometimes cause issues with splitting multi-alts
    allele_frequency_field = models.TextField(null=True)
    ref_depth_field = models.TextField(null=True)
    alt_depth_field = models.TextField(null=True)
    read_depth_field = models.TextField(null=True)
    genotype_quality_field = models.TextField(null=True)
    phred_likelihood_field = models.TextField(null=True)
    sample_filters_field = models.TextField(null=True)
    allele_frequency_percent = models.BooleanField(default=False)  # Legacy data used AF as percent
    # We don't want some VCFs to add to variant zygosity count (see VCFSourceSettings)
    variant_zygosity_count = models.BooleanField(default=True)

    @lazy
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

    def can_view(self, user):
        read_perm = DjangoPermission.perm(VCF, DjangoPermission.READ)
        return user.has_perm(read_perm, self)

    def can_write(self, user):
        write_perm = DjangoPermission.perm(VCF, DjangoPermission.WRITE)
        return user.has_perm(write_perm, self)

    def save(self, **kwargs):
        super().save(**kwargs)
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

    def get_warnings(self) -> List:
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

    @lazy
    def samples_by_vcf_name(self) -> Dict[str, 'Sample']:
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
    tag = models.ForeignKey(Tag, on_delete=CASCADE)
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)

    def __str__(self):
        return f"{self.tag}:{self.vcf}"


class Sample(SortByPKMixin, models.Model):
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
    bam_file_path = models.TextField(null=True, blank=True)
    variants_type = models.CharField(max_length=1, choices=VariantsType.choices, default=VariantsType.UNKNOWN)

    @property
    def genome_build(self):
        return self.vcf.genome_build

    @property
    def has_genotype(self):
        return self.vcf.has_genotype

    @property
    def is_somatic(self):
        return self.variants_type in VariantsType.SOMATIC_TYPES

    def can_view(self, user):
        vcf_read_perm = DjangoPermission.perm(self.vcf, DjangoPermission.READ)
        sample_read_perm = DjangoPermission.perm(self, DjangoPermission.READ)
        return user.has_perm(vcf_read_perm, self.vcf) or user.has_perm(sample_read_perm, self)

    def check_can_view(self, user):
        if not self.can_view(user):
            msg = f"You do not have permission to access sample {self.pk} (vcf {self.vcf.pk})"
            raise PermissionDenied(msg)

    def can_write(self, user):
        write_perm = DjangoPermission.perm(self, DjangoPermission.WRITE)
        return self.vcf.can_write(user) or user.has_perm(write_perm, self)

    def check_can_write(self, user):
        if not self.can_write(user):
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

    @lazy
    def cohort_genotype_collection(self):
        return self.vcf.cohort.cohort_genotype_collection

    def get_genotype(self, variant: Variant) -> 'SampleGenotype':
        try:
            cg = self.cohort_genotype_collection.cohortgenotype_set.get(variant=variant)
            sample_genotype = cg.get_sample_genotype(self)
        except ObjectDoesNotExist:
            sample_genotype = None
        return sample_genotype

    def __str__(self):
        return f"{self.name} ({self.vcf})"

    def get_absolute_url(self):
        return reverse('view_sample', kwargs={"sample_id": self.pk})

    @property
    def zygosity_alias(self):
        return f"sample_{self.pk}"

    def get_annotation_kwargs(self) -> Dict:
        """ For annotating Variant queries """
        cgc = self.cohort_genotype_collection
        i = cgc.get_sql_index_for_sample_id(self.pk)
        sample_zygosity = Substr(f"{cgc.cohortgenotype_alias}__samples_zygosity", i, length=1)
        return {self.zygosity_alias: sample_zygosity}

    def get_cohort_genotype_field(self, field_name):
        if not field_name.startswith("samples_"):
            field_name = f"samples_{field_name}"

        # samples_zygosity is a string not array, and should be added to w/get_annotation_kwargs
        if field_name == "samples_zygosity":
            field = self.zygosity_alias
        else:
            cgc = self.cohort_genotype_collection
            i = cgc.get_array_index_for_sample_id(self.pk)
            field = f"{cgc.cohortgenotype_alias}__{field_name}__{i}"
        return field

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

    @lazy
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


class AbstractVariantStats(models.Model):
    """ Base class used for Cohort/Sample stats
        @see also annotation.models.models_sample_stats """
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
            pbi = None
        return pbi

    @staticmethod
    def get_with_enrichment_kit_for_sample(sample):
        try:
            enrichment_kit = sample.samplefromsequencingsample.sequencing_sample.enrichment_kit
            pbi = VCFBedIntersection.get_for_vcf_and_enrichment_kit(sample.vcf, enrichment_kit)
            return pbi, enrichment_kit
        except:
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
