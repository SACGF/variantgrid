import logging
from functools import cached_property
from typing import Tuple

from django.conf import settings
from django.db import models
from django.db.models import F, FilteredRelation, Q, QuerySet
from django.db.models.deletion import CASCADE
from django_extensions.db.models import TimeStampedModel

from library.django_utils.django_partition import RelatedModelsPartitionModel
from snpdb.models import Variant
from snpdb.models.models_vcf import VCF, Sample


class VariantZygosityCountCollection(RelatedModelsPartitionModel):
    """ We store ref/het/hom counts for each variant - in a collection (own partition) """
    RECORDS_BASE_TABLE_NAMES = ["snpdb_variantzygositycount"]
    RECORDS_FK_FIELD_TO_THIS_MODEL = "collection_id"
    PARTITION_LABEL_TEXT = "collection"

    GLOBAL_ALIAS = "global_variant_zygosity"

    name = models.TextField(unique=True)
    description = models.TextField()

    @cached_property
    def alias(self) -> str:
        if self.name == settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION:
            return self.GLOBAL_ALIAS
        return f"zygosity_count_{self.pk}"

    @cached_property
    def ref_alias(self) -> str:
        return f"{self.alias}__ref_count"

    @cached_property
    def het_alias(self) -> str:
        return f"{self.alias}__het_count"

    @cached_property
    def hom_alias(self) -> str:
        return f"{self.alias}__hom_count"

    @cached_property
    def germline_counts_alias(self) -> str:
        return f"germline_counts_{self.pk}"

    @cached_property
    def all_zygosity_counts_alias(self) -> str:
        return f"all_zygosity_counts_{self.pk}"

    def get_annotation_kwargs(self, **kwargs):
        q_collection = Q(variantzygositycount__collection=self)
        return {self.alias: FilteredRelation('variantzygositycount', condition=q_collection),
                self.germline_counts_alias: F(self.het_alias) + F(self.hom_alias),
                self.all_zygosity_counts_alias: F(self.ref_alias) + F(self.het_alias) + F(self.hom_alias)}

    @staticmethod
    def get_global_germline_counts() -> 'VariantZygosityCountCollection':
        try:
            return VariantZygosityCountCollection.objects.get(name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)
        except VariantZygosityCountCollection.DoesNotExist:
            logging.error("Could not find collection for settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION = %s",
                          settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)
            raise

    @staticmethod
    def annotate_global_germline_counts(qs: QuerySet) -> Tuple[QuerySet, 'VariantZygosityCountCollection']:
        vzcc = VariantZygosityCountCollection.get_global_germline_counts()
        return vzcc.annotate_all_germline_counts(qs), vzcc

    def annotate_all_germline_counts(self, qs: QuerySet) -> QuerySet:
        """ returns annotated_qs, column_name """
        kwargs = self.get_annotation_kwargs()
        return qs.annotate(**kwargs)

    def __str__(self) -> str:
        return self.name


class VariantZygosityCount(models.Model):
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    collection = models.ForeignKey(VariantZygosityCountCollection, on_delete=CASCADE)
    ref_count = models.IntegerField(default=0)  # HOM_REF - some reads, no zyg call (usually somatic)
    het_count = models.IntegerField(default=0)
    hom_count = models.IntegerField(default=0)  # hom_alt
    unk_count = models.IntegerField(default=0)  # Unknown (ie ./.)

    class Meta:
        unique_together = ("variant", "collection")


class AbstractVariantZygosityCountRecord(TimeStampedModel):
    collection = models.ForeignKey(VariantZygosityCountCollection, on_delete=CASCADE)
    count_complete = models.DateTimeField(null=True)
    deleted = models.DateTimeField(null=True)

    def check_can_delete(self):
        count_object = self.get_count_object()
        if self.count_complete is None:
            msg = f"GlobalZygosityCount for {count_object}: count_complete was None (never completed) - aborting"
            raise ValueError(msg)

        if self.deleted is not None:
            msg = f"GlobalZygosityCount for {count_object}: Attempting to delete record already deleted on {self.deleted} - aborting"
            raise ValueError(msg)

    def get_count_object(self):
        raise NotImplementedError()

    class Meta:
        abstract = True


class VariantZygosityCountForVCF(AbstractVariantZygosityCountRecord):
    collection = models.ForeignKey(VariantZygosityCountCollection, on_delete=CASCADE)
    vcf = models.ForeignKey(VCF, on_delete=CASCADE)
    is_split_to_sample_counts = models.BooleanField(default=False)

    class Meta:
        unique_together = ("collection", "vcf")

    def get_count_object(self):
        return self.vcf

    def split_to_sample_counts(self):
        """ May need to do this if you delete a single from a sample which was counted via CohortGenotype """
        msg = f"Splitting VariantZygosityCountForVCF ({self.vcf}) into per-sample counts"
        logging.info(msg)

        for sample in self.vcf.sample_set.all():
            VariantZygosityCountForSample.objects.create(sample=sample,
                                                         count_complete=self.count_complete,
                                                         deleted=self.deleted)

        self.is_split_to_sample_counts = True
        self.save()


# This is used to keep track of whether GVZC has been updated for a sample
class VariantZygosityCountForSample(AbstractVariantZygosityCountRecord):
    collection = models.ForeignKey(VariantZygosityCountCollection, on_delete=CASCADE)
    sample = models.ForeignKey(Sample, on_delete=CASCADE)

    class Meta:
        unique_together = ("collection", "sample")

    def get_count_object(self):
        return self.sample
