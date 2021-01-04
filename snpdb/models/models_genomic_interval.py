import logging

from django.contrib.auth.models import User
from django.db import models
from django.db.models import CASCADE
from django.urls.base import reverse
import os

from lazy import lazy

from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from library.genomics.bed_file import BedFileReader
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_enums import ImportStatus


class GenomicIntervalsCategory(models.Model):
    name = models.TextField()
    description = models.TextField()

    def __str__(self):
        return self.name


class GenomicIntervalsCollection(GuardianPermissionsAutoInitialSaveMixin, models.Model):
    """ Can be a BED file (processed_file will be set - won't have related GenomicIntervals)
        or standard collection (w/related GenomicIntervals) """
    name = models.TextField()
    category = models.ForeignKey(GenomicIntervalsCategory, on_delete=CASCADE)
    genome_build = models.ForeignKey(GenomeBuild, null=True, blank=True, on_delete=CASCADE)
    processed_file = models.TextField(null=True, blank=True)
    processed_records = models.IntegerField(null=True, blank=True)
    user = models.ForeignKey(User, on_delete=CASCADE)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)

    @lazy
    def num_intervals(self):
        if self.processed_file:
            return self.processed_records
        return self.genomicinterval_set.count()

    def genomic_interval_iterator(self):
        """ returns iterator of GenomicInterval (collection) or HTSeq.GenomicInterval (bed file)
            both classes share chrom/start/end fields """
        if self.processed_file is not None:
            with open(self.processed_file) as f:
                for feature in BedFileReader(f):
                    yield feature.iv
        else:
            for gi in self.genomicinterval_set.all():
                yield gi

    def delete(self, using=None, keep_parents=False):
        if self.processed_file and os.path.exists(self.processed_file):
            logging.debug("Deleting '%s'", self.processed_file)
            os.remove(self.processed_file)
        super().delete(using=using, keep_parents=keep_parents)

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('view_genomic_intervals', kwargs={"genomic_intervals_collection_id": self.pk})

    @classmethod
    def get_listing_url(cls):
        return reverse('data')


class GenomicInterval(models.Model):
    genomic_intervals_collection = models.ForeignKey(GenomicIntervalsCollection, null=True, on_delete=CASCADE)
    chrom = models.TextField()
    start = models.IntegerField()
    end = models.IntegerField()

    def clone(self):
        copy = self
        copy.pk = None
        copy.save()
        return copy

    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end}"
