"""
Annotation version diffs: counts of variants added / modified / removed / unchanged between two
VariantAnnotationVersions, with per-column value changes. Nothing in VariantGrid creates or shows these
any more - the models remain so historical rows can still be read by paper/scripts/paper_data_mining.py.
"""
from django.db import models
from django.db.models.deletion import CASCADE
from django.utils.timesince import timesince
from model_utils.managers import InheritanceManager

from annotation.models.models import VariantAnnotationVersion
from snpdb.models import VariantGridColumn


class VersionDiff(models.Model):
    objects = InheritanceManager()
    num_added = models.IntegerField()
    num_modified = models.IntegerField()
    num_removed = models.IntegerField()
    num_unchanged = models.IntegerField()

    @staticmethod
    def description(version_from, version_to):
        from_date = version_from.annotation_date.date()
        to_date = version_to.annotation_date.date()
        time_between = timesince(version_to.annotation_date,
                                 version_from.annotation_date)
        return f"{from_date} to {to_date} ({time_between})"


class VersionDiffFromToResult(models.Model):
    version_diff = models.ForeignKey(VersionDiff, on_delete=CASCADE)
    vg_column = models.ForeignKey(VariantGridColumn, on_delete=CASCADE)
    value_from = models.TextField(null=True)
    value_to = models.TextField(null=True)
    count = models.IntegerField()

    def __str__(self):
        name = f"{self.version_diff} {self.vg_column}"
        return f"{name}: {self.value_from} => {self.value_to}: {self.count}"


class VersionDiffChangeCountResult(models.Model):
    version_diff = models.ForeignKey(VersionDiff, on_delete=CASCADE)
    vg_column = models.ForeignKey(VariantGridColumn, on_delete=CASCADE)
    count = models.IntegerField()

    def __str__(self):
        return f"{self.version_diff} {self.vg_column}: changed {self.count}"


class VariantAnnotationVersionDiff(VersionDiff):
    version_from = models.ForeignKey(VariantAnnotationVersion, related_name='version_diff_from', on_delete=CASCADE)
    version_to = models.ForeignKey(VariantAnnotationVersion, related_name='version_diff_to', on_delete=CASCADE)

    def __str__(self):
        return self.description(self.version_from, self.version_to)
