from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.db.models import Q
from django.db.models.deletion import CASCADE
from guardian.shortcuts import get_objects_for_group
from lazy import lazy
from model_utils.models import TimeStampedModel

from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from library.guardian_utils import public_group
from snpdb.models.models_enums import ColumnAnnotationLevel, VCFInfoTypes


class VariantGridColumn(models.Model):
    """ Holds info about columns used in VariantGrid analyses

        "variant_column" is a django path (eg variantannotation__) passed to a Variant values queryset to return data """
    grid_column_name = models.TextField(primary_key=True)
    variant_column = models.TextField()
    annotation_level = models.CharField(max_length=1, choices=ColumnAnnotationLevel.choices, null=True)
    width = models.IntegerField(null=True)
    label = models.TextField()
    description = models.TextField()
    model_field = models.BooleanField(default=True)  # Standard field, can use Meta inspection to determine colmodel
    queryset_field = models.BooleanField(default=True)  # In queryset.values() (field or alias)

    def get_css_classes(self):
        css_classes = ["user-column"]
        if self.annotation_level:
            annotation_level_class = ColumnAnnotationLevel(self.annotation_level).label
            css_classes.append("%s-column" % annotation_level_class.lower())
        return " ".join(css_classes)

    @lazy
    def columns_version_description(self) -> str:
        q = Q(max_vep_columns_version__isnull=False) | Q(max_vep_columns_version__isnull=False)
        if cvf := self.columnvepfield_set.filter(q).first():
            return cvf.columns_version_description
        return ""

    def __str__(self):
        return self.grid_column_name


class ColumnVCFInfo(models.Model):
    """ Used to export columns to VCF (vcf_export_utils) """
    info_id = models.TextField(primary_key=True)
    column = models.OneToOneField(VariantGridColumn, on_delete=CASCADE)
    number = models.IntegerField(null=True, blank=True)
    type = models.CharField(max_length=1, choices=VCFInfoTypes.choices)
    description = models.TextField(null=True)

    def __str__(self):
        number = self.number or '.'
        return f"ID={self.info_id},number={number},type={self.type},descr: {self.description}"


class CustomColumnsCollection(GuardianPermissionsAutoInitialSaveMixin, TimeStampedModel):
    MANDATORY_COLUMNS = ["variant"]
    name = models.TextField()
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)  # null = Public
    version_id = models.IntegerField(null=False, default=0)

    def save(self, **kwargs):
        # Add mandatory columns to all custom columns
        initial_save = not self.pk
        super().save(**kwargs)
        if initial_save:
            for i, grid_column_name in enumerate(CustomColumnsCollection.MANDATORY_COLUMNS):
                column = VariantGridColumn.objects.get(grid_column_name=grid_column_name)
                CustomColumn.objects.create(custom_columns_collection=self, column=column, sort_order=i)

    @staticmethod
    def get_system_default():
        cc, _ = CustomColumnsCollection.objects.get_or_create(user=None, name=settings.DEFAULT_COLUMNS_NAME)
        return cc

    @staticmethod
    def get_system_default_id():
        return CustomColumnsCollection.get_system_default().pk

    def increment_version(self):
        self.version_id += 1
        self.save()

    def clone_for_user(self, user):
        name = f"{user}'s copy of {self.name}"
        clone_cc = CustomColumnsCollection(name=name, user=user)
        clone_cc.save()

        # Mandatory columns are already inserted
        for cc in self.customcolumn_set.exclude(column__grid_column_name__in=CustomColumnsCollection.MANDATORY_COLUMNS):
            cc.pk = None
            cc.custom_columns_collection = clone_cc
            cc.save()

        return clone_cc

    @classmethod
    def filter_public(cls):
        return get_objects_for_group(public_group(), cls.get_read_perm(), klass=cls)

    def __str__(self):
        who = self.user or 'global'
        return f"({who}): {self.name}"


class CustomColumn(models.Model):
    custom_columns_collection = models.ForeignKey(CustomColumnsCollection, on_delete=CASCADE)
    column = models.ForeignKey(VariantGridColumn, on_delete=CASCADE)
    sort_order = models.IntegerField()

    class Meta:
        unique_together = ("custom_columns_collection", "column")

    def save(self, **kwargs):
        self.custom_columns_collection.increment_version()
        super().save(**kwargs)

    def delete(self, **kwargs):
        self.custom_columns_collection.increment_version()
        super().delete(**kwargs)

    @property
    def variant_column(self):
        return self.column.variant_column

    def __str__(self):
        return self.column.__str__()
