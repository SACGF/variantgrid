from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.db import models
from django.db.models.deletion import CASCADE
from django.db.models.query_utils import Q

from snpdb.models.models_enums import ColumnAnnotationLevel


class VariantGridColumn(models.Model):
    """ Used to populate analysis VariantGrid with annotation.

        Normally, we take a Variant queryset and get values queryset using "variant_column" """
    grid_column_name = models.TextField(primary_key=True)
    variant_column = models.TextField()
    annotation_level = models.CharField(max_length=1, choices=ColumnAnnotationLevel.CHOICES, null=True)
    width = models.IntegerField(null=True)
    label = models.TextField()
    description = models.TextField()
    model_field = models.BooleanField(default=True)  # Standard field, can use Meta inspection to determine colmodel
    queryset_field = models.BooleanField(default=True)  # In queryset.values() (field or alias)

    def get_css_classes(self):
        css_classes = ["user-column"]
        annotation_dict = dict(ColumnAnnotationLevel.CHOICES)
        annotation_level_class = annotation_dict.get(self.annotation_level)
        if annotation_level_class:
            css_classes.append("%s-column" % annotation_level_class.lower())
        return " ".join(css_classes)

    def __str__(self):
        return self.grid_column_name


class CustomColumnsCollection(models.Model):
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

    @staticmethod
    def filter_for_user(user):
        public = Q(user__isnull=True)
        mine = Q(user=user)
        return CustomColumnsCollection.objects.filter(public | mine)

    @staticmethod
    def get_permission_check(user, custom_columns_collection_id, mode='r'):
        ccc = CustomColumnsCollection.objects.get(pk=custom_columns_collection_id)
        mine = ccc.user == user
        if mode == 'r':
            valid = mine or ccc.user is None
        else:
            valid = mine

        if not valid:
            msg = f"You don't have permission to view CustomColumnsCollection {ccc.pk} (mode='{mode}')"
            raise PermissionDenied(msg)
        return ccc

    def clone_for_user(self, user):
        name = f"{user}'s copy of {self.name}"
        clone_cc = CustomColumnsCollection(name=name,
                                           user=user)
        clone_cc.save()

        # Mandatory columns are already inserted
        for cc in self.customcolumn_set.exclude(column__grid_column_name__in=CustomColumnsCollection.MANDATORY_COLUMNS):
            cc.pk = None
            cc.custom_columns_collection = clone_cc
            cc.save()

        return clone_cc

    def can_write(self, user):
        return self.user == user

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
