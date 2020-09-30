from annotation import models
from django.contrib import admin


admin.site.register(models.AnnotationRun)
admin.site.register(models.AnnotationVersion)
admin.site.register(models.ColumnVCFInfo)
admin.site.register(models.ClinVar)
admin.site.register(models.VariantAnnotationVersion)
