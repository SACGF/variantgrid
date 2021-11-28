from django.contrib import admin

from annotation import models

admin.site.register(models.AnnotationRun)
admin.site.register(models.AnnotationVersion)
admin.site.register(models.ClinVar)
admin.site.register(models.VariantAnnotationVersion)
