# -*- coding: utf-8 -*-
from django.contrib import admin

from snpdb.admin import ModelAdminBasics
from variantclassification.models import EvidenceKey
from variantclassification.models.admin_forms import EvidenceKeyAdmin, VariantClassificationAdmin, \
    ClinicalContextAdmin
from variantclassification.models.clinical_context_models import ClinicalContext
from variantclassification.models.discordance_models import DiscordanceReportClassification,\
    DiscordanceReport
from variantclassification.models.variant_classification import VariantClassification


# Register your models here.
admin.site.register(EvidenceKey, EvidenceKeyAdmin)
admin.site.register(VariantClassification, VariantClassificationAdmin)
admin.site.register(ClinicalContext, ClinicalContextAdmin)

admin.site.register(DiscordanceReport, ModelAdminBasics)
admin.site.register(DiscordanceReportClassification, ModelAdminBasics)
