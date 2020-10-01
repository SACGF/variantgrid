# -*- coding: utf-8 -*-
from django.contrib import admin

from snpdb.admin import ModelAdminBasics
from classification.models import EvidenceKey
from classification.models.admin_forms import EvidenceKeyAdmin, VariantClassificationAdmin, \
    ClinicalContextAdmin
from classification.models.clinical_context_models import ClinicalContext
from classification.models.discordance_models import DiscordanceReportClassification,\
    DiscordanceReport
from classification.models.variant_classification import VariantClassification


# Register your models here.
admin.site.register(EvidenceKey, EvidenceKeyAdmin)
admin.site.register(VariantClassification, VariantClassificationAdmin)
admin.site.register(ClinicalContext, ClinicalContextAdmin)

admin.site.register(DiscordanceReport, ModelAdminBasics)
admin.site.register(DiscordanceReportClassification, ModelAdminBasics)
