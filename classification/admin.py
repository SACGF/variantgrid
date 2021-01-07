# -*- coding: utf-8 -*-
from django.contrib import admin

from snpdb.admin import ModelAdminBasics
from classification.models import EvidenceKey, ConditionAlias, ConditionAliasAdmin, ClassificationReportTemplate, \
    ClassificationReportTemplateAdmin
from classification.models.admin_forms import EvidenceKeyAdmin, ClassificationAdmin, \
    ClinicalContextAdmin
from classification.models.clinical_context_models import ClinicalContext
from classification.models.discordance_models import DiscordanceReportClassification,\
    DiscordanceReport
from classification.models.classification import Classification


# Register your models here.

admin.site.register(EvidenceKey, EvidenceKeyAdmin)
admin.site.register(Classification, ClassificationAdmin)
admin.site.register(ClinicalContext, ClinicalContextAdmin)
admin.site.register(ConditionAlias, ConditionAliasAdmin)
admin.site.register(DiscordanceReport, ModelAdminBasics)
admin.site.register(DiscordanceReportClassification, ModelAdminBasics)
admin.site.register(ClassificationReportTemplate, ClassificationReportTemplateAdmin)