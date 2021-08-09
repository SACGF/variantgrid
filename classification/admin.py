# -*- coding: utf-8 -*-
from django.contrib import admin

from classification.models import ConditionText, ConditionTextMatch, ClinVarExport, ClinVarExportBatch, \
    ClinVarAllele
from classification.models import EvidenceKey, ClassificationReportTemplate
from classification.models.admin_forms import ClinicalContextAdmin, DiscordanceReportAdmin, ConditionTextMatchAdmin, \
    EvidenceKeyAdmin, ConditionTextAdmin, ClassificationReportTemplateAdmin, ClassificationAdmin, \
    ClinVarExportAdmin, ClinVarExportBatchAdmin, ClinVarAlleleAdmin
from classification.models.classification import Classification
from classification.models.clinical_context_models import ClinicalContext
from classification.models.discordance_models import DiscordanceReportClassification, DiscordanceReport
from snpdb.admin import ModelAdminBasics

# Register your models here.

admin.site.register(EvidenceKey, EvidenceKeyAdmin)
admin.site.register(Classification, ClassificationAdmin)
admin.site.register(ClinVarAllele, ClinVarAlleleAdmin)
admin.site.register(ClinVarExport, ClinVarExportAdmin)
admin.site.register(ClinVarExportBatch, ClinVarExportBatchAdmin)
admin.site.register(ConditionText, ConditionTextAdmin)
admin.site.register(ConditionTextMatch, ConditionTextMatchAdmin)
admin.site.register(ClinicalContext, ClinicalContextAdmin)
admin.site.register(DiscordanceReport, DiscordanceReportAdmin)
admin.site.register(ClassificationReportTemplate, ClassificationReportTemplateAdmin)
# FIXME add all the new ClinVarExport models