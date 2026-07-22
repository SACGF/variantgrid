from django.contrib import admin

from mme.models import MMESubmission, MMEMatchResult, MMEInboundQuery
from snpdb.admin_utils import ModelAdminBasics


@admin.register(MMESubmission)
class MMESubmissionAdmin(ModelAdminBasics):
    list_display = ("pk", "classification", "node_id", "external_patient_id", "status",
                    "created", "submitted", "submitted_by")
    list_filter = ("status", "node_id")
    search_fields = ("external_patient_id",)


@admin.register(MMEMatchResult)
class MMEMatchResultAdmin(ModelAdminBasics):
    list_display = ("pk", "submission", "score", "matched_patient_id", "contact_name", "created")
    search_fields = ("matched_patient_id", "contact_name")


@admin.register(MMEInboundQuery)
class MMEInboundQueryAdmin(ModelAdminBasics):
    list_display = ("pk", "num_results", "created")
