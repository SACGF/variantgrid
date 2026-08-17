from django.contrib import admin

from mme.models import MMEInboundMatch, MMEInboundQuery, MMEMatchResult, MMESubmission
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
    list_display = ("pk", "peer_node_id", "num_results", "created")
    list_filter = ("peer_node_id",)


@admin.register(MMEInboundMatch)
class MMEInboundMatchAdmin(ModelAdminBasics):
    list_display = ("pk", "inbound_query", "classification", "score", "remote_patient_id",
                    "notified", "created")
    search_fields = ("remote_patient_id",)
