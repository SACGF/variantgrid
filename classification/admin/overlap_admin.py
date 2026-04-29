from django.contrib.admin import TabularInline
from django.db.models import QuerySet

from classification.models import Overlap, OverlapContribution, OverlapContributionSkew
from classification.services.overlaps_services import OverlapServices
from snpdb.admin_utils import ModelAdminBasics, admin_action, admin_list_column
from django.contrib import admin

from snpdb.models import AlleleOrigin, Allele


# class OverlapContributionInline(admin.TabularInline):
#     model = OverlapContribution
#     fields = ['source', 'classification_grouping', 'value', 'effective_date']
#
#     def has_add_permission(self, request, obj):
#         return False
#
#     def has_change_permission(self, request, obj=None):
#         return False
#
#     def has_delete_permission(self, request, obj=None):
#         return False

@admin.register(OverlapContribution)
class OverlapContributionAdmin(ModelAdminBasics):
    show_auditlog_history_link = True
    list_display = ['source', 'allele', 'classification_grouping', 'value_type', 'value', 'testing_context_bucket', 'effective_date__date', 'classification_grouping__lab']
    list_filter = ('source', 'value_type', 'testing_context_bucket', 'classification_grouping__lab')

    @admin_list_column(short_description="Effective Date", order_field="effective_date__date")
    def effective_date__date(self, obj: OverlapContribution):
        return obj.effective_date.date


@admin.register(OverlapContributionSkew)
class OverlapContributionSkewAdmin(ModelAdminBasics):
    list_display = ('overlap', 'contribution', 'skew_perspective', 'contribution__classification_grouping__lab')
    list_filter = ('contribution__classification_grouping__lab', 'overlap__overlap_status', 'skew_perspective')


@admin.register(Overlap)
class OverlapAdmin(ModelAdminBasics):
    list_display = ('overlap_status', 'valid', 'overlap_type', 'value_type', 'allele', 'testing_context_bucket', 'tumor_type_category', 'contributions_list', 'modified_detailed')
    # inlines = (OverlapContributionInline, )
    search_fields = ('pk', 'allele__id')
    list_filter = ('overlap_status', 'valid', 'overlap_type', 'value_type')

    @admin_list_column()
    def contributions_list(self, obj: Overlap):
        return ", ".join(str(contribution) for contribution in obj.contributions.all())

    @admin_list_column(short_description="Modified", order_field="modified")
    def modified_detailed(self, obj: Overlap):
        return self.format_datetime(obj.modified)

    # @admin_list_column(short_description="overlap_status_display", order_field="overlap_status")
    # def overlap_status_display(self, obj: Overlap):
    #     return obj.get_overlap_status_display()

    @admin_action("Refresh Overlap")
    def refresh_overlap(self, request, queryset: QuerySet[Overlap]):
        for overlap in queryset:
            OverlapServices.recalc_overlap(overlap)
            OverlapServices.update_skews(overlap)
