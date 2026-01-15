from django.contrib.admin import TabularInline
from django.db.models import QuerySet

from classification.models import Overlap, OverlapContribution
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
    list_display = ['source', 'allele', 'classification_grouping', 'value_type', 'value', 'testing_context_bucket', 'effective_date']
    list_filter = ('source', 'value_type', 'testing_context_bucket')


@admin.register(Overlap)
class OverlapAdmin(ModelAdminBasics):
    list_display = ('overlap_status_display', 'valid', 'overlap_type', 'value_type', 'allele', 'testing_contexts', 'tumor_type_category', 'contributions_list', 'modified_detailed')
    # inlines = (OverlapContributionInline, )
    search_fields = ('pk', 'allele__id')
    list_filter = ('overlap_status', 'valid', )

    @admin_list_column()
    def contributions_list(self, obj: Overlap):
        return ", ".join(str(contribution) for contribution in obj.contributions.all())

    @admin_action("Re-calculate allele's overlaps")
    def recalculate_alleles_overlaps(self, request, queryset: QuerySet[Overlap]):
        allele_ids = queryset.order_by('allele').distinct('allele').values_list('allele', flat=True)
        alleles = Allele.objects.filter(pk__in=allele_ids)
        for allele in alleles:
            OverlapServices().calculate_and_apply_overlaps_for_allele(allele)
        #     update_overlaps_for_allele(allele)

    @admin_list_column(short_description="Modified", order_field="modified")
    def modified_detailed(self, obj: Overlap):
        return self.format_datetime(obj.modified)

    @admin_list_column(short_description="overlap_status_display", order_field="overlap_status")
    def overlap_status_display(self, obj: Overlap):
        return obj.get_overlap_status_display()