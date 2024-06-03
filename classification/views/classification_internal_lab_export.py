from django.contrib import messages
from django.http import JsonResponse
from django.shortcuts import redirect, render
from django.urls import reverse

from classification.enums import ShareLevel

from classification.views.exports import ClassificationExportFormatterCSV
from classification.views.exports.classification_export_filter import ClassificationFilter
from classification.views.exports.classification_export_formatter_csv import FormatDetailsCSV

from snpdb.models import GenomeBuild, Lab, AlleleOriginFilterDefault


def internal_lab_download(request):
    if request.method == 'POST':
        user = request.user
        share_level = request.POST.get('share_level')
        build = request.POST.get('genome_build')
        allele_origin = AlleleOriginFilterDefault(request.POST.get('allele-origin-toggle'))
        genomic_locations = request.POST.get('record_filters')
        genomic_build = GenomeBuild.get_from_fuzzy_string(build)
        form_data = {
            'share_level': share_level,
            'genome_build': build,
            'allele_origin': allele_origin,
            'record_filters': genomic_locations,
            'active_tab': 'internal_download:texport-my-data'
        }
        share_level_map = {
            'any': ShareLevel.ALL_USERS,
            'public': ShareLevel.PUBLIC
        }
        share_level = share_level_map.get(share_level, share_level)
        try:
            user_labs = set(Lab.valid_labs_qs(user, admin_check=True))
            if not user_labs:
                messages.error(request, 'You are not associated with any labs.')
                return render(request, 'classification/classification_export.html',
                              {"active_tab": 'internal_download:texport-my-data', "form_data": form_data})
            filter_data = ClassificationFilter(
                    user=user,
                    genome_build=genomic_build,
                    allele_origin_filter=allele_origin,
                    min_share_level=share_level,
                    record_filters=genomic_locations,
                    file_prefix=f"Internal_lab_report",
                    include_sources=user_labs,
                )
            if filter_data.cms_qs.count() == 0:
                messages.error(request, 'No records found for the selected filters.')
                return render(request, 'classification/classification_export.html',
                              {"active_tab": 'internal_download:texport-my-data', 'form_data': form_data})
            if filter_data.excluded_record_filters:
                excluded_records = ', '.join(filter_data.excluded_record_filters)
                messages.error(request, f'No records found for the filters: {excluded_records}')
                return render(request, 'classification/classification_export.html',
                              {"active_tab": 'internal_download:texport-my-data', 'form_data': form_data})
            response = ClassificationExportFormatterCSV(
                filter_data,
                FormatDetailsCSV(exclude_discordances=True, exclude_transient=True)
            ).serve()
            return response
        except Exception as e:
            print(e)
            messages.error(request, 'Unable to generate report')
            return render(request, 'classification/classification_export.html',
                          {"active_tab": 'internal_download:texport-my-data', 'form_data': form_data})
    else:
        return redirect(reverse('classification_export') + '?activeTab=internal_download%3Atexport-my-data')
