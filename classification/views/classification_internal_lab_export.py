import operator
import re
from functools import reduce
from typing import Optional

from django.contrib import messages
from django.db.models import Q
from django.shortcuts import redirect, render
from django.urls import reverse

from classification.enums import ShareLevel

from classification.views.exports import ClassificationExportFormatterCSV
from classification.views.exports.classification_export_filter import ClassificationFilter, \
    classification_export_user_string_to_q
from classification.views.exports.classification_export_formatter_csv import FormatDetailsCSV

from snpdb.models import GenomeBuild, Lab, AlleleOriginFilterDefault


def internal_lab_download(request):
    if request.method == 'POST':
        user = request.user
        share_level = request.POST.get('share_level')
        build = request.POST.get('genome_build')
        allele_origin = AlleleOriginFilterDefault(request.POST.get('allele-origin-toggle'))
        record_filter_str = request.POST.get('record_filters')
        genomic_build = GenomeBuild.get_from_fuzzy_string(build)
        form_data = {
            'share_level': share_level,
            'genome_build': build,
            'allele_origin': allele_origin,
            'record_filters': record_filter_str,
            'active_tab': 'internal_download:texport-my-data'
        }
        share_level_map = {
            'any': ShareLevel.ALL_USERS,
            'public': ShareLevel.PUBLIC
        }
        share_level = share_level_map.get(share_level, share_level)
        context = {
            "active_tab": "internal_download:texport-my-data",
            "form_data": form_data
        }
        try:
            user_labs = set(Lab.valid_labs_qs(user, admin_check=True))
            if not user_labs:
                messages.error(request, 'You are not associated with any labs.')
                return render(request, 'classification/classification_export.html', context)

            extra_filter_qs: Optional[Q] = None
            if record_filter_str:
                all_qs: list[Q] = []
                has_errors = False
                for part in re.split('[,\n\t]', record_filter_str):
                    # handle each part here so we can validate that there's 1+ record found for each filter
                    part = part.strip()
                    try:
                        part_q = classification_export_user_string_to_q(part, build)
                        filter_data = ClassificationFilter(
                            user=user,
                            genome_build=genomic_build,
                            allele_origin_filter=allele_origin,
                            min_share_level=share_level,
                            file_prefix=f"Internal_lab_report",
                            include_sources=user_labs,
                            extra_filter=part_q
                        )
                        if not filter_data.cms_qs.exists():
                            messages.error(request, f"\"{part}\" - no results found")
                            has_errors = True
                        else:
                            all_qs.append(part_q)

                    except ValueError as ve:
                        has_errors = True
                        messages.error(request, f"\"{part}\" - could not be turned into a classification filter")

                if has_errors:
                    return render(request, 'classification/classification_export.html', context)

                if all_qs:
                    extra_filter_qs = reduce(operator.or_, all_qs)

            filter_data = ClassificationFilter(
                user=user,
                genome_build=genomic_build,
                allele_origin_filter=allele_origin,
                min_share_level=share_level,
                file_prefix=f"Internal_lab_report",
                include_sources=user_labs,
                extra_filter=extra_filter_qs
            )

            response = ClassificationExportFormatterCSV(
                filter_data,
                FormatDetailsCSV(exclude_discordances=True, exclude_transient=True)
            ).serve()
            return response
        except Exception as e:
            print(e)
            messages.error(request, 'Unable to generate report')
            return render(request, 'classification/classification_export.html', context)
    else:
        return redirect(reverse('classification_export') + '?activeTab=internal_download%3Atexport-my-data')
