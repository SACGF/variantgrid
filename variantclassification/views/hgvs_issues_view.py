import csv
from dataclasses import dataclass
from typing import List, io

from django.conf import settings
from django.contrib.auth.decorators import user_passes_test
from django.db.models import QuerySet
from django.http import StreamingHttpResponse
from django.http.request import HttpRequest
from django.shortcuts import render
from pytz import timezone
from django.utils.timezone import now
from requests.models import Response

from flags.models import Flag
from flags.models.models import FlagCollection, FlagTypeContext
from library.django_utils import get_url_from_view_path
from library.guardian_utils import is_superuser
from library.utils import delimited_row
from snpdb.models import VariantAllele
from snpdb.models.models_variant import Allele
from variantclassification.models.flag_types import variant_classification_flag_types
from variantclassification.models.variant_classification import VariantClassification,\
    VariantClassificationModification

# Only consider alleles that have variant classifications
# as there are way too many alleles to consider otherwise
def _alleles_with_variants_qs() -> QuerySet:
    # find the variant for ALL variant classifications, and keep a dict of variant id to classification id
    classification_variant_ids_qs = VariantClassification.objects.exclude(withdrawn=True).values_list('variant__id', flat=True)
    allele_ids = VariantAllele.objects.filter(variant_id__in=classification_variant_ids_qs).values_list('allele_id', flat=True)
    return Allele.objects.filter(id__in=allele_ids)


@user_passes_test(is_superuser)
def view_hgvs_issues(request: HttpRequest) -> Response:

    vcqs = VariantClassification.objects
    flagged_classifications = FlagCollection.filter_for_open_flags(
        qs=vcqs,
        flag_types=[
            variant_classification_flag_types.matching_variant_flag,
            variant_classification_flag_types.matching_variant_warning_flag,
            variant_classification_flag_types.transcript_version_change_flag
        ]).filter(withdrawn=False)

    flagged_vcms = VariantClassificationModification.objects.filter(
        variant_classification__in=flagged_classifications,
        is_last_published=True
    ).select_related('variant_classification', 'variant_classification__lab')

    variant_alleles_qs = _alleles_with_variants_qs()
    alleles_qs = FlagCollection.filter_for_open_flags(qs=variant_alleles_qs)

    # Only show missing ClinGen alleles if that's our only option
    # for now always show alleles missing Clingen Allele ID so we can easily review NCBI lifted over variants
    show_missing_clingen_alleles = True  # settings.LIFTOVER_NCBI_REMAP_ENABLED is False
    if show_missing_clingen_alleles:
        missing_clingen_alleles = variant_alleles_qs.filter(clingen_allele_id__isnull=True)
        alleles_qs = alleles_qs.union(missing_clingen_alleles)
    alleles_qs = alleles_qs.order_by('id').distinct()

    context = {
        "classifications": flagged_vcms,
        "alleles": alleles_qs,
        "filter_flags": ' '.join([
            variant_classification_flag_types.matching_variant_flag.id,
            variant_classification_flag_types.matching_variant_warning_flag.id,
            variant_classification_flag_types.transcript_version_change_flag.id,
        ] +
        [ft.id for ft in FlagTypeContext.objects.get(pk='allele').flagtype_set.all()]
    )
    }
    return render(request, "variantclassification/hgvs_issues.html", context)


@user_passes_test(is_superuser)
def download_hgvs_issues(request: HttpRequest) -> Response:

    def row_generator():
        yield delimited_row([
            "Allele ID",
            "Allele URL",
            "Clingen Allele ID",
            "GRCh37",
            "GRCh38",
            "Issue Type",
            "Text"
        ])
        alleles_qs = FlagCollection.filter_for_open_flags(qs=_alleles_with_variants_qs())
        allele: Allele
        for allele in alleles_qs:
            open_flags = Flag.objects.filter(collection=allele.flag_collection).filter(FlagCollection.Q_OPEN_FLAGS)
            flag: Flag
            for flag in open_flags:

                yield delimited_row([
                    allele.id,
                    get_url_from_view_path(allele.get_absolute_url()),
                    allele.clingen_allele_id,
                    str(allele.grch37),
                    str(allele.grch38),
                    flag.flag_type.label,
                    flag.flagcomment_set.first().text
                ])

    response = StreamingHttpResponse(row_generator(), content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename="allele_issues_{now().astimezone(tz=timezone(settings.TIME_ZONE)).strftime("%Y-%m-%d")}.csv"'
    return response
