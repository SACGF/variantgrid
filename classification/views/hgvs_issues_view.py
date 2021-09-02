from typing import Dict, Any, Optional

from django.conf import settings
from django.contrib.auth.decorators import user_passes_test
from django.db.models import QuerySet
from django.http import StreamingHttpResponse
from django.http.request import HttpRequest
from django.shortcuts import render
from pytz import timezone
from django.utils.timezone import now
from requests.models import Response

from classification.models import Classification
from flags.models import Flag
from flags.models.models import FlagCollection
from library.guardian_utils import is_superuser
from library.utils import delimited_row
from snpdb.models import VariantAllele, allele_flag_types, GenomeBuild, Variant
from snpdb.models.models_variant import Allele
from classification.models.flag_types import classification_flag_types
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


def _alleles_with_classifications_qs() -> QuerySet:
    """
    Only consider alleles that have variant classifications
    as there are way too many alleles to consider otherwise
    """
    classification_variant_ids_qs = Classification.objects.exclude(withdrawn=True).values_list('variant__id', flat=True)
    allele_ids = VariantAllele.objects.filter(variant_id__in=classification_variant_ids_qs).values_list('allele_id', flat=True)
    return Allele.objects.filter(id__in=allele_ids)


@user_passes_test(is_superuser)
def view_hgvs_issues(request: HttpRequest) -> Response:

    alleles_qs = _alleles_with_classifications_qs().order_by('-id')
    allele_missing_rep_qs = FlagCollection.filter_for_open_flags(qs=alleles_qs, flag_types=[allele_flag_types.missing_37, allele_flag_types.missing_38])

    allele_missing_clingen = alleles_qs.filter(clingen_allele__isnull=True)
    allele_37_not_38 = FlagCollection.filter_for_open_flags(qs=alleles_qs, flag_types=[allele_flag_types.allele_37_not_38])

    vcqs = Classification.objects.filter(withdrawn=False)
    matching_variant = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[classification_flag_types.matching_variant_flag])
    matching_variant_warning = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[classification_flag_types.matching_variant_warning_flag])
    matching_variant_transcript = FlagCollection.filter_for_open_flags(qs=vcqs, flag_types=[classification_flag_types.transcript_version_change_flag])

    issue_counts = {
        "allele": {
            "missing_rep": allele_missing_rep_qs.count(),
            "missing_clingen": allele_missing_clingen.count(),
            "chgvs_37_not_38": allele_37_not_38.count()
        },
        "classifications": {
            "matching_variant": matching_variant.count(),
            "matching_variant_warning": matching_variant_warning.count(),
            "matching_variant_transcript": matching_variant_transcript.count()
        }
    }
    from classification.views.classification_datatables import ClassificationColumns
    context = {
        "counts": issue_counts,
        "allele_columns": AlleleColumns(request),
        "classification_columns": ClassificationColumns(request)
    }

    return render(request, "classification/hgvs_issues.html", context)


class AlleleColumns(DatatableConfig):

    def get_allele(self, allele_id: int) -> Allele:
        if last_allele := self.last_allele:
            if last_allele.id == allele_id:
                return last_allele
        last_allele = Allele.objects.get(id=allele_id)
        return last_allele

    def variant_for(self, row: Dict[str, Any], genome_build: GenomeBuild) -> str:
        values = list()
        allele = self.get_allele(row["id"])
        variant: Variant
        try:
            variant = allele.variant_for_build(genome_build=genome_build, best_attempt=False)
        except ValueError:
            return "-"
        values.append(f'<span class="text-secondary">{variant}</span>')
        if variant_annotation := variant.get_best_variant_transcript_annotation(genome_build):
            if gene := variant_annotation.gene:
                gene_symbol = gene.get_gene_symbol(genome_build)
                values.append(gene_symbol.symbol)

        return "<br/>".join(values)

    def variant_37(self, row: Dict[str, Any]) -> Optional[str]:
        return self.variant_for(row, GenomeBuild.grch37())

    def variant_38(self, row: Dict[str, Any]) -> Optional[str]:
        return self.variant_for(row, GenomeBuild.grch38())

    def __init__(self, request):
        super().__init__(request)
        self.last_allele: Optional[Allele] = None

        self.rich_columns = [
            RichColumn(key="id", label='ID', client_renderer='alleleIdRender', orderable=True, default_sort=SortOrder.DESC),
            RichColumn(key="clingen_allele__id", client_renderer='clingenIdRenderer', label='ClinGen Allele', orderable=True),
            RichColumn(name="variant_37", label='37 Variant', renderer=self.variant_37, orderable=False),
            RichColumn(name="variant_38", label='38 Variant', renderer=self.variant_38, orderable=False),
            RichColumn(key="flag_collection_id", label="Flags", client_renderer='TableFormat.flags')
        ]

    def get_initial_queryset(self):
        # exclude where we've auto matched and have 0 outstanding left
        return _alleles_with_classifications_qs()

    def filter_queryset(self, qs: QuerySet) -> QuerySet:
        if flag_filter := self.get_query_param('flag'):
            if flag_filter == 'missing_rep':
                qs = FlagCollection.filter_for_open_flags(qs=qs, flag_types=[allele_flag_types.missing_37,
                                                                             allele_flag_types.missing_38])
            elif flag_filter == 'missing_clingen':
                qs = qs.filter(clingen_allele__isnull=True)
            elif flag_filter == 'chgvs_37_not_38':
                qs = FlagCollection.filter_for_open_flags(qs=qs, flag_types=[allele_flag_types.allele_37_not_38])

        return qs


@user_passes_test(is_superuser)
def download_hgvs_issues(request: HttpRequest) -> StreamingHttpResponse:

    def row_generator():
        yield delimited_row([
            "Allele ID",
            "Allele URL",
            "Clingen Allele ID",
            "Gene",
            "GRCh37",
            "GRCh38",
            "Issue Type",
            "Text"
        ])
        alleles_qs = FlagCollection.filter_for_open_flags(qs=_alleles_with_classifications_qs())
        allele: Allele
        for allele in alleles_qs:
            open_flags = Flag.objects.filter(collection=allele.flag_collection).filter(FlagCollection.Q_OPEN_FLAGS)
            flag: Flag
            for flag in open_flags:

                gene_symbol_str = None
                for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
                    try:
                        variant = allele.variant_for_build(genome_build, best_attempt=False)
                        if variant_annotation := variant.get_best_variant_transcript_annotation(genome_build):
                            if gene := variant_annotation.gene:
                                if gene_symbol := gene.get_gene_symbol(genome_build):
                                    gene_symbol_str = gene_symbol.symbol
                                    break
                    except ValueError:
                        pass

                yield delimited_row([
                    allele.id,
                    request.build_absolute_uri(allele.get_absolute_url()),
                    allele.clingen_allele_id,
                    gene_symbol_str,
                    str(allele.grch37) if allele.grch37 else "",
                    str(allele.grch38) if allele.grch38 else "",
                    flag.flag_type.label,
                    flag.flagcomment_set.first().text
                ])

    response = StreamingHttpResponse(row_generator(), content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename="allele_issues_{now().astimezone(tz=timezone(settings.TIME_ZONE)).strftime("%Y-%m-%d")}.csv"'
    return response
