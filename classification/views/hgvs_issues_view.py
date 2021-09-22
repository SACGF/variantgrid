from dataclasses import dataclass
from datetime import datetime
from typing import Dict, Any, Optional, Tuple, List
from django.contrib.auth.decorators import user_passes_test
from django.db.models import QuerySet
from django.http import StreamingHttpResponse
from django.http.request import HttpRequest
from django.shortcuts import render
from django.urls import reverse
from lazy import lazy
from requests.models import Response

from classification.enums import SpecialEKeys
from classification.models import Classification, classification_flag_types
from flags.models import Flag, FlagStatus, FlagType
from flags.models.models import FlagCollection
from genes.hgvs import CHGVS, chgvs_diff_description
from library.django_utils import get_url_from_view_path
from library.guardian_utils import is_superuser
from library.utils import ExportRow, export_column
from snpdb.models import VariantAllele, allele_flag_types, GenomeBuild, Variant, ClinGenAllele
from snpdb.models.models_variant import Allele
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder


def _alleles_with_classifications_qs() -> QuerySet[Allele]:
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
    context = {
        "counts": issue_counts
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


@dataclass(frozen=True)
class FlagReport(ExportRow):
    flag: Flag

    @lazy
    def data(self) -> Tuple[str, datetime]:
        """
        Returns the username of who closed the flag (if the flag is closed)
        And then the closed date of the flag, or created date if not closed
        """
        if last_status_comment := self.flag.flagcomment_set.order_by('-created').filter(resolution__isnull=False).first():
            if last_status_comment.resolution.status == FlagStatus.CLOSED:
                return last_status_comment.created, last_status_comment.user.username
            else:
                return self.flag.created, None

    @export_column("Date")
    def date(self) -> str:
        if date := self.data[0]:
            return date.strftime("%Y-%m-%d")

    @export_column("Closed By")
    def username(self) -> str:
        return self.data[1]

    @export_column("Note")
    def note(self) -> str:
        return self.flag.flagcomment_set.order_by('-created').first().text


@dataclass(frozen=True)
class ProblemHgvs(ExportRow):
    classification: Classification

    def flag_formatter(self, flag_type: FlagType, data: Dict[str, Any] = None):
        qs: QuerySet[Flag]
        if flag_type.context_id == 'classification':
            qs = self.classification.flags_of_type(flag_type=flag_type)
        elif flag_type.context_id == 'allele':
            if allele := self.allele:
                qs = allele.flags_of_type(flag_type=flag_type)
            else:
                return
        else:
            raise ValueError(f"Unexpected flag context {flag_type.context_id}")

        qs = qs.order_by('-created')
        if data:
            for key, value in data.items():
                qs = qs.filter(**{f'data__{key}': value})

        if last_flag := qs.first():
            return FlagReport(last_flag)

    @lazy
    def allele(self):
        if variant := self.classification.variant:
            return variant.allele

    @export_column("Classification URL")
    def classification_url(self):
        return get_url_from_view_path(self.classification.get_absolute_url())

    @export_column("Withdrawn")
    def withdrawn(self):
        return self.classification.withdrawn

    @export_column("Imported Build")
    def imported_build(self):
        return self.classification.get(SpecialEKeys.GENOME_BUILD)

    @export_column("Imported 37")
    def imported_c_hgvs_37(self):
        if (imported_build := self.imported_build()) and '37' in imported_build:
            return self.classification.get(SpecialEKeys.C_HGVS)

    @export_column("Resolved 37")
    def resolved_37(self):
        return self.classification.chgvs_grch37

    @export_column("Imported 38")
    def imported_c_hgvs_38(self):
        if (imported_build := self.imported_build()) and '38' in imported_build:
            return self.classification.get(SpecialEKeys.C_HGVS)

    @export_column("Resolved 38")
    def resolved_38(self):
        return self.classification.chgvs_grch38

    @export_column("Flag Matching", sub_data=FlagReport)
    def flag_variant_matching(self):
        return self.flag_formatter(classification_flag_types.matching_variant_flag)

    @export_column("Flag Matching Warning", sub_data=FlagReport)
    def flag_matching_warning(self):
        return self.flag_formatter(classification_flag_types.matching_variant_warning_flag)

    @export_column("Transcript Change", sub_data=FlagReport)
    def flag_transcript_change(self):
        return self.flag_formatter(classification_flag_types.transcript_version_change_flag)

    @export_column("Allele URL")
    def allele_url(self):
        if allele := self.allele:
            return get_url_from_view_path(allele.get_absolute_url())

    @export_column("37 != 38", sub_data=FlagReport)
    def flag_37_not_38(self):
        return self.flag_formatter(allele_flag_types.allele_37_not_38, data={"transcript": self.classification.transcript})

    @export_column("!37", sub_data=FlagReport)
    def flag_37(self):
        return self.flag_formatter(allele_flag_types.missing_37)

    @export_column("!38", sub_data=FlagReport)
    def flag_38(self):
        return self.flag_formatter(allele_flag_types.missing_38)


@dataclass(frozen=True)
class ClassificationResolution(ExportRow):
    _imported_genome_build: str
    _chgvs_imported: str
    _chgvs_grch37: str
    _chgvs_grch38: str
    _variant_id: Optional[int]
    _url: str

    @export_column("Import Key")
    def key_import(self):
        return (self._imported_genome_build or "") + "#" + (self._chgvs_imported or "")

    @export_column("Export Key")
    def key_export(self):
        return (self._chgvs_grch37 or "") + "#" + (self._chgvs_grch38 or "")

    @export_column("URL")
    def url(self):
        return self._url

    @export_column("Imported Build")
    def imported_build(self):
        return self._imported_genome_build

    @export_column("c.HGVS Imported")
    def imported_chgvs(self):
        return self._chgvs_imported

    @export_column("c.HGVS 37")
    def resolved_37(self):
        return self._chgvs_grch37

    @export_column("c.HGVS 38")
    def resolved_38(self):
        return self._chgvs_grch38

    @property
    def is_valid(self):
        return self._imported_genome_build is not None and self._chgvs_imported is not None

    @staticmethod
    def diff_description(c_hgvs1: str, c_hgvs2: str):
        if not c_hgvs1 or not c_hgvs2:
            return "Failed"
        if c_hgvs1 != c_hgvs2 and (diff_descripts := chgvs_diff_description(CHGVS(c_hgvs1).diff(CHGVS(c_hgvs2)))):
            return ", ".join(diff_descripts)
        return ""  # blank is going to make it easier to view the spreadsheet

    @export_column("Normal/Liftover")
    def diffs(self):

        normal = ""
        if not self.is_valid:
            normal = "Invalid"
        else:
            if "38" in self.imported_build():
                normalised_c = self._chgvs_grch38
                lifted_over_c = self._chgvs_grch37
            else:
                normalised_c = self._chgvs_grch37
                lifted_over_c = self._chgvs_grch38

            normal = ClassificationResolution.diff_description(self._chgvs_imported, normalised_c)
            liftover = ClassificationResolution.diff_description(self._chgvs_imported, lifted_over_c)

            if normal or liftover:
                return f"{normal}#{liftover}"
            else:
                return ""

    @lazy
    def allele(self) -> Optional[Allele]:
        if variant_id := self._variant_id:
            v = Variant.objects.get(pk=variant_id)
            return v.allele
        return None

    @export_column("Variant")
    def variant(self):
        allele: Allele  # @lazy screws up type hints :(
        if allele := self.allele:
            parts: List[str] = list()
            clingen_id = ""
            if clingen_allele_id := allele.clingen_allele_id:
                clingen_id = ClinGenAllele.format_clingen_allele(clingen_allele_id)

            parts.append(clingen_id or "")
            for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
                coord = ""
                try:
                    v = allele.variant_for_build(genome_build)
                    coord = str(v)
                except ValueError:
                    pass
                parts.append(coord)
            return "#".join(parts)
        else:
            return ""


@user_passes_test(is_superuser)
def download_hgvs_issues(request: HttpRequest) -> StreamingHttpResponse:

    # source data as all classifications with a flag
    # and all classifications attached to an allele with a flag
    alleles_qs = FlagCollection.filter_for_flags(Allele.objects.all())
    classification_qs = FlagCollection.filter_for_flags(Classification.objects.all(), flag_types=[
        classification_flag_types.matching_variant_warning_flag,
        classification_flag_types.transcript_version_change_flag,
        classification_flag_types.matching_variant_flag
    ])
    complete_qs = classification_qs.union(Classification.objects.filter(variant__variantallele__allele__in=alleles_qs))
    complete_qs = complete_qs.order_by('-variant', '-pk')

    return ProblemHgvs.streaming(request, complete_qs, "hgvs_issues")


@user_passes_test(is_superuser)
def download_hgvs_resolution(request: HttpRequest) -> StreamingHttpResponse:

    imported_genome_build_col = 'evidence__genome_build__value'
    c_hgvs_imported_col = 'evidence__c_hgvs__value'
    c_hgvs_37_col = 'chgvs_grch37'
    c_hgvs_38_col = 'chgvs_grch38'

    qs = Classification.objects.order_by(imported_genome_build_col, c_hgvs_imported_col, c_hgvs_37_col, c_hgvs_38_col, 'variant', 'pk').values_list(
        imported_genome_build_col, c_hgvs_imported_col, c_hgvs_37_col, c_hgvs_38_col, 'variant', 'pk'
    )

    last_record = list()
    def mapper(data):
        record = [data[0], data[1], data[2], data[3]]
        if record != last_record and data[0] and data[1]:
            variant_id = data[4]
            classification_id = data[5]
            url = get_url_from_view_path(
                reverse('view_classification', kwargs={'record_id': data[5]}),
            )
            # url = f"https://shariant.org.au/classification/classification/{data[4]}"
            return ClassificationResolution(
                _imported_genome_build=data[0],
                _chgvs_imported=data[1],
                _chgvs_grch37=data[2],
                _chgvs_grch38=data[3],
                _variant_id=variant_id,
                _url=url
            )

    stream = (mapper(row) for row in qs)
    return ClassificationResolution.streaming(request, stream, "c_hgvs_resolution")