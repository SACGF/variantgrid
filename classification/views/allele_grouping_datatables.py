from functools import cached_property
from typing import Any

from django.db.models import QuerySet
from django.http import HttpRequest

from classification.models import AlleleGrouping, ClassificationGrouping
from genes.hgvs import CHGVS
from library.cache import timed_cache
from library.utils import JsonDataType
from snpdb.models import GenomeBuild, UserSettings
from snpdb.views.datatable_view import DatatableConfig, RichColumn, CellData


@timed_cache(size_limit=1, ttl=60)
def _classification_groups_for_allele_group(allele_group_id: int) -> list[ClassificationGrouping]:
    return list(ClassificationGrouping.objects.filter(allele_grouping__pk=allele_group_id))


class AlleleGroupingColumns(DatatableConfig[AlleleGrouping]):

    @cached_property
    def genome_build_prefs(self) -> list[GenomeBuild]:
        user_settings = UserSettings.get_for_user(self.user)
        return GenomeBuild.builds_with_annotation_priority(user_settings.default_genome_build)

    def get_initial_queryset(self) -> QuerySet[AlleleGrouping]:
        # TODO, consider making the groups GuardianPermission rather than this manual security check
        return AlleleGrouping.objects.all()

    def c_hgvs_for(self, cg: ClassificationGrouping) -> CHGVS:
        is_preferred_genome_build = True
        allele_info = cg.latest_classification.classification.allele_info
        for gb in self.genome_build_prefs:
            if ri := allele_info[gb]:
                if ri.c_hgvs_obj:
                    return ri.c_hgvs_obj
            is_preferred_genome_build = False

        c_hgvs = allele_info.imported_c_hgvs_obj
        c_hgvs.is_normalised = False
        return c_hgvs

    def render_allele(self, row: CellData) -> JsonDataType:
        cgs = _classification_groups_for_allele_group(row.get("allele_id"))
        all_chgvs = list(sorted({self.c_hgvs_for(cg) for cg in cgs}))
        return all_chgvs[0].to_json()
        # format_hgvs

    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(key="allele", extra_columns=["allele_id"], renderer=self.c_hgvs_for, client_renderer='VCTable.hgvs'),
            RichColumn(key="overlap_status")
        ]