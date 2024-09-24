import json
from functools import cached_property
from typing import Any, Optional

from django.db.models import QuerySet, OuterRef, Subquery
from django.http import HttpRequest

from classification.enums import AlleleOriginBucket
from classification.models import AlleleOriginGrouping, ClassificationGrouping, OverlapStatus, AlleleGrouping
from classification.templatetags.classification_tags import classification
from genes.hgvs import CHGVS
from library.cache import timed_cache
from library.utils import JsonDataType
from snpdb.lab_picker import LabPickerData
from snpdb.models import GenomeBuild, UserSettings, Lab
from snpdb.views.datatable_view import DatatableConfig, RichColumn, CellData, SortOrder


@timed_cache(size_limit=1, ttl=60)
def _allele_group(allele_id: int) -> AlleleGrouping:
    return AlleleGrouping.objects.get(allele=allele_id)


class AlleleGroupingColumns(DatatableConfig[AlleleGrouping]):

    @cached_property
    def genome_build_prefs(self) -> list[GenomeBuild]:
        user_settings = UserSettings.get_for_user(self.user)
        return GenomeBuild.builds_with_annotation_priority(user_settings.default_genome_build)

    @cached_property
    def lab_picker(self) -> LabPickerData:
        return LabPickerData.for_user(user=self.user, selection=self.get_query_param("lab_id"))

    def get_initial_queryset(self) -> QuerySet[AlleleGrouping]:
        # TODO, consider making the groups GuardianPermission rather than this manual security check
        qs = AlleleGrouping.objects.all()

        germline_subquery = AlleleOriginGrouping.objects.filter(
            allele_grouping=OuterRef('pk'),
            allele_origin_bucket=AlleleOriginBucket.GERMLINE
        ).values_list("overlap_status", flat=True)
        somatic_subquery = AlleleOriginGrouping.objects.filter(
            allele_grouping=OuterRef('pk'),
            allele_origin_bucket=AlleleOriginBucket.SOMATIC
        ).values_list("overlap_status", flat=True)

        qs = qs.annotate(germline_overlap_status=Subquery(germline_subquery))
        qs = qs.annotate(somatic_overlap_status=Subquery(somatic_subquery))

        lab_picker = self.lab_picker
        if not lab_picker.is_admin_mode:
            linked_class_groupings = ClassificationGrouping.objects.filter(
                lab__in=lab_picker.selected_labs
            ).values_list('allele_origin_grouping__allele_grouping', flat=True)
            qs = qs.filter(pk__in=linked_class_groupings)

        return qs

    def render_labs(self, row: CellData) -> set[Lab]:
        ag = _allele_group(row.get("allele"))
        labs: set[Lab] = set()
        for ag in ag.allele_origin_dict.values():
            for cg in ag.classificationgrouping_set.all():
                labs.add(cg.lab)
        labs_list = list(sorted(labs))
        return "".join("<div>" + str(lab) +"</div>" for lab in labs_list)

    def c_hgvs_for(self, cg: ClassificationGrouping) -> CHGVS:
        is_preferred_genome_build = True
        allele_info = cg.latest_classification_modification.classification.allele_info
        for gb in self.genome_build_prefs:
            if ri := allele_info[gb]:
                if c_hgvs := ri.c_hgvs_obj:
                    c_hgvs.is_desired_build = is_preferred_genome_build
                    return c_hgvs
            is_preferred_genome_build = False

        c_hgvs = allele_info.imported_c_hgvs_obj
        c_hgvs.is_normalised = False
        return c_hgvs

    def render_allele(self, row: CellData) -> JsonDataType:
        allele_group = _allele_group(row.get("allele"))
        # FIXME cache this
        cgs = ClassificationGrouping.objects.filter(allele_origin_grouping__allele_grouping=allele_group.pk, dirty=False)
        all_chgvs = list(sorted({self.c_hgvs_for(cg) for cg in cgs}))
        c_hgvs_json: JsonDataType
        if all_chgvs:
            c_hgvs_json = all_chgvs[0].to_json()
        else:
            c_hgvs_json = {}
        c_hgvs_json["allele_id"] = row.get("allele")
        return c_hgvs_json
        # format_hgvs

    def _render_bucket_status(self, row: CellData, bucket: AlleleOriginBucket) -> JsonDataType:
        aog = _allele_group(row.get("allele"))
        if germline := aog.allele_origin_grouping(bucket):
            return OverlapStatus(germline.overlap_status)
        else:
            # FIXME make No shared records status
            return OverlapStatus.NO_SHARED_RECORDS

    def render_germline_status(self, row: CellData) -> JsonDataType:
        germline_status = self._render_bucket_status(row, AlleleOriginBucket.GERMLINE)
        aog = _allele_group(row.get("allele"))
        classification_values: Optional[list[str]] = None
        if germline := aog.allele_origin_grouping(AlleleOriginBucket.GERMLINE):
            classification_values = germline.classification_values

        return [germline_status, classification_values]

    def render_somatic_status(self, row: CellData) -> JsonDataType:
        somatic_status = self._render_bucket_status(row, AlleleOriginBucket.SOMATIC)
        aog = _allele_group(row.get("allele"))
        classification_values: Optional[list[str]] = None
        if somatic := aog.allele_origin_grouping(AlleleOriginBucket.SOMATIC):
            classification_values = somatic.somatic_clinical_significance_values
        return [somatic_status, classification_values]


    # def render_details(self, row: CellData) -> JsonDataType:
    #     cgs = _allele_group(row.get("id"))
    #     return json.dumps([cg.to_json() for cg in cgs])

    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.expand_client_renderer = DatatableConfig._row_expand_ajax('allele_grouping_detail',
                                                                       expected_height=108)

        self.rich_columns = [
            RichColumn(key="allele", renderer=self.render_allele, client_renderer='VCTable.hgvs'),
            RichColumn(
                name="germline_overlap",
                renderer=self.render_germline_status,
                client_renderer=RichColumn.client_renderer_combine([
                    RichColumn.choices_client_renderer(OverlapStatus.choices),
                    RichColumn.client_renderer_repeat({"formatter": 'VCTable.classification'})
                ]),
                order_sequence=[SortOrder.DESC, SortOrder.ASC],
                default_sort=SortOrder.DESC,
                sort_keys=["germline_overlap_status"],
                extra_columns=["allele"]
            ),
            RichColumn(
                name="somatic_overlap",
                renderer=self.render_somatic_status,
                client_renderer=RichColumn.client_renderer_combine([
                    RichColumn.choices_client_renderer(OverlapStatus.choices),
                    RichColumn.client_renderer_repeat({"formatter": "VCTable.somatic_clinical_significance"})
                ]),
                order_sequence=[SortOrder.DESC, SortOrder.ASC],
                sort_keys=["somatic_overlap_status"],
                extra_columns=["allele"]
            ),
            RichColumn(name="labs", renderer=self.render_labs, extra_columns=["allele"]),
            RichColumn(key="id", visible=False)
            # RichColumn(name="details", label="Details", extra_columns=["allele"], renderer=self.render_details),
        ]