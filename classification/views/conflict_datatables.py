from dataclasses import dataclass
from functools import cached_property
from typing import Optional

from django.db.models import QuerySet
from django.http import HttpRequest
from django.template.loader import render_to_string

from classification.enums import ConflictSeverity, TestingContextBucket, AlleleOriginBucket, ConflictType
from classification.models import Conflict, ConflictLab, DiscordanceReportTriageStatus
from classification.services.conflict_services import ConflictDataRow
from genes.hgvs import CHGVS
from snpdb.lab_picker import LabPickerData
from snpdb.models import Allele, Lab
from snpdb.views.datatable_view import DatatableConfig, RichColumn, DatatableConfigQuerySetMode, SortOrder, CellData, DC


@dataclass(frozen=True)
class ConflictRowWithStatus:
    conflict_row: ConflictDataRow
    status: DiscordanceReportTriageStatus

class ConflictColumns(DatatableConfig[Conflict]):

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        # self.expand_client_renderer = DatatableConfig._row_expand_ajax('conflict_detail', expected_height=108)

        self.allele_map: dict[int, Allele] = {}

        self.rich_columns = [
            RichColumn(
                name="context",
                label="Context",
                extra_columns=["allele_origin_bucket", "testing_context_bucket", "conflict_type"],
                sort_keys=["allele_origin_bucket", "testing_context_bucket", "conflict_type"],
                renderer=self.render_context,
                client_renderer='renderContext',
            ),
            RichColumn(
                name="c_hgvs",
                label="c.HGVS",
                extra_columns=["allele", "meta_data"],
                renderer=self.render_c_hgvs,
                client_renderer=RichColumn.client_renderer_repeat({"formatter": 'VCTable.hgvs', "wrapperCSS": "display-block vertical-sep"})
            ),
            RichColumn(
                "severity",
                renderer=self.render_severity,
                extra_columns=["pk"],
                order_sequence=[SortOrder.DESC, SortOrder.ASC],
                default_sort=SortOrder.DESC,
                client_renderer='renderSeverity',
            ),
            RichColumn(name="involved_labs", label="Involved Labs", extra_columns=["data", "conflict_type", "pk"], renderer=self.involved_labs),
            RichColumn("id", visible=False)
        ]

    @cached_property
    def lab_map(self) -> dict[int, Lab]:
        all_labs = Lab.objects.all()
        return {lab.pk: lab for lab in all_labs}

    def pre_render(self, qs: QuerySet[DC]):
        alleles = Allele.objects.filter(pk__in=qs.values_list("allele__pk", flat=True))
        self.allele_map = {allele.pk: allele for allele in alleles}

    def render_context(self, row_data: CellData) -> dict:
        allele_origin_bucket = AlleleOriginBucket(row_data.get("allele_origin_bucket"))
        testing_context_bucket = TestingContextBucket(row_data.get("testing_context_bucket"))
        conflict_type = ConflictType(row_data.get("conflict_type"))
        conflict_type_label = conflict_type.label
        if conflict_type == ConflictType.ONCPATH:
            if allele_origin_bucket == AlleleOriginBucket.GERMLINE:
                conflict_type_label = "Pathogenicity"
            elif allele_origin_bucket == AlleleOriginBucket.SOMATIC:
                conflict_type_label = "Oncogenicity"
            else:
                conflict_type_label = "Onco-Path"

        return {
            "allele_origin_bucket": allele_origin_bucket.value,
            "allele_origin_bucket_label": allele_origin_bucket.label,
            "testing_context_bucket": testing_context_bucket.value,
            "testing_context_bucket_label": testing_context_bucket.label,
            "conflict_type": conflict_type.value,
            "conflict_type_label": conflict_type_label
        }

    def render_c_hgvs(self, row_data: CellData):
        allele = self.allele_map.get(row_data.get("allele"))
        c_hgvses = [CHGVS.from_json_short(hgvs) for hgvs in row_data.get("meta_data").get("c_hgvs")]
        hgvs_list = []
        for c_hgvs in c_hgvses:
            c_hgvs_json = c_hgvs.to_json()
            c_hgvs_json["allele_id"] = allele.id
            c_hgvs_json["always_show_genome_build"] = True
            hgvs_list.append(c_hgvs_json)

        return hgvs_list

        # return {"text": f"A{allele.pk} {allele:CA}", "url": allele.get_absolute_url()}

    def render_severity(self, row_data: CellData):
        if row_data.value:
            cs = ConflictSeverity(row_data.value)
            return {"code": cs.value, "label": cs.label, "conflict_id": row_data.get("pk")}
        return "-"

    def involved_labs(self, row_data: CellData):
        # if data := row_data.get("data"):
        #     if involved_labs := data.get("involved_labs"):
        #         labs = list(sorted(self.lab_map[lab_id] for lab_id in involved_labs))
        #         return ", ".join(str(lab) for lab in labs)
        # return "-"
        # TODO sort
        conflict_rows = [ConflictDataRow.from_json(row) for row in row_data.get("data").get("rows")]
        conflict_labs = list(ConflictLab.objects.select_related("lab").filter(conflict=row_data.get("pk")))
        conflict_labs_dict = {cl.lab: cl for cl in conflict_labs}
        combined_rows = []
        for conflict_row in conflict_rows:
            conflict_lab: Optional[ConflictLab] = None
            if found_conflict_lab := conflict_labs_dict.get(conflict_row.lab):
                conflict_labs_dict.pop(conflict_row.lab)
                conflict_lab = found_conflict_lab

            combined_rows.append(ConflictRowWithStatus(
                conflict_row,
                conflict_lab
            ))

        return render_to_string('classification/conflict_summary_cell.html', {
            "conflict_id": row_data.get("pk"),
            "conflict_type": ConflictType(row_data.get("conflict_type")),
            "combined_rows": combined_rows
        }).strip()

    def get_initial_queryset(self) -> QuerySet[Conflict]:
        lab_ids = LabPickerData.for_user(self.user, self.get_query_param("lab")).lab_ids
        # TODO optimise so if admin with no selection we don't filter like this
        return Conflict.objects.all().for_labs(lab_ids)

    def filter_queryset(self, qs: QuerySet[Conflict]) -> QuerySet[Conflict]:
        if self.get_query_param("exclude_unknown") == "true":
            qs = qs.exclude(allele_origin_bucket=AlleleOriginBucket.UNKNOWN).exclude(testing_context_bucket=TestingContextBucket.UNKNOWN)
        if self.get_query_param("multiple_submitters") == "true":
            qs = qs.filter(severity__gt=ConflictSeverity.SINGLE_SUBMISSION)
        return qs
