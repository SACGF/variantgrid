import operator
from functools import cached_property, reduce
from typing import Any, Dict, List, Optional
from django.conf import settings
from django.db.models import QuerySet, Q
from django.http import HttpRequest
from more_itertools import first
from classification.enums import AlleleOriginBucket, EvidenceCategory, SpecialEKeys
from classification.models import ClassificationGrouping, ImportedAlleleInfo, ClassificationGroupingSearchTerm, \
    ClassificationGroupingSearchTermType, EvidenceKeyMap, ClassificationModification, ClassificationGroupingEntry, \
    Classification, DiscordanceReport, DiscordanceReportClassification
from genes.hgvs import CHGVS
from genes.models import GeneSymbol, TranscriptVersion
from library.utils import JsonDataType
from ontology.models import OntologyTerm, OntologyTermRelation, OntologySnake
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import UserSettings, GenomeBuild, Variant
from snpdb.views.datatable_view import DatatableConfig, RichColumn, DC, SortOrder, CellData


class ClassificationGroupingColumns(DatatableConfig[ClassificationGrouping]):
    """
    File to display the ClassificationGrouping as data tables.
    This is taking over from ClassificationGroup (in memory merging of rows that then have to be all rendered client side)
    ClassificationColumns - the one classification per row

    Filters are done between a combination of get_initial_query (when on the gene symbol page and we want the maximum
    number of results to be how many groupings for that gene symbol, rather than presenting the user with
    showing 6 out of 30,000 records)
    and then on filter_query_set when it's a filter that the user can do above and beyond what the page initially loads
    e.g. allele origin filter
    """

    def render_row_header(self, row: CellData) -> JsonDataType:

        matches: Optional[Dict[str, str]] = None
        search: Optional[str] = None

        if settings.CLASSIFICATION_ID_FILTER:
            if id_filter := self.get_query_param("id_filter"):
                search = id_filter
                id_filter = id_filter.lower()
                matches = {}

                for cm in ClassificationGrouping.objects.get(pk=row.get("id")).classification_modifications:
                    for id_key in self.id_columns:
                        if (value := cm.get(id_key)) and id_filter in value.lower():
                            matches[id_key] = value

        return {
            "dirty": row.get("dirty"),
            "id": row.get("id"),
            "classification_count": row.get('classification_count'),
            "org_name": row.get('lab__organization__short_name') or row.get('lab__organization__name'),
            "lab_name": row.get('lab__name'),
            "share_level": row.get('share_level'),
            "allele_origin_bucket": row.get('allele_origin_bucket'),
            "matches": matches,
            "search": search
        }

    # def render_latest_curation_date(self, row: CellData) -> JsonDataType:
    #     return {
    #         "curation_date": row["latest_curation_date"],
    #         "classification_id": row["latest_classification_modification__classification_id"]
    #     }

    def render_somatic(self, row: CellData) -> JsonDataType:
        if row["allele_origin_bucket"] != "G":
            diff_value = row["somatic_difference"]
            if somatic_dict := row["latest_classification_modification__classification__summary__somatic"]:
                somatic_dict["diff"] = diff_value
                return somatic_dict

    def render_pathogenic(self, row: CellData) -> JsonDataType:
        diff_value = row["pathogenic_difference"]
        result_dict = row["latest_classification_modification__classification__summary__pathogenicity"] or {}
        result_dict["diff"] = diff_value

        if dr := self.discordance_report:
            if drc := DiscordanceReportClassification.objects.filter(report_id=dr.pk,
                                                                     classification_original__classification=row[
                                                                         "latest_classification_modification__classification_id"]).first():
                old_cs = drc.classification_original.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
                if result_dict and result_dict.get("classification") != old_cs:
                    result_dict["old"] = old_cs
                if "pending" not in result_dict:
                    effective_cs = drc.classification_effective.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
                    if effective_cs != result_dict.get("classification"):
                        result_dict["new"] = result_dict.get("classification")
                        result_dict["classification"] = effective_cs

        return result_dict

    def _render_date(self, row: CellData) -> JsonDataType:
        return {
            "classification_id": row["latest_classification_modification__classification_id"],
            "curation_date": row.get_nested_json("latest_classification_modification__classification__summary__date",
                                                 "value"),
            "date_type": row.get_nested_json("latest_classification_modification__classification__summary__date",
                                             "type")
        }

    @cached_property
    def genome_build_prefs(self) -> List[GenomeBuild]:
        return GenomeBuild.builds_with_annotation_priority(GenomeBuildManager.get_current_genome_build())

    def render_c_hgvs(self, row: CellData) -> JsonDataType:
        def get_preferred_chgvs_json() -> Dict:
            nonlocal row
            for index, genome_build in enumerate(self.genome_build_prefs):
                if c_hgvs_string := row.get(ImportedAlleleInfo.column_name_for_build(genome_build, "latest_allele_info")):
                    c_hgvs = CHGVS(c_hgvs_string)
                    c_hgvs.genome_build = genome_build
                    c_hgvs.is_desired_build = index == 0
                    return c_hgvs.to_json()

            # May still have linked to an allele without having the c_hgvs on either build
            # TODO check imported g_hgvs or other importable columns
            c_hgvs = CHGVS(row["latest_allele_info__imported_c_hgvs"])
            c_hgvs.genome_build = GenomeBuild.get_name_or_alias(row["latest_allele_info__imported_genome_build_patch_version__genome_build"])
            c_hgvs.is_normalised = False
            return c_hgvs.to_json()

        response = get_preferred_chgvs_json()
        if settings.CLASSIFICATION_GRID_SHOW_PHGVS:
            if p_hgvs := row['latest_classification_modification__published_evidence__p_hgvs__value']:
                p_dot = p_hgvs.find('p.')
                if p_dot != -1:
                    p_hgvs = p_hgvs[p_dot::]
            response['p_hgvs'] = p_hgvs

        response['allele_id'] = row.get('latest_allele_info__allele_id')
        response['allele_info_id'] = row.get('latest_allele_info__allele_info__id')
        if warning_icon := ImportedAlleleInfo.icon_for(
                status=row.get('latest_allele_info__status'),
                include=row.get('latest_allele_info__latest_validation__include')
        ):
            response.update(warning_icon.as_json())

        return response

    def get_initial_queryset(self) -> QuerySet[DC]:
        qs = ClassificationGrouping.filter_for_user(self.user, ClassificationGrouping.objects.all())

        page = self.get_query_param('page_id')

        filters: List[Q] = []

        # run the filters that are perma-applied on certain pages

        if allele_id := self.get_query_param('allele_id'):
            filters.append(Q(allele_origin_grouping__allele_grouping__allele_id=int(allele_id)))

        if condition := self.get_query_param('ontology_term_id'):
            if c_filter := self.condition_filter(condition):
                filters.append(c_filter)

        if page == "gene_symbol":
            if gene_symbol_str := self.get_query_param("gene_symbol"):
                if gs_filter := self.gene_symbol_filter(gene_symbol_str):
                    filters.append(gs_filter)

        if dr := self.discordance_report:
            classification_ids = [cm.classification_id for cm in dr.all_classification_modifications]
            group_ids = ClassificationGroupingEntry.objects.filter(
                classification_id__in=classification_ids).values_list('grouping', flat=True)
            filters.append(Q(pk__in=group_ids))

        if filters:
            return qs.filter(*filters)
        else:
            return qs

    def filter_queryset(self, qs: QuerySet[ClassificationGrouping]) -> QuerySet[ClassificationGrouping]:
        page = self.get_query_param('page_id')

        # run the filters that are optionally applied

        filters: List[Q] = []
        if lab_id := self.get_query_param('lab'):
            lab_list = lab_id.split(",")
            filters.append(Q(lab_id__in=lab_list))

        if allele_origin := self.get_query_param("allele_origin"):
            if allele_origin != "A":
                filters.append(Q(allele_origin_bucket__in=[allele_origin, AlleleOriginBucket.UNKNOWN]))

        # for view gene symbol
        if protein_position := self.get_query_param("protein_position"):
            protein_position_transcript_version_id = self.get_query_param("protein_position_transcript_version_id")
            transcript_version = TranscriptVersion.objects.get(pk=protein_position_transcript_version_id)
            variant_qs = Variant.objects.filter(varianttranscriptannotation__transcript_version=transcript_version,
                                                varianttranscriptannotation__protein_position__icontains=protein_position)
            # Join through allele so it works across genome builds
            filters.append(Q(allele_origin_grouping__allele_grouing__allele__variantallele__variant__in=variant_qs))

        if page != "gene_symbol":
            if gene_symbol_str := self.get_query_param("gene_symbol"):
                if gs_filter := self.gene_symbol_filter(gene_symbol_str):
                    filters.append(gs_filter)

        if settings.CLASSIFICATION_ID_FILTER:
            if filter_text := self.get_query_param("id_filter"):
                filters.append(self.id_filter(filter_text))

        if settings.CLASSIFICATION_GRID_SHOW_USERNAME:
            if user_id := self.get_query_param('user'):
                filters.append(self.classification_filter_to_grouping(Q(user__pk=user_id)))

        return qs.filter(*filters)

    @cached_property
    def discordance_report(self) -> Optional[DiscordanceReport]:
        if discordance_report_id := self.get_query_param('discordance_report'):
            dr = DiscordanceReport.objects.get(pk=discordance_report_id)
            dr.check_can_view(self.user)
            return dr

    @cached_property
    def id_columns(self) -> List[str]:
        keys = EvidenceKeyMap.instance()
        return [e_key.key for e_key in keys.all_keys if '_id' in e_key.key and
                e_key.evidence_category in (
                    EvidenceCategory.HEADER_PATIENT, EvidenceCategory.HEADER_TEST, EvidenceCategory.SIGN_OFF)]

    def classification_modification_filter_to_grouping(self, cm_q: Q) -> Q:
        return Q(
            pk__in=ClassificationGroupingEntry.objects.filter(
                    classification__in=ClassificationModification.objects.filter(is_last_published=True).filter(
                        cm_q).values_list('classification_id', flat=True)
                ).values_list('grouping_id', flat=True)
        )

    def classification_filter_to_grouping(self, cm_q: Q) -> Q:
        return Q(
            pk__in=ClassificationGroupingEntry.objects.filter(
                    classification__in=Classification.objects.filter(cm_q).values_list('pk', flat=True)
                ).values_list('grouping_id', flat=True)
        )

    def id_filter(self, text: str):
        ids_contain_q_list: List[Q] = []
        for id_key in self.id_columns:
            ids_contain_q_list.append(Q(**{f'published_evidence__{id_key}__value__icontains': text}))
        id_filter_q = reduce(operator.or_, ids_contain_q_list)
        return Q(
            pk__in=ClassificationGroupingEntry.objects.filter(
                    classification__in=ClassificationModification.objects.filter(is_last_published=True).filter(
                        id_filter_q).values_list('classification_id', flat=True)
                ).values_list('grouping_id', flat=True)
        )

    def gene_symbol_filter(self, gene_symbol: str):
        if gene_symbol := GeneSymbol.objects.filter(symbol=gene_symbol).first():
            all_strs = [gene_symbol.symbol] + gene_symbol.alias_meta.alias_symbol_strs
            all_strs = [gs.upper() for gs in all_strs]
            return ClassificationGroupingSearchTerm.filter_q(ClassificationGroupingSearchTermType.GENE_SYMBOL, all_strs)
        # FIXME add support for gene symbol alias

    def scv_filter(self, scv: str):
        if scv.startswith("SCV"):
            return ClassificationGroupingSearchTerm.filter_q(ClassificationGroupingSearchTermType.CLINVAR_SCV, scv)

    def condition_filter(self, text, must_exist: bool = False):
        try:
            term = OntologyTerm.get_or_stub(text)
            if must_exist and term.is_stub:
                return False

            all_terms = {term}
            if mondo_term := OntologyTermRelation.as_mondo(term):
                # if we have (or can translate) into a mondo term
                # get all the direct parent and all the direct children terms
                # and then the OMIM equiv of them
                mondo_terms = {mondo_term}
                mondo_terms |= OntologySnake.get_children(mondo_term)
                mondo_terms |= OntologySnake.get_parents(mondo_term)
                all_terms |= mondo_terms
                for term in mondo_terms:
                    if omim_term := OntologyTermRelation.as_omim(term):
                        all_terms.add(omim_term)

            all_strs = [term.id.upper() for term in all_terms]
            return ClassificationGroupingSearchTerm.filter_q(ClassificationGroupingSearchTermType.CONDITION_ID,
                                                             all_strs)
        except ValueError:
            pass

    def __init__(self, request: HttpRequest):
        super().__init__(request)

        genome_build_preferred = first(self.genome_build_prefs)

        self.expand_client_renderer = DatatableConfig._row_expand_ajax('classification_grouping_detail',
                                                                       expected_height=108)

        self.rich_columns = [
            RichColumn(
                key='lab',
                # share_level_sort annotated column
                sort_keys=[
                    'allele_origin_bucket',
                    'lab__organization__name',
                    'lab__name',
                    'share_level'
                ],
                name='id',
                label='Lab',
                order_sequence=[SortOrder.ASC, SortOrder.DESC],
                renderer=self.render_row_header,
                client_renderer='VCTable.groupIdentifier',
                extra_columns=[
                    'id',
                    'classification_count',
                    'lab__organization__short_name',
                    'lab__organization__name',
                    'lab__name',
                    'allele_origin_bucket',
                    'share_level',
                    'dirty'
                ]
            ),
            RichColumn(
                key=ImportedAlleleInfo.column_name_for_build(genome_build_preferred, "latest_allele_info"),
                # sort_keys=['variant_sort', 'c_hgvs'],  # annotated column
                sort_keys=[ImportedAlleleInfo.column_name_for_build(genome_build_preferred, "latest_allele_info",
                                                                    "genomic_sort")],
                name='c_hgvs',
                label=f'HGVS <span class="text-secondary">{genome_build_preferred.name}</span>',
                renderer=self.render_c_hgvs,
                client_renderer='VCTable.hgvs',
                order_sequence=[SortOrder.ASC, SortOrder.DESC],
                extra_columns=[
                    "latest_allele_info__grch37__c_hgvs",
                    "latest_allele_info__grch38__c_hgvs",
                    'latest_allele_info__id',
                    'latest_allele_info__allele_id',
                    'latest_allele_info__latest_validation__include',
                    'latest_allele_info__status',
                    'latest_allele_info__imported_c_hgvs',
                    'latest_allele_info__imported_genome_build_patch_version__genome_build',
                    'latest_classification_modification__published_evidence__p_hgvs__value'  # TODO move this to allele info
                ]
            ),
            RichColumn(
                name="Classification",
                sort_keys=[
                    "latest_classification_modification__classification__summary__pathogenicity__sort",
                    "latest_classification_modification__classification__summary__somatic__sort"
                ],
                client_renderer='VCTable.classification',
                renderer=self.render_pathogenic,
                order_sequence=[SortOrder.DESC, SortOrder.ASC],
                extra_columns=[
                    "latest_classification_modification__classification_id",
                    "latest_classification_modification__classification__summary__pathogenicity",
                    "pathogenic_difference"
                ]
            ),
            RichColumn(
                name='somatic_clinical_significances',
                label='Somatic Clinical<br/>Significance',
                client_renderer="VCTable.somatic_clinical_significance",
                sort_keys=[
                    'latest_classification_modification__classification__summary__somatic__sort',
                    'latest_classification_modification__classification__summary__pathogenicity__sort'
                ],
                order_sequence=[SortOrder.DESC, SortOrder.ASC],
                renderer=self.render_somatic,
                extra_columns=[
                    "latest_classification_modification__classification__summary__somatic",
                    "allele_origin_bucket",
                    "somatic_difference"
                ]
            ),
            RichColumn(
                key="latest_classification_modification__classification__summary__criteria_labels",
                name="latest_criteria",
                label='<span class="text-secondary">Latest</span><br/>Criteria',
                client_renderer='TableFormat.list_codes'
            ),
            RichColumn(
                key='conditions',
                name='conditions',
                label='Conditions',
                sort_keys=['conditions__sort_text'],
                client_renderer='VCTable.condition',
                order_sequence=[SortOrder.ASC, SortOrder.DESC],
            ),
            RichColumn(
                key="latest_classification_modification__classification__summary__date",
                name="latest_curation_date",
                label='<span class="text-secondary">Latest</span><br/>Curated',
                sort_keys=["latest_classification_modification__classification__summary__date__value"],
                client_renderer="VCTable.latest_curation_and_link",
                renderer=self._render_date,
                extra_columns=[
                    "latest_classification_modification__classification_id"
                ],
                order_sequence=[SortOrder.DESC, SortOrder.ASC],
                default_sort=SortOrder.DESC
            )
        ]
