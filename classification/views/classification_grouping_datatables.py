import operator
from functools import cached_property, reduce
from typing import Any, Dict, List, Optional
from django.conf import settings
from django.db.models import QuerySet, Q
from django.http import HttpRequest
from more_itertools import first
from classification.enums import ShareLevel, AlleleOriginBucket, EvidenceCategory, SpecialEKeys
from classification.models import ClassificationGrouping, ImportedAlleleInfo, ClassificationGroupingSearchTerm, \
    ClassificationGroupingSearchTermType, EvidenceKeyMap, ClassificationModification, ClassificationGroupingEntry, \
    Classification
from genes.hgvs import CHGVS
from genes.models import GeneSymbol, TranscriptVersion
from library.utils import JsonDataType
from ontology.models import OntologyTerm, OntologyTermRelation, OntologySnake
from snpdb.models import UserSettings, GenomeBuild, Lab, Variant
from snpdb.views.datatable_view import DatatableConfig, RichColumn, DC, SortOrder, CellData


class ClassificationGroupingColumns(DatatableConfig[ClassificationGrouping]):

    def render_row_header(self, row: Dict[str, Any]) -> JsonDataType:

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

    def render_lastest_curation_date(self, row: Dict[str, Any]) -> JsonDataType:
        return {
            "curation_date": row["latest_curation_date"],
            "classification_id": row["latest_classification_modification__classification_id"]
        }

    def render_somatic(self, row: CellData) -> JsonDataType:
        if row["allele_origin_bucket"] != "G":
            if somatic_dict := row["latest_classification_modification__classification__summary__somatic"]:
                return somatic_dict

    def _render_date(self, row: CellData) -> JsonDataType:
        # TODO list date type
        return {
            "classification_id": row["latest_classification_modification__classification_id"],
            "curation_date": row.get_nested_json("latest_classification_modification__classification__summary__date", "value"),
            "date_type": row.get_nested_json("latest_classification_modification__classification__summary__date", "type")
        }

    @cached_property
    def genome_build_prefs(self) -> List[GenomeBuild]:
        user_settings = UserSettings.get_for_user(self.user)
        return GenomeBuild.builds_with_annotation_priority(user_settings.default_genome_build)

    def render_c_hgvs(self, row: Dict[str, Any]) -> JsonDataType:
        def get_preferred_chgvs_json() -> Dict:
            nonlocal row
            for index, genome_build in enumerate(self.genome_build_prefs):
                if c_hgvs_string := row.get(ImportedAlleleInfo.column_name_for_build(genome_build, "latest_allele_info")):
                    c_hgvs = CHGVS(c_hgvs_string)
                    c_hgvs.genome_build = genome_build
                    c_hgvs.is_desired_build = index == 0
                    return c_hgvs.to_json()

            # FIXME need to fallback onto a plain text condition with genome build
            return {}

            # TODO - is there a need for the imported c.HGVS? (Yes but in edge cases)
            # c_hgvs = CHGVS(row.get('published_evidence__c_hgvs__value'))
            # c_hgvs.is_normalised = False
            # json_data = c_hgvs.to_json()
            # # use this rather than genome build object so we can get a patch version
            # json_data['genome_build'] = row.get('published_evidence__genome_build__value')
            # return json_data

        response = get_preferred_chgvs_json()
        # if settings.CLASSIFICATION_GRID_SHOW_PHGVS:
        #     if p_hgvs := row.get('published_evidence__p_hgvs__value'):
        #         p_dot = p_hgvs.find('p.')
        #         if p_dot != -1:
        #             p_hgvs = p_hgvs[p_dot::]
        #     response['p_hgvs'] = p_hgvs

        response['allele_id'] = row.get('latest_allele_info__allele_id')
        response['allele_info_id'] = row.get('latest_allele_info__allele_info__id')
        if warning_icon := ImportedAlleleInfo.icon_for(
            status=row.get('latest_allele_infoo__status'),
            include=row.get('latest_allele_info__latest_validation__include')
        ):
            response.update(warning_icon.as_json())

        return response

    def get_initial_queryset(self) -> QuerySet[DC]:
        # TODO, consider making the groups GuardianPermission rather than this manual security check
        base_qs = ClassificationGrouping.objects.all()
        if not self.user.is_superuser:
            permission_q: list[Q] = []
            # super user can see everyone
            labs = Lab.valid_labs_qs(self.user, admin_check=False)
            orgs = {lab.organization for lab in labs}
            permission_q.append(Q(share_level=ShareLevel.LAB) & Q(lab__in=labs))
            permission_q.append(Q(share_level=ShareLevel.INSTITUTION) & Q(lab__organization__in=orgs))
            permission_q.append(Q(share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS))
            base_qs = base_qs.filter(reduce(operator.or_, permission_q))
        return base_qs

    def filter_queryset(self, qs: QuerySet[ClassificationGrouping]) -> QuerySet[ClassificationGrouping]:
        filters: List[Q] = []
        if lab_id := self.get_query_param('lab'):
            lab_list = lab_id.split(",")
            filters.append(Q(lab_id__in=lab_list))

        if allele_id := self.get_query_param('allele_id'):
            filters.append(Q(allele_origin_grouping__allele_grouping__allele_id=int(allele_id)))

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

        if condition := self.get_query_param('ontology_term_id'):
            if c_filter := self.condition_filter(condition):
                filters.append(c_filter)

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
    def id_columns(self) -> List[str]:
        keys = EvidenceKeyMap.instance()
        return [e_key.key for e_key in keys.all_keys if '_id' in e_key.key and
                e_key.evidence_category in (EvidenceCategory.HEADER_PATIENT, EvidenceCategory.HEADER_TEST, EvidenceCategory.SIGN_OFF)]

    def classification_modification_filter_to_grouping(self, cm_q: Q) -> Q:
        return Q(
            pk__in=\
                ClassificationGroupingEntry.objects.filter(
                    classification__in=ClassificationModification.objects.filter(is_last_published=True).filter(cm_q).values_list('classification_id', flat=True)
                ).values_list('grouping_id', flat=True)
        )

    def classification_filter_to_grouping(self, cm_q: Q) -> Q:
        return Q(
            pk__in=\
                ClassificationGroupingEntry.objects.filter(
                    classification__in=Classification.objects.filter(cm_q).values_list('pk', flat=True)
                ).values_list('grouping_id', flat=True)
        )

    def id_filter(self, text: str):
        ids_contain_q_list: List[Q] = []
        for id_key in self.id_columns:
            ids_contain_q_list.append(Q(**{f'published_evidence__{id_key}__value__icontains': text}))
        id_filter_q = reduce(operator.or_, ids_contain_q_list)
        return Q(
            pk__in=\
                ClassificationGroupingEntry.objects.filter(
                    classification__in=ClassificationModification.objects.filter(is_last_published=True).filter(id_filter_q).values_list('classification_id', flat=True)
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

            all_terms = set([term])
            if mondo_term := OntologyTermRelation.as_mondo(term):
                # if we have (or can translate) into a mondo term
                # get all the direct parent and all the direct children terms
                # and then the OMIM equiv of them
                mondo_terms = set([mondo_term])
                mondo_terms |= OntologySnake.get_children(mondo_term)
                mondo_terms |= OntologySnake.get_parents(mondo_term)
                all_terms |= mondo_terms
                for term in mondo_terms:
                    if omim_term := OntologyTermRelation.as_omim(term):
                        all_terms.add(omim_term)

            all_strs = [term.id.upper() for term in all_terms]
            return ClassificationGroupingSearchTerm.filter_q(ClassificationGroupingSearchTermType.CONDITION_ID, all_strs)
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
                sort_keys=['lab__organization__name', 'lab__name'],
                name='id',
                label='ID',
                orderable=True,
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
                label=f'HGVS ({genome_build_preferred.name})',
                renderer=self.render_c_hgvs,
                client_renderer='VCTable.hgvs',
                orderable=True,
                extra_columns=[
                    "latest_allele_info__grch37__c_hgvs",
                    "latest_allele_info__grch38__c_hgvs",
                    'latest_allele_info__id',
                    'latest_allele_info__allele_id',
                    'latest_allele_info__latest_validation__include',
                    'latest_allele_info__status'
                ]
            ),
            RichColumn(
                key="latest_classification_modification__classification__summary__classification_value",
                name="Classification",
                sort_keys=["latest_classification_modification__classification__summary__classification_sort"],
                client_renderer='VCTable.classification',
                order_sequence=[SortOrder.DESC, SortOrder.ASC]
            ),
            RichColumn(
                name='somatic_clinical_significances',
                label='Somatic Clinical<br/>Significance',
                client_renderer="VCTable.somatic_clinical_significance",
                sort_keys=['latest_classification_modification__classification__summary__pathogenicity__sort'],
                order_sequence=[SortOrder.DESC, SortOrder.ASC],
                renderer=self.render_somatic,
                extra_columns=[
                    "latest_classification_modification__classification__summary__somatic",
                    "allele_origin_bucket"
                ]
            ),
            RichColumn(
                key="latest_classification_modification__classification__summary__criteria_labels",
                name="latest_criteria",
                label="Latest Criteria",
                client_renderer='TableFormat.list_codes'
            ),
            RichColumn(
                key='conditions',
                name='conditions',
                label='Conditions',
                sort_keys=['conditions__sort_text'],
                client_renderer='VCTable.condition',
                orderable=True
            ),
            RichColumn(
                key="latest_classification_modification__classification__summary__date",
                name="latest_curation_date",
                label="Latest Curated",
                sort_keys=["latest_classification_modification__classification__summary__date__value"],
                client_renderer="VCTable.latest_curation_and_link",
                renderer=self._render_date,
                extra_columns=[
                    "latest_classification_modification__classification_id"
                ],
                default_sort=SortOrder.DESC
            ),
            # RichColumn(
            #     key='lab',
            #     # share_level_sort annotated column
            #     sort_keys=['lab__organization__name', 'lab__name'],
            #     name='id',
            #     label='ID',
            #     orderable=True,
            #     renderer=self.render_row_header,
            #     client_renderer='VCTable.groupIdentifier',
            #     extra_columns=[
            #         'id',
            #         'classification_count',
            #         'lab__organization__short_name',
            #         'lab__organization__name',
            #         'lab__name',
            #         'allele_origin_bucket',
            #         'share_level',
            #         'dirty'
            #     ]
            # ),
            # RichColumn(
            #     key=ImportedAlleleInfo.column_name_for_build(genome_build_preferred, "latest_allele_info"),
            #     # sort_keys=['variant_sort', 'c_hgvs'],  # annotated column
            #     sort_keys=[ImportedAlleleInfo.column_name_for_build(genome_build_preferred, "latest_allele_info", "genomic_sort")],
            #     name='c_hgvs',
            #     label=f'HGVS ({genome_build_preferred.name})',
            #     renderer=self.render_c_hgvs,
            #     client_renderer='VCTable.hgvs',
            #     orderable=True,
            #     extra_columns=[
            #         "latest_allele_info__grch37__c_hgvs",
            #         "latest_allele_info__grch38__c_hgvs",
            #         #'published_evidence__c_hgvs__value',
            #         #'published_evidence__p_hgvs__value',
            #         #'published_evidence__genome_build__value',
            #         'latest_allele_info__id',
            #         'latest_allele_info__allele_id',
            #         'latest_allele_info__latest_validation__include',
            #         'latest_allele_info__status'
            #     ]
            # ),
            #
            # RichColumn(
            #     key='classification_values',
            #     name='classifications',
            #     label='Classifications',
            #     client_renderer=RichColumn.client_renderer_repeat({"formatter": 'VCTable.classification'}),
            #     sort_keys=['classification_sort_value'], # FIXME add a sort column
            #     order_sequence=[SortOrder.DESC, SortOrder.ASC]
            # ),
            #
            # RichColumn(
            #     key='somatic_clinical_significance_values',
            #     name='somatic_clinical_significances',
            #     label='Somatic Clinical<br/>Significance',
            #     client_renderer=RichColumn.client_renderer_repeat({"formatter": 'VCTable.somatic_clinical_significance'}),
            #     sort_keys=['somatic_clinical_significance_sort'],
            #     order_sequence=[SortOrder.DESC, SortOrder.ASC]
            # ),
            #
            # RichColumn(
            #     key='conditions',
            #     name='conditions',
            #     label='Conditions',
            #     sort_keys=['conditions__sort_text'],
            #     client_renderer='VCTable.condition',
            #     orderable=True
            # ),
            #
            # RichColumn(
            #     key='latest_criteria',
            #     name='latest_criteria',
            #     label='Latest Criteria',
            #     client_renderer='TableFormat.list_codes',
            #     sort_keys=['latest_criteria'],
            #     orderable=True
            # ),
            #
            # RichColumn(
            #     key='latest_classification_modification__classification__summary__date_value',
            #     name='latest_classification',
            #     label="Last Curated",
            #     renderer=self.render_lastest_curation_date,
            #     client_renderer='VCTable.latest_curation_and_link',
            #     extra_columns=[
            #         "latest_classification_modification__classification_id",
            #         "latest_classification_modification__classification__summary"
            #     ],
            #     sort_keys=["latest_classification_modification__classification__summary__"],
            #     order_sequence=[SortOrder.DESC, SortOrder.ASC]
            # )
        ]