import operator
from functools import reduce
from typing import Dict, Any, List, Optional
from django.conf import settings
from django.db.models import Q, Subquery, When, Case, TextField, Value, IntegerField, QuerySet
from django.db.models.fields.json import KeyTextTransform, KeyTransform
from django.db.models.functions import Lower, Cast
from django.http import HttpRequest
from lazy import lazy
from classification.models.classification_utils import classification_gene_symbol_filter
from flags.models import FlagCollection, FlagStatus
from genes.models import TranscriptVersion
from ontology.models import OntologyTerm
from snpdb.models import UserSettings, GenomeBuild, Variant, Lab
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder
from classification.enums import SpecialEKeys, EvidenceCategory, ShareLevel
from classification.models import ClassificationModification, classification_flag_types, Classification, EvidenceKeyMap
from uicore.json.json_types import JsonDataType

ALLELE_GERMLINE_VALUES = ['germline', 'likely_germline']
ALLELE_SOMATIC_VALUES = ['somatic', 'likely_somatic']
ALLELE_KNOWN_VALUES = ALLELE_GERMLINE_VALUES + ALLELE_SOMATIC_VALUES


class ClassificationColumns(DatatableConfig[ClassificationModification]):

    def render_c_hgvs(self, row: Dict[str, Any]) -> JsonDataType:
        values = []
        for index, genome_build in enumerate(self.genome_build_prefs):
            values.append({
                'type': 'normal-pref' if index == 0 else 'normal',
                'build': genome_build.pk,
                'value': row.get(ClassificationModification.column_name_for_build(genome_build))
            })
        values.append({
            'type': 'imported',
            'build': row.get('published_evidence__genome_build__value'),
            'value': row.get('published_evidence__c_hgvs__value'),
        })
        response = {
            'values': values,
            'variant_id': row.get('classification__variant_id')
        }

        if settings.VARIANT_CLASSIFICATION_GRID_SHOW_PHGVS:
            if p_hgvs := row.get('published_evidence__p_hgvs__value'):
                p_dot = p_hgvs.find('p.')
                if p_dot != -1:
                    p_hgvs = p_hgvs[p_dot::]
            response['p_hgvs'] = p_hgvs

        return response

    def render_condition(self, row: Dict[str, Any]) -> JsonDataType:
        if cr := row['classification__condition_resolution']:
            return cr
        else:
            return {"display_text": row['published_evidence__condition__value']}

    def classification_id(self, row: Dict[str, Any]) -> JsonDataType:
        matches: Optional[Dict[str, str]] = None
        if id_filter := self.get_query_param("id_filter"):
            matches = dict()
            id_keys = self.id_columns
            for key in id_keys:
                value = row.get(f'published_evidence__{key}__value')
                if value:
                    value = str(value)
                    if id_filter.lower() in value.lower():
                        matches[key] = value

        return {
            "id": row.get('classification__id'),
            "lab_name": row.get('classification__lab__name'),
            "lab_record_id": row.get('classification__lab_record_id'),
            "share_level": row.get('classification__share_level'),
            "matches": matches,
            "search": id_filter
        }

    @lazy
    def genome_build_prefs(self) -> List[GenomeBuild]:
        user_settings = UserSettings.get_for_user(self.user)
        return GenomeBuild.builds_with_annotation_priority(user_settings.default_genome_build)

    def __init__(self, request: HttpRequest):
        self.term_cache: Dict[str, OntologyTerm] = dict()
        super().__init__(request)

        user_settings = UserSettings.get_for_user(self.user)

        self.rich_columns = [
            RichColumn(
                key='classification__id',
                # share_level_sort annotated column
                sort_keys=['share_level_sort', 'classification__lab__name', 'classification__lab_record_id'],
                name='id',
                label='ID',
                orderable=True,
                renderer=self.classification_id,
                client_renderer='VCTable.identifier',
                extra_columns=[
                    'classification__id',
                    'classification__lab__name',
                    'classification__lab_record_id',
                    'classification__share_level'
                ]
            ),
            RichColumn(
                key='published_evidence__gene_symbol__value',
                name='gene_symbol',
                label='Gene Symbol',
                client_renderer=f'VCTable.evidence_key.bind(null, "{ SpecialEKeys.GENE_SYMBOL }")',
                orderable=True
            ),
            RichColumn(
                key=ClassificationModification.column_name_for_build(self.genome_build_prefs[0]),
                sort_keys=['c_hgvs'],  # annotated column
                name='c_hgvs',
                label=f'HGVS ({user_settings.default_genome_build.name})',
                renderer=self.render_c_hgvs,
                client_renderer='VCTable.c_hgvs',
                orderable=True,
                extra_columns=[
                    "classification__chgvs_grch37",
                    "classification__chgvs_grch38",
                    'published_evidence__c_hgvs__value',
                    'published_evidence__p_hgvs__value',
                    'published_evidence__genome_build__value',
                    'classification__variant_id',
                ]
            ),
            RichColumn(
                key='published_evidence__clinical_significance__value',
                name='clinical_significance',
                label='Clinical Significance',
                client_renderer='VCTable.clinical_significance',
                client_renderer_td='VCTable.clinical_significance_td',
                sort_keys=['clinical_significance', 'clin_sig_sort'],
                orderable=True
            ),
            RichColumn(
                key='published_evidence__allele_origin__value',
                name='allele_origin',
                label='Allele Origin',
                client_renderer=f'VCTable.evidence_key.bind(null, "{ SpecialEKeys.ALLELE_ORIGIN }")',
                enabled=settings.VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN,
                orderable=True
            ),
            RichColumn(
                key='published_evidence__condition__value',
                name='condition',
                label='Condition',
                sort_keys=['condition_sort'],
                renderer=self.render_condition,
                client_renderer='VCTable.condition',
                orderable=True,
                extra_columns=["classification__condition_resolution"]
            ),
            RichColumn(
                key='classification__user__username',
                name='user',
                label='User',
                enabled=settings.VARIANT_CLASSIFICATION_GRID_SHOW_USERNAME,
                orderable=True
            ),
            RichColumn(
                key='classification__created',
                name='created',
                label='Created',
                orderable=True,
                client_renderer='TableFormat.timestamp',
                default_sort=SortOrder.DESC
            ),
            RichColumn(
                key='classification__flag_collection_id',
                name='flags',
                label='Flags',
                client_renderer='TableFormat.flags'
            )
        ]

    def get_initial_queryset(self) -> QuerySet[ClassificationModification]:

        exclude_withdrawn = True
        if flags := self.get_query_json("flags"):
            if "classification_withdrawn" in flags:
                exclude_withdrawn = False

        initial_qs = ClassificationModification.latest_for_user(
            user=self.user,
            published=True,
            exclude_withdrawn=exclude_withdrawn)

        # filtering to your labs only is done on the
        if my_labs := self.get_query_param('my_labs'):
            if my_labs == "true":
                initial_qs = initial_qs.filter(classification__lab__in=Lab.valid_labs_qs(user=self.user, admin_check=True))

        # Make an annotated column c_hgvs which is the first non null value of
        # user's normalised preference (e.g. 37), the alternative normalised (e.g. 38) the imported c.hgvs
        whens = []
        for genome_build in self.genome_build_prefs:
            column_name = ClassificationModification.column_name_for_build(genome_build)
            whens.append(When(**{f'{column_name}__isnull': False, 'then': column_name}))
        case = Case(*whens, default=KeyTextTransform('value', KeyTransform('c_hgvs', 'published_evidence')),
                    output_field=TextField())
        initial_qs = initial_qs.annotate(c_hgvs=case)

        # So VUS-A is actually VUS-Pathogenic, VUS-C is VUS-Benign
        # so when sorting within VUS we want to do so in reverse alphabetical order
        # but when sorting outside of VUS we want alphabetical order
        # Note: Clinical Significance will first be sorted by the clinical significant int and then fall back to this
        whens = [
            When(published_evidence__clinical_significance__value='VUS_A', then=Value('VUS_3')),
            When(published_evidence__clinical_significance__value='VUS_B', then=Value('VUS_2')),
            When(published_evidence__clinical_significance__value='VUS', then=Value('VUS_1')),
            When(published_evidence__clinical_significance__value='VUS_C', then=Value('VUS_0')),
        ]
        case = Case(*whens, default=KeyTextTransform('value', KeyTransform('clinical_significance', 'published_evidence')),
                    output_field=TextField())
        initial_qs = initial_qs.annotate(clin_sig_sort=case)

        whens = [
            When(classification__share_level=ShareLevel.LAB.value, then=Value(1)),
            When(classification__share_level=ShareLevel.INSTITUTION.value, then=Value(2)),
            When(classification__share_level=ShareLevel.ALL_USERS.value, then=Value(3)),
            When(classification__share_level=ShareLevel.PUBLIC.value, then=Value(4))
        ]
        case = Case(*whens, default=Value(0), output_field=IntegerField())
        initial_qs = initial_qs.annotate(share_level_sort=case)

        case = Case(
            When(
                classification__condition_resolution__sort_text__isnull=False,
                then=Cast(KeyTextTransform('sort_text', 'classification__condition_resolution'), TextField())
            ),
            default=Lower(Cast(KeyTextTransform('value', KeyTransform('condition', 'published_evidence')), TextField())),
            output_field=TextField()
        )
        initial_qs = initial_qs.annotate(condition_sort=case)

        return initial_qs

    @lazy
    def id_columns(self) -> List[str]:
        keys = EvidenceKeyMap.instance()
        return [e_key.key for e_key in keys.all_keys if '_id' in e_key.key and
                e_key.evidence_category in (EvidenceCategory.HEADER_PATIENT, EvidenceCategory.HEADER_TEST, EvidenceCategory.SIGN_OFF)]

    def value_columns(self) -> List[str]:
        all_columns = super().value_columns()
        if self.get_query_param("id_filter"):
            id_keys = self.id_columns
            for id_key in id_keys:
                all_columns.append(f'published_evidence__{id_key}__value')
        return all_columns

    def filter_queryset(self, qs: QuerySet[ClassificationModification]) -> QuerySet[ClassificationModification]:
        filters: List[Q] = list()
        if settings.VARIANT_CLASSIFICATION_GRID_SHOW_USERNAME:
            if user_id := self.get_query_param('user'):
                filters.append(Q(classification__user__pk=user_id))

        if lab_id := self.get_query_param('lab'):
            filters.append(Q(classification__lab__pk=lab_id))

        if flags := self.get_query_json("flags"):
            # Use inner query vs join to only return unique results
            flag_collections = FlagCollection.objects.filter(flag__flag_type__in=flags,
                                                             flag__resolution__status=FlagStatus.OPEN)
            filters.append(Q(classification__flag_collection__in=flag_collections))

        if id_filter := self.get_query_param("id_filter"):
            id_keys = self.id_columns
            ids_contain_q_list: List[Q] = list()
            for id_key in id_keys:
                ids_contain_q_list.append(Q(**{f'published_evidence__{id_key}__value__icontains': id_filter}))
            id_filter_q = reduce(operator.or_, ids_contain_q_list)
            filters.append(id_filter_q)

        if gene_symbol_str := self.get_query_param("gene_symbol"):
            if gene_filter := classification_gene_symbol_filter(gene_symbol_str):
                filters.append(gene_filter)
            else:
                return qs.none()

        if settings.VARIANT_CLASSIFICATION_GRID_SHOW_ORIGIN:
            allele_origin = self.get_query_param("allele_origin")
            if allele_origin:
                if allele_origin == 'germline':
                    filters.append(Q(published_evidence__allele_origin__value__in=ALLELE_GERMLINE_VALUES))
                elif allele_origin == 'somatic':
                    filters.append(Q(published_evidence__allele_origin__value__in=ALLELE_SOMATIC_VALUES))
                elif allele_origin == 'other':
                    filters.append(
                        ~Q(published_evidence__allele_origin__value__in=ALLELE_KNOWN_VALUES) |
                        Q(published_evidence__allele_origin__value__isnull=True)
                    )

        if issues_type := self.get_query_param("issues"):
            flag_types = None
            if issues_type == 'data':
                flag_types = [
                    classification_flag_types.matching_variant_flag,
                    classification_flag_types.matching_variant_warning_flag,
                    classification_flag_types.unshared_flag,
                    classification_flag_types.classification_outstanding_edits]
            elif issues_type == 'discussion':
                flag_types = [
                    classification_flag_types.discordant,
                    classification_flag_types.internal_review,
                    classification_flag_types.significance_change,
                    classification_flag_types.suggestion
                ]
            elif issues_type == 'withdrawn':
                flag_types = [
                    classification_flag_types.classification_withdrawn
                ]

            vcqs = FlagCollection.filter_for_open_flags(
                Classification.filter_for_user(user=self.user),
                flag_types=flag_types  # if flag types was not set it will be all flag types
            ).order_by('-created')

            filters.append(Q(
                classification__in=Subquery(vcqs.values('pk')),
            ))

        # ClassificationListGrid is also used on patient/sample page so we
        # need to filter by samples
        if sample_ids := self.get_query_json("sample_ids"):
            filters.append(Q(classification__sample__in=sample_ids))

        if protein_position := self.get_query_param("protein_position"):
            protein_position_transcript_version_id = self.get_query_param("protein_position_transcript_version_id")
            transcript_version = TranscriptVersion.objects.get(pk=protein_position_transcript_version_id)
            variant_qs = Variant.objects.filter(varianttranscriptannotation__transcript_version=transcript_version,
                                                varianttranscriptannotation__protein_position__icontains=protein_position)
            # Join through allele so it works across genome builds
            filters.append(Q(classification__variant__variantallele__allele__variantallele__variant__in=variant_qs))

        if analysis_id := self.get_query_json("analysis_id"):
            filters.append(Q(classification__analysisclassification__analysis_id=analysis_id))

        if term_id := self.get_query_param("ontology_term_id"):
            filters.append(Q(classification__condition_resolution__resolved_terms__contains=[{"term_id": term_id}]))

        if filters:
            q = reduce(operator.and_, filters)
            qs = qs.filter(q)

        return super().filter_queryset(qs)
