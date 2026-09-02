import operator
from functools import cached_property, reduce
from typing import Any, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Case, F, IntegerField, OuterRef, QuerySet, StringAgg, Subquery, Value, When
from django.db.models.aggregates import Count, Max
from django.db.models.fields import TextField
from django.db.models.query_utils import Q
from django.http import HttpRequest
from django.shortcuts import get_object_or_404
from django.urls import reverse
from guardian.shortcuts import get_objects_for_user

from annotation.annotation_version_querysets import get_queryset_for_latest_annotation_version
from annotation.models import PATIENT_ONTOLOGY_TERM_PATH, AnnotationVersion, ManualVariantEntryCollection
from annotation.models.models_enums import ClinVarReviewStatus
from library.genomics.vcf_enums import INFO_LIFTOVER_SWAPPED_REF_ALT
from library.unit_percent import get_allele_frequency_formatter
from library.utils import JsonDataType, JsonObjType, calculate_age
from ontology.models import OntologyService
from patients.models_enums import Sex
from snpdb.grid_columns.custom_columns import get_variant_grid_columns, get_variantgrid_extra_annotate
from snpdb.models import (
    VCF,
    Allele,
    AlleleConversionTool,
    AlleleLiftover,
    ClinGenAllele,
    Cohort,
    CohortGenotypeStats,
    CustomColumnsCollection,
    GenomeBuild,
    GenomicIntervalsCollection,
    ImportSource,
    ImportStatus,
    LiftoverRun,
    ProcessingStatus,
    Quad,
    Sample,
    SuperPopulationCode,
    TagColorsCollection,
    Trio,
    UserGridConfig,
    UserSettings,
    Variant,
    VariantsType,
    VariantZygosityCountCollection,
)
from snpdb.sample_filters import get_sample_ontology_q, get_sample_qc_gene_list_gene_symbol_q
from snpdb.views.datatable_view import DC, CellData, DatatableConfig, RichColumn, SortOrder
from uicore.templatetags.js_tags import jsonify_for_js
from variantgrid.perm_path import get_visible_url_names


def url_if_visible(url_name: str, **kwargs) -> Optional[str]:
    """ Deployments without patients (eg Shariant) unregister these urls entirely """
    if get_visible_url_names().get(url_name):
        return reverse(url_name, kwargs=kwargs)
    return None


class VCFListColumns(DatatableConfig[VCF]):
    server_csv_download = True

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.scroll_x = True

        self.rich_columns = [
            RichColumn(key="id", visible=False),
            RichColumn(key="name", label="Name", orderable=True,
                       renderer=self.view_primary_key, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="date", label="Date", orderable=True, default_sort=SortOrder.DESC,
                       css_class="text-nowrap", client_renderer='TableFormat.timestamp'),
            RichColumn(key="import_status", label="Import Status", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(ImportStatus.choices)),
            RichColumn(key="data_archived_date", label="Archived", orderable=True,
                       css_class="text-nowrap", client_renderer='TableFormat.timestamp'),
            RichColumn(key="genome_build__name", label="Genome Build", orderable=True),
            RichColumn(key="user__username", label="Uploaded by", orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn(key="source", label="VCF source", orderable=True),
            RichColumn(key="uploadedvcf__file_upload__import_source", label="Import Source", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(ImportSource.choices)),
            RichColumn(key="genotype_samples", label="Genotype Samples", orderable=True),
            RichColumn(key="project__name", label="Project", orderable=True),
            RichColumn(key="uploadedvcf__vcf_importer__name", label="VCF Importer", orderable=True,
                       enabled=self.user.is_superuser),
            RichColumn(key="uploadedvcf__vcf_importer__version", label="VCF Importer Version", orderable=True,
                       enabled=self.user.is_superuser),
            RichColumn(key="id", name="delete", label="", orderable=False,
                       renderer=self.render_delete, client_renderer='TableFormat.deleteRow'),
        ]

    def get_initial_queryset(self) -> QuerySet[VCF]:
        user_grid_config = UserGridConfig.get(self.user, 'VCFs')
        return VCF.filter_for_user(self.user, group_data=user_grid_config.show_group_data)

    def filter_queryset(self, qs: QuerySet[VCF]) -> QuerySet[VCF]:
        if project := self.get_query_param("project"):
            qs = qs.filter(project=project)
        if genome_build_name := self.get_query_param("genome_build_name"):
            genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
            qs = qs.filter(genome_build=genome_build)
        return qs


class SamplesListColumns(DatatableConfig[Sample]):
    server_csv_download = True
    # The unfiltered count is over a correlated subquery and a group by, and only feeds the
    # "(filtered from N total)" text
    count_unfiltered = False

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.scroll_x = True

        # Only show columns that have data somewhere in what this user can see
        qs = self._sample_queryset
        has_mutational_signature = qs.filter(mutationalsignature__isnull=False).exists()
        has_somalier_ancestry = qs.filter(somaliersampleextract__somalierancestry__isnull=False).exists()
        has_sample_gene_lists = qs.filter(samplegenelist__isnull=False).exists()

        if settings.PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB:
            dob_label = "Age"
            dob_renderer = self._render_age
            dob_client_renderer = None
        else:
            dob_label = "D.O.B."
            dob_renderer = None
            dob_client_renderer = 'TableFormat.timestamp'

        self.rich_columns = [
            RichColumn(key="id", visible=False),
            RichColumn(key="name", label="Name", orderable=True,
                       renderer=self.view_primary_key, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="het_hom_count", label="Het/Hom Count", orderable=True),
            RichColumn(key="vcf__date", label="Date", orderable=True, default_sort=SortOrder.DESC,
                       css_class="text-nowrap", client_renderer='TableFormat.timestamp'),
            RichColumn(key="import_status", label="Import Status", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(ImportStatus.choices)),
            RichColumn(key="vcf__genome_build__name", label="Genome Build", orderable=True),
            RichColumn(key="variants_type", label="Variants Type", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(VariantsType.choices)),
            RichColumn(key="vcf__user__username", label="Uploaded by", orderable=True,
                       extra_columns=["vcf__user__id"], renderer=self.render_user),
            RichColumn(key="vcf__source", label="VCF source", orderable=True),
            RichColumn(key="vcf__name", label="VCF Name", orderable=True, extra_columns=["vcf__id"],
                       renderer=self._render_vcf, client_renderer='renderOptionalLink'),
            RichColumn(key="vcf__project__name", label="Project", orderable=True),
            RichColumn(key="vcf__uploadedvcf__file_upload__import_source", label="Import Source", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(ImportSource.choices)),
            RichColumn(key="sample_gene_list_count", label="# Sample GeneLists", orderable=True,
                       extra_columns=["activesamplegenelist__id"], enabled=has_sample_gene_lists,
                       renderer=self._render_sample_gene_list_count,
                       client_renderer='renderSampleGeneListCount'),
            RichColumn(key="mutationalsignature__summary", label="Mutational Signature", orderable=True,
                       extra_columns=["mutationalsignature__id"], enabled=has_mutational_signature,
                       renderer=self._render_mutational_signature, client_renderer='renderOptionalLink'),
            RichColumn(key="somaliersampleextract__somalierancestry__predicted_ancestry",
                       label="Predicted Ancestry", orderable=True, enabled=has_somalier_ancestry,
                       client_renderer=RichColumn.choices_client_renderer(SuperPopulationCode.choices)),
            RichColumn(key="patient__patient_code", label="Patient Code", orderable=True),
            RichColumn(key="patient__first_name", label="First Name", orderable=True),
            RichColumn(key="patient__last_name", label="Last Name", orderable=True),
            RichColumn(key="patient__sex", label="Sex", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(Sex.choices)),
            RichColumn(key="patient__date_of_birth", label=dob_label, orderable=True,
                       renderer=dob_renderer, client_renderer=dob_client_renderer),
            RichColumn(key="extraction__specimen__reference_id", label="Specimen", orderable=True,
                       extra_columns=["extraction__specimen__id"],
                       renderer=self._render_specimen, client_renderer='renderOptionalLink'),
            # The DNA and RNA arms of one block share a specimen, so the specimen reference alone can't
            # tell those rows apart
            RichColumn(key="extraction__reference_id", label="Extraction", orderable=True,
                       extra_columns=["extraction__id"],
                       renderer=self._render_extraction, client_renderer='renderOptionalLink'),
            RichColumn(key="extraction__specimen__tissue__name", label="Tissue", orderable=True),
            RichColumn(key="extraction__specimen__collection_date", label="Collected", orderable=True,
                       css_class="text-nowrap", client_renderer='TableFormat.timestamp'),
            RichColumn(key="id", name="delete", label="", orderable=False,
                       renderer=self.render_delete, client_renderer='TableFormat.deleteRow'),
        ]

    @staticmethod
    def _render_age(cell: CellData) -> JsonDataType:
        return calculate_age(cell.value)

    @staticmethod
    def _render_vcf(cell: CellData) -> JsonDataType:
        return {"text": cell.value, "url": url_if_visible("view_vcf", vcf_id=cell["vcf__id"])}

    @staticmethod
    def _render_sample_gene_list_count(cell: CellData) -> JsonDataType:
        return {"count": cell.value, "active": bool(cell["activesamplegenelist__id"])}

    @staticmethod
    def _render_mutational_signature(cell: CellData) -> JsonDataType:
        url = None
        if mutational_signature_id := cell["mutationalsignature__id"]:
            url = url_if_visible("view_mutational_signature",
                                  mutational_signature_id=mutational_signature_id)
        return {"text": cell.value, "url": url}

    @staticmethod
    def _render_specimen(cell: CellData) -> JsonDataType:
        return SamplesListColumns._optional_link(cell, "view_specimen", "extraction__specimen__id",
                                                 "specimen_id")

    @staticmethod
    def _render_extraction(cell: CellData) -> JsonDataType:
        return SamplesListColumns._optional_link(cell, "view_extraction", "extraction__id", "extraction_id")

    @staticmethod
    def _optional_link(cell: CellData, url_name: str, pk_column: str, url_kwarg: str) -> JsonDataType:
        """ A sample doesn't have to have come through an extraction, and an extraction's own
            reference is optional - so only draw a link when there's something to link to """
        pk = cell[pk_column]
        if not pk:
            return {"text": cell.value}
        return {"text": cell.value or f"({pk})", "url": url_if_visible(url_name, **{url_kwarg: pk})}

    @cached_property
    def _sample_queryset(self) -> QuerySet[Sample]:
        user_grid_config = UserGridConfig.get(self.user, 'Samples')
        return Sample.filter_for_user(self.user, group_data=user_grid_config.show_group_data)

    def get_initial_queryset(self) -> QuerySet[Sample]:
        # het_hom_count comes from the per-sample CohortGenotypeStats row
        # (sample IS NOT NULL, filter_key NULL, passing_filter=False).
        cgs_subquery = (CohortGenotypeStats.objects
                        .filter(sample=OuterRef("pk"),
                                filter_key__isnull=True, passing_filter=False)
                        .annotate(het_plus_hom=F("het_count") + F("hom_count"))
                        .values("het_plus_hom")[:1])
        return self._sample_queryset.annotate(
            sample_gene_list_count=Count("samplegenelist", distinct=True),
            het_hom_count=Subquery(cgs_subquery, output_field=IntegerField()))

    def filter_queryset(self, qs: QuerySet[Sample]) -> QuerySet[Sample]:
        if project := self.get_query_param("project"):
            qs = qs.filter(vcf__project=project)
        if genome_build_name := self.get_query_param("genome_build_name"):
            genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
            qs = qs.filter(vcf__genome_build=genome_build)
        if (variants_type := self.get_query_json("variants_type")) is not None:
            qs = qs.filter(variants_type__in=variants_type)
        return qs


class AbstractSkippedAnnotationColumns(DatatableConfig[Variant]):
    """ Shows Variants that VEP was unable to annotate (variantannotation__vep_skipped_reason set).
        Subclasses provide a variant source (VCF/Sample - anything with get_variant_qs) via
        _get_variant_source(). """
    grid_name = "Skipped Annotation"

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="id", visible=False),
            RichColumn(key="variant_string", label="Variant", orderable=True, default_sort=SortOrder.ASC,
                       renderer=self._render_variant, extra_columns=["id"],
                       client_renderer='renderOptionalLink'),
            RichColumn(key="variantannotation__vep_skipped_reason", label="Skipped Reason", orderable=True),
            RichColumn(key="variantannotation__annotation_run_id", label="Annotation Run", orderable=True,
                       renderer=self._render_annotation_run, client_renderer='renderOptionalLink'),
        ]

    def _get_variant_source(self) -> tuple[Any, GenomeBuild]:
        """ (object with get_variant_qs, the build its variants are in) """
        raise NotImplementedError()

    @staticmethod
    def _render_variant(cell: CellData) -> JsonDataType:
        return {"text": cell.value, "url": reverse("view_variant", kwargs={"variant_id": cell["id"]})}

    @staticmethod
    def _render_annotation_run(cell: CellData) -> JsonDataType:
        if annotation_run_id := cell.value:
            return {"text": f"AnnotationRun {annotation_run_id}",
                    "url": reverse("view_annotation_run", kwargs={"annotation_run_id": annotation_run_id})}
        return None

    def get_initial_queryset(self) -> QuerySet[Variant]:
        variant_source, genome_build = self._get_variant_source()
        qs = get_queryset_for_latest_annotation_version(Variant, genome_build)
        qs = variant_source.get_variant_qs(qs).filter(variantannotation__vep_skipped_reason__isnull=False)
        return Variant.annotate_variant_string(qs)


class SampleSkippedAnnotationColumns(AbstractSkippedAnnotationColumns):
    def _get_variant_source(self) -> tuple[Any, GenomeBuild]:
        sample = Sample.get_for_user(self.user, self.get_query_param("sample_id"))
        return sample, sample.genome_build


class VCFSkippedAnnotationColumns(AbstractSkippedAnnotationColumns):
    def _get_variant_source(self) -> tuple[Any, GenomeBuild]:
        vcf = VCF.get_for_user(self.user, self.get_query_param("vcf_id"))
        return vcf, vcf.genome_build


class CohortSampleListColumns(DatatableConfig[Sample]):
    """ Sample picker on the cohort page - either the cohort's samples, or the ones it could add """
    grid_name = "Cohort Samples"
    SHOW_COHORT = "show_cohort"
    EXCLUDE_COHORT = "exclude_cohort"

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key="id", label="ID", orderable=True, default_sort=SortOrder.DESC,
                       renderer=self.view_primary_key, client_renderer='renderCohortSampleCheckbox'),
            RichColumn(key="name", label="Name", orderable=True),
            RichColumn(key="vcf__name", label="VCF", orderable=True),
            RichColumn(key="patient__family_code", label="Family Code", orderable=True),
            RichColumn(key="patient__patient_code", label="Patient Code", orderable=True),
            RichColumn(key="patient__first_name", label="First Name", orderable=True),
            RichColumn(key="patient__last_name", label="Last Name", orderable=True),
            RichColumn(key="patient__sex", label="Sex", orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(Sex.choices)),
            RichColumn(key="patient__date_of_birth", label="D.O.B.", orderable=True,
                       client_renderer='TableFormat.timestamp'),
        ]

    def get_initial_queryset(self) -> QuerySet[Sample]:
        cohort = Cohort.get_for_user(self.user, self.get_query_param("cohort_id"))
        qs = Sample.filter_for_user(self.user)
        qs = qs.filter(vcf__genome_build=cohort.genome_build, import_status=ImportStatus.SUCCESS)

        cohort_op = self.get_query_param("cohort_op") or self.EXCLUDE_COHORT
        cohort_q = Q(cohortsample__cohort=cohort)
        if cohort_op == self.SHOW_COHORT:
            pass
        elif cohort_op == self.EXCLUDE_COHORT:
            cohort_q = ~cohort_q
        else:
            raise ValueError(f"Unknown cohort_op: '{cohort_op}'")
        return qs.filter(cohort_q)


class CohortListColumns(DatatableConfig[Cohort]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='import_status', label='Status', orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(ImportStatus.choices)),
            RichColumn(key="user__username", label='User', orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn(key='modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
            RichColumn(key='sample_count', label='Sample Count', orderable=True),
            RichColumn(key='id', name='delete', label='', orderable=False,
                       renderer=self.render_delete,
                       client_renderer='TableFormat.deleteRow'),
        ]

    def get_initial_queryset(self) -> QuerySet[Cohort]:
        return Cohort.filter_for_user(self.user, success_status_only=False).filter(vcf__isnull=True)

    def filter_queryset(self, qs: QuerySet[Cohort]) -> QuerySet[Cohort]:
        user_grid_config = UserGridConfig.get(self.user, 'Cohorts')
        if not user_grid_config.show_group_data:
            qs = qs.filter(user=self.user)
        return qs


class FamilyGroupListColumns(DatatableConfig[DC]):
    """ Trios/Quads listing - same grid bar the extra family members """
    MODEL: type[DC]
    GRID_NAME: str
    # (field prefix, label, has an affected column)
    FAMILY_MEMBERS = [("mother", "Mother", True), ("father", "Father", True), ("proband", "Proband", False)]

    def __init__(self, request: HttpRequest):
        super().__init__(request)
        member_columns = []
        for member, label, has_affected in self.FAMILY_MEMBERS:
            member_columns.append(RichColumn(key=f'{member}__sample__name', label=label, orderable=True))
            if has_affected:
                member_columns.append(RichColumn(key=f'{member}_affected', label=f'{label} Affected',
                                                 orderable=True))
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key="user__username", label='User', orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn(key='modified', client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
            *member_columns,
            RichColumn(key='id', name='delete', label='', orderable=False,
                       renderer=self.render_delete,
                       client_renderer='TableFormat.deleteRow'),
        ]

    def get_initial_queryset(self) -> QuerySet[DC]:
        return self.MODEL.filter_for_user(self.user)

    def filter_queryset(self, qs: QuerySet[DC]) -> QuerySet[DC]:
        user_grid_config = UserGridConfig.get(self.user, self.GRID_NAME)
        if not user_grid_config.show_group_data:
            qs = qs.filter(user=self.user)
        return qs


class TriosListColumns(FamilyGroupListColumns[Trio]):
    MODEL = Trio
    GRID_NAME = 'Trios'


class QuadsListColumns(FamilyGroupListColumns[Quad]):
    MODEL = Quad
    GRID_NAME = 'Quads'
    FAMILY_MEMBERS = [*FamilyGroupListColumns.FAMILY_MEMBERS, ("sibling", "Sibling", True)]


class GenomicIntervalsListColumns(DatatableConfig[GenomicIntervalsCollection]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)
        self.rich_columns = [
            RichColumn(key='id', visible=False),
            RichColumn(key='name', label='Name', orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key='import_status', label='Status', orderable=True,
                       client_renderer=RichColumn.choices_client_renderer(ImportStatus.choices)),
            RichColumn(key='genome_build__name', label='Genome Build', orderable=True),
            RichColumn(key="user__username", label='Uploaded by', orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn(key='id', name='delete', label='', orderable=False,
                       renderer=self.render_delete,
                       client_renderer='TableFormat.deleteRow'),
        ]

    def get_initial_queryset(self) -> QuerySet[GenomicIntervalsCollection]:
        return get_objects_for_user(self.user, 'snpdb.view_genomicintervalscollection', accept_global_perms=False)


class NamedCollectionColumns(DatatableConfig[DC]):
    """ A user's named collections - nothing to show but who made it and when """
    MODEL: type[DC]

    def __init__(self, request):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn(key="id", visible=False),
            RichColumn(key="name", label="Name", orderable=True,
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key="user__username", label="User", orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="modified", client_renderer='TableFormat.timestamp', orderable=True,
                       default_sort=SortOrder.DESC),
        ]

    def get_initial_queryset(self) -> QuerySet[DC]:
        return self.MODEL.filter_for_user(self.user)


class CustomColumnsCollectionColumns(NamedCollectionColumns[CustomColumnsCollection]):
    MODEL = CustomColumnsCollection


def render_clingen_allele(cell: CellData) -> JsonDataType:
    if ca_id := cell.value:
        ca_id = ClinGenAllele.format_clingen_allele(ca_id)
    return ca_id


def render_exon_and_intron(cell: CellData) -> JsonDataType:
    """ MS Excel will turn '8/11' into a date :( """
    if val := cell.value:
        val = val.replace("/", " of ")
    return val


def render_annotsv_pathogenic_overlaps(cell: CellData) -> JsonDataType:
    """ AnnotSV's nested {event type: {source/phen/hpo/coord}} JSON -> one entry per event type """
    text = None
    if overlaps := cell.value:
        events = []
        for event, data in overlaps.items():
            details = " / ".join(str(v) for v in data.values())
            events.append(f"{event}: {details}")
        text = ", ".join(events)
    return text


class AbstractVariantGrid(DatatableConfig[Variant]):
    """ The variant grids - the analysis node grid and the standalone Variant tables. Their columns are
        built per user from a CustomColumnsCollection (@see snpdb.grid_columns.custom_columns) rather
        than declared, and the same config serves the grid page and the CSV/VCF exports. """
    VARIANT_COLUMN_NAME = "id"  # the mandatory Variant column
    # These grids are wide (dozens of columns) and hold cells with hundreds of links, so the client lays
    # them out table-layout: fixed and clips each cell to one line. That needs every column to carry a
    # width. @see .variantgrid-datatable in global.scss
    table_class = "variantgrid-datatable"
    default_column_width = 150
    # Counting the unfiltered queryset is as expensive as the page itself
    count_unfiltered = False
    compact_controls = True  # every pixel of chrome is a row the user can't see
    ajax_type = 'GET'  # keeps @cache_page on a data endpoint working
    cache_stable_params = True
    scroll_x = True
    server_csv_download = True
    filter_builder = True
    approximate_count_enabled = True
    expand_prefetch = False  # rows are full of links; the arrow is affordance enough
    # The variant grids call the column 'id' - say what it's the id of
    csv_label_overrides = {VARIANT_COLUMN_NAME: "variant_id"}
    csv_name = "Variant"
    analysis_tags = False

    def __init__(self, request: HttpRequest, af_show_in_percent: Optional[bool] = None):
        super().__init__(request)
        if af_show_in_percent is None:
            af_show_in_percent = settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT
        self.af_show_in_percent = af_show_in_percent
        self.rich_columns = self._get_rich_columns()
        # Row expansion (@see DataTableDefinition.setupClientExpend) - variantGridRowDetail in grid.js
        # fetches variant_grid_row_detail for the annotation version this grid is showing
        self.expand_client_renderer = \
            f"variantGridRowDetail.bind(null, {self._get_annotation_version().pk})"

    def _get_rich_columns(self) -> list[RichColumn]:
        rich_columns, _sample_columns_position = self._build_variant_grid_columns()
        return rich_columns

    def _build_variant_grid_columns(self) -> tuple[list[RichColumn], Optional[int]]:
        """ (this user's columns, where sample columns go) - @see snpdb.grid_columns.custom_columns """
        return get_variant_grid_columns(self._get_custom_columns_collection(),
                                        self._get_annotation_version(),
                                        self._get_standard_overrides(self.af_show_in_percent),
                                        analysis_tags=self.analysis_tags)

    def _get_custom_columns_collection(self) -> CustomColumnsCollection:
        return UserSettings.get_for_user(self.user).columns

    def _get_standard_overrides(self, af_show_in_percent: bool) -> dict[str, dict]:
        overrides = {
            # Note:     client side renderers should only be used for adding links etc, never conversion of data, such as
            #           unit to percent, as the CSV downloads (w/o JS renderers) won't match the grid.
            # The representative variant cell: expand arrow, select checkbox, cascade label (details link)
            'id': {'width': 280, 'client_renderer': 'VariantGridFormat.representativeVariant'},
            'classifications': {
                'css_class': 'no-word-wrap',
                'client_renderer': 'VariantGridFormat.classifications',
            },
            'tags_global': {
                'model_field': False, 'css_class': 'no-word-wrap', 'orderable': False,
                'client_renderer': 'VariantGridFormat.tagsGlobal',
            },
            'clinvar__clinvar_variation_id': {'width': 60, 'client_renderer': 'VariantGridFormat.clinvarLink'},
            'variantallele__allele__clingen_allele__id': {
                'width': 90,
                'renderer': render_clingen_allele, 'csv_rendered': True,
                'client_renderer': 'VariantGridFormat.clinGenAlleleId',
            },
            'variantannotation__cosmic_id': {'width': 130, 'client_renderer': 'VariantGridFormat.cosmicLink'},
            'variantannotation__cosmic_legacy_id': {'width': 130, 'client_renderer': 'VariantGridFormat.cosmicLink'},
            'variantannotation__dbsnp_rs_id': {'width': 130, 'client_renderer': 'VariantGridFormat.dbsnp'},
            'variantannotation__pubmed': {'client_renderer': 'VariantGridFormat.pubMed'},
            'variantannotation__gene__geneannotation__hpo_terms': {'client_renderer': 'VariantGridFormat.ontologyTerms'},
            'variantannotation__gene__geneannotation__mondo_terms': {'client_renderer': 'VariantGridFormat.ontologyTerms'},
            'variantannotation__gene__geneannotation__omim_terms': {'client_renderer': 'VariantGridFormat.ontologyTerms'},
            'variantannotation__transcript_version__gene_version__gene_symbol__symbol': {
                'client_renderer': 'VariantGridFormat.geneSymbolLink'},
            'variantannotation__overlapping_symbols': {'client_renderer': 'VariantGridFormat.geneSymbolNewWindowLink'},
            'variantannotation__transcript_version__gene_version__hgnc__omim_ids': {
                'width': 60, 'client_renderer': 'VariantGridFormat.omimLink'},
            # A member shown standalone keeps its own link - and the composite cells reuse these
            'variantannotation__gnomad_filtered': {'client_renderer': 'VariantGridFormat.gnomadFiltered'},
            # Composite cells, keyed by the composite column's own name - the members they draw
            # ride along hidden and their sort menus come from CompositeColumnMember
            'consequence_impact': {'client_renderer': 'VariantGridFormat.impactConsequence'},
            'gnomad': {'client_renderer': 'VariantGridFormat.gnomad'},
            'spliceai': {'client_renderer': 'VariantGridFormat.spliceai'},
            'maxentscan': {'client_renderer': 'VariantGridFormat.maxentscan'},
            'mastermind': {'client_renderer': 'VariantGridFormat.mastermind'},
            'aloft': {'client_renderer': 'VariantGridFormat.aloft'},
            'predictions': {'client_renderer': 'VariantGridFormat.predictions'},
            # The same cell (and the same menu) the cohort node draws with its own counts
            # @see CohortNode._get_node_extra_columns
            'db_zygosity': {'client_renderer': 'VariantGridFormat.dbZygosityCounts'},
            'variantannotation__exon': {'renderer': render_exon_and_intron, 'csv_rendered': True},
            'variantannotation__intron': {'renderer': render_exon_and_intron, 'csv_rendered': True},
            'variantannotation__mastermind_mmid3': {'client_renderer': 'VariantGridFormat.masterMind'},
            'variantannotation__mavedb_urn': {'client_renderer': 'VariantGridFormat.mavedbUrn'},
            'variantannotation__annotsv_pathogenic_overlaps': {
                'renderer': render_annotsv_pathogenic_overlaps, 'csv_rendered': True,
            },
        }

        if af_show_in_percent:
            # gnomAD etc are all stored as AF in DB - want to show as percentage on grid
            # But need to be able to turn it off to export VCF as AF
            render_unit_af = get_allele_frequency_formatter(source_in_percent=False,
                                                            dest_in_percent=af_show_in_percent)
            af_override = {'renderer': render_unit_af, 'csv_rendered': True}
            af_columns = [
                'variantannotation__af_1kg',
                'variantannotation__af_uk10k',
                'variantannotation__gnomad2_liftover_af',
                'variantannotation__gnomad_af',
                'variantannotation__gnomad_afr_af',
                'variantannotation__gnomad_amr_af',
                'variantannotation__gnomad_asj_af',
                'variantannotation__gnomad_eas_af',
                'variantannotation__gnomad_fin_af',
                'variantannotation__gnomad_nfe_af',
                'variantannotation__gnomad_oth_af',
                'variantannotation__gnomad_popmax_af',
                'variantannotation__gnomad_sas_af',
                'variantannotation__topmed_af',
            ]
            for column in af_columns:
                overrides.setdefault(column, {}).update(af_override)
        return overrides

    def _get_base_queryset(self) -> QuerySet:
        raise NotImplementedError()

    def get_extra(self) -> JsonObjType:
        # gnomAD links are per genome build, and the client renderers have no other way to know it
        extra = {"genomeBuild": self.genome_build.name,
                 "clinvarStars": dict(ClinVarReviewStatus.STARS),
                 # What counts as a bad GQ/PL in the sample genotype cell
                 "genotypeQuality": settings.VARIANT_GRID_GENOTYPE_QUALITY_THRESHOLDS}
        # The AF the import 'common' filter uses, in the units the grid shows AFs in - the gnomAD
        # cell mutes at or above it so rare variants keep the reader's full attention
        if cf_data := settings.VCF_IMPORT_COMMON_FILTERS.get(self.genome_build.name):
            common_af = cf_data["gnomad_af_min"]
            if settings.VARIANT_ALLELE_FREQUENCY_CLIENT_SIDE_PERCENT:
                common_af *= 100
            extra["commonGnomadAf"] = common_af
        return extra

    def get_table_classes(self) -> list[str]:
        """ Two line rows are a per-user setting. The second line is in the markup either way -
            CSS decides whether it shows, so switching it doesn't need the rows re-rendering """
        classes = super().get_table_classes()
        if UserSettings.get_for_user(self.user).variant_grid_two_line_rows:
            classes.append("two-line-rows")
        return classes

    def _get_annotation_version(self) -> AnnotationVersion:
        return self.annotation_version

    def initial_order(self) -> Optional[list]:
        """ Only a column that asked to be the initial sort - unlike other tables, falling back to
            the first column would sort the whole variant table on every first load """
        for rc in self.enabled_columns:
            if rc.default_sort and rc.visible and rc.orderable:
                return [[self.column_index(rc), "desc" if rc.default_sort == SortOrder.DESC else "asc"]]
        return None

    def _genomic_order_by(self, qs: QuerySet, descending: bool) -> QuerySet:
        """ The Variant column sorts in genome build order whatever it displays. Contig order is a CASE over
            the build's standard contigs rather than a join through GenomeBuildContig - MT is shared between
            builds and the join would return the row once per build. Non-standard contigs sort after, by id """
        whens = [When(locus__contig_id=contig_id, then=Value(i))
                 for i, contig_id in enumerate(self.genome_build.standard_contigs.values_list("pk", flat=True))]
        qs = qs.annotate(_contig_order=Case(*whens, default=Value(len(whens)), output_field=IntegerField()))
        fields = ["_contig_order", "locus__contig_id", "locus__position", "locus__ref__seq", "alt__seq", "pk"]
        prefix = "-" if descending else ""
        return qs.order_by(*[f"{prefix}{f}" for f in fields])

    def ordering(self, qs: QuerySet) -> QuerySet:
        for rich_column, desc in self.requested_ordering():
            if rich_column.name == self.VARIANT_COLUMN_NAME:
                return self._genomic_order_by(qs, descending=desc)
        return super().ordering(qs)

    def _get_permission_user(self) -> User:
        return self.user

    def get_initial_queryset(self) -> QuerySet[Variant]:
        qs = self._get_base_queryset()
        # Restrict the variantallele join to this grid's genome build. Some contigs are shared between builds
        # (e.g. MT / NC_012920 is shared by GRCh37 & GRCh38) so the same Variant has a VariantAllele per build -
        # without this, grid columns that join through 'variantallele' (e.g. the ClinGen Allele ID) return that
        # variant once per build (duplicate rows / inflated counts vs the node count, which doesn't make that join).
        # @see https://github.com/SACGF/variantgrid/issues/1626
        qs = qs.filter(Q(variantallele__isnull=True) | Q(variantallele__genome_build=self.genome_build))
        # Annotate so we can use global_variant_zygosity in grid columns
        qs, _ = VariantZygosityCountCollection.annotate_global_germline_counts(qs)
        # Column filter rules are applied by apply_filters on our result - doing it here as well
        # adds a second JOIN per filtered relation, which multiplies rows
        if q := self._get_q():
            qs = qs.filter(q)
        return qs.annotate(**self._get_grid_only_annotation_kwargs())

    def _get_grid_only_annotation_kwargs(self) -> dict:
        """ Things not used in counts etc - only to display grid """
        return get_variantgrid_extra_annotate(self._get_permission_user())

    def _get_q(self) -> Optional[Q]:
        return None

class TagColorsCollectionColumns(NamedCollectionColumns[TagColorsCollection]):
    MODEL = TagColorsCollection


class LiftoverRunColumns(DatatableConfig[LiftoverRun]):

    def __init__(self, request):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn(key="id", renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            # Don't bother with user as always admin/admin bot
            # RichColumn(key="user__username", label="User", orderable=True),
            RichColumn(key="created", client_renderer='TableFormat.timestamp',
                       orderable=True, default_sort=SortOrder.DESC),
            # Don't bother with modified as always same time as created
            # RichColumn(key='allele_source', orderable=True),
            RichColumn(key='conversion_tool', renderer=self.render_conversion_tool, orderable=True),
            RichColumn(key='num_alleles', label='Num Alleles', orderable=True, css_class="num"),
            RichColumn(key='source_vcf', orderable=True, css_class="formatted-text"),
            RichColumn(key='source_genome_build', label='Source Build', orderable=True),
            RichColumn(key='genome_build', label='Dest Build', orderable=True),
            RichColumn(key="uploadedliftover__file_upload__uploadpipeline__status",
                       label='Status', renderer=self.render_import_status, orderable=True),
            RichColumn(key="uploadedliftover__file_upload__uploadpipeline__items_processed",
                       label='Processed', orderable=True, css_class="num"),
        ]

    def render_conversion_tool(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        if conversion_tool := row['conversion_tool']:
            act = AlleleConversionTool(conversion_tool)
            label = act.label
        return label

    def render_import_status(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        if status := row['uploadedliftover__file_upload__uploadpipeline__status']:
            processing_status = ProcessingStatus(status)
            label = processing_status.label
        return label

    def get_initial_queryset(self) -> QuerySet[LiftoverRun]:
        qs = LiftoverRun.objects.all()
        return qs.annotate(num_alleles=Count("alleleliftover"))


class AbstractAlleleLiftoverColumns(DatatableConfig[AlleleLiftover]):

    def __init__(self, request):
        super().__init__(request)
        self.user = request
        self.liftover_run = None

        self.rich_columns = [
            RichColumn(key="allele", orderable=True,
                       renderer=self.render_allele, client_renderer='TableFormat.linkUrl'),
            RichColumn(key="status", label="Status", renderer=self.render_status, orderable=True),
            RichColumn(key="data", label="Data", renderer=self.render_data_json, orderable=True),
            RichColumn(key="error", label="Error", renderer=self.render_error_json, orderable=True),
        ]

    def render_allele(self, row: dict[str, Any]) -> JsonDataType:
        data = {}
        if allele_id := row['allele']:
            allele = get_object_or_404(Allele, id=allele_id)
            data = {
                "text": str(allele),
                "url": allele.get_absolute_url(),
            }
        return data

    def render_status(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        current = ""
        if status := row['status']:
            processing_status = ProcessingStatus(status)
            label = processing_status.label

        has_build_to_icon = {
            True: "✅",
            False: "❌"
        }

        if allele_id := row['allele']:
            allele: Allele
            if allele := Allele.objects.filter(id=allele_id).first():
                has_37 = bool(allele.variant_for_build_optional(GenomeBuild.grch37()))
                has_38 = bool(allele.variant_for_build_optional(GenomeBuild.grch38()))
                label += f" (Current: {has_build_to_icon[has_37]} GRCh37, {has_build_to_icon[has_38]} GRCh38)"

        return label

    def render_data_json(self, row: dict[str, Any]) -> JsonDataType:
        if js := row["data"]:
            if INFO_LIFTOVER_SWAPPED_REF_ALT in str(js):
                return "Swapped Ref/Alt due to SWAP=1"
        if js is None:
            return "-"
        return jsonify_for_js(js, pretty=True)

    def render_error_json(self, row: dict[str, Any]) -> JsonDataType:
        if js := row["error"]:
            if "message" in js and len(js.keys()) == 1:
                return js.get("message")
        if js is None:
            return "-"
        return jsonify_for_js(js, pretty=True)


class LiftoverRunAlleleLiftoverColumns(AbstractAlleleLiftoverColumns):
    def get_initial_queryset(self) -> QuerySet[AlleleLiftover]:
        liftover_run_id = self.get_query_param("liftover_run_id")
        if liftover_run_id is None:
            raise ValueError("liftover_run_id not provided")
        liftover_run = get_object_or_404(LiftoverRun, pk=liftover_run_id)
        return AlleleLiftover.objects.filter(liftover=liftover_run)


class AlleleLiftoverFailureColumns(AbstractAlleleLiftoverColumns):
    def __init__(self, request):
        super().__init__(request)
        ct_column = RichColumn(key="liftover__conversion_tool",
                               label="Conversion Tool", renderer=self.render_conversion_tool,
                               orderable=True)
        self.rich_columns.insert(1, ct_column)

    def render_conversion_tool(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        if status := row['liftover__conversion_tool']:
            conversion_tool = AlleleConversionTool(status)
            label = conversion_tool.label
        return label

    def get_initial_queryset(self) -> QuerySet[AlleleLiftover]:
        genome_build_name = self.get_query_param("genome_build_name")
        if genome_build_name is None:
            raise ValueError("genome_build_name not provided")
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

        qs = AlleleLiftover.objects.filter(liftover__genome_build=genome_build,
                                           status=ProcessingStatus.ERROR)
        qs = qs.annotate(
            max_id=Max('allele__alleleliftover__id')
        ).filter(
            id=F('max_id')
        )
        return qs


class ManualVariantEntryCollectionColumns(DatatableConfig[ManualVariantEntryCollection]):
    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.expand_client_renderer = DatatableConfig._row_expand_ajax('manual_variant_entry_collection_detail', expected_height=120)
        self.rich_columns = [
            RichColumn('id', orderable=True, default_sort=SortOrder.DESC),
            RichColumn('created', client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="user__username", orderable=True,
                       extra_columns=["user__id"], renderer=self.render_user),
            RichColumn('import_status', orderable=True, renderer=self._render_import_status),
            RichColumn('genome_build', orderable=True),
        ]

    def _render_import_status(self, row: dict[str, Any]) -> JsonDataType:
        label = ""
        if status := row['import_status']:
            import_status = ImportStatus(status)
            label = import_status.label
        return label

    def get_initial_queryset(self) -> QuerySet[ManualVariantEntryCollection]:
        qs = ManualVariantEntryCollection.objects.all()
        if not self.user.is_staff:
            qs = qs.filter(user=self.user)
        return qs


class SampleColumns(DatatableConfig[Sample]):
    """ This is currently only used on """
    def __init__(self, request):
        super().__init__(request)
        self.user = request.user

        self.rich_columns = [
            RichColumn(key="id",
                       renderer=self.view_primary_key,
                       client_renderer='TableFormat.linkUrl'),
            RichColumn(key="name", label="Name", orderable=True),
            RichColumn(key="vcf__name", label="VCF", orderable=True),
            RichColumn(key="vcf__date", client_renderer='TableFormat.timestamp', orderable=True),
            RichColumn(key="vcf__project__name", label="Project", orderable=True),
            RichColumn(key=OntologyService.OMIM, orderable=True),
            RichColumn(key=OntologyService.HPO, orderable=True),
            RichColumn(key=OntologyService.MONDO, orderable=True),
        ]

    def get_initial_queryset(self) -> QuerySet[Sample]:
        qs = Sample.filter_for_user(self.user)
        sample_patient_ontology_path = f"patient__{PATIENT_ONTOLOGY_TERM_PATH}"
        ontology_path = f"{sample_patient_ontology_path}__name"
        annotation_kwargs = {}
        for ot in [OntologyService.OMIM, OntologyService.HPO, OntologyService.MONDO]:
            q_ot = Q(**{f"{sample_patient_ontology_path}__ontology_service": ot})
            annotation_kwargs[ot.label] = StringAgg(ontology_path, Value('|'),
                                                    filter=q_ot, distinct=True,
                                                    output_field=TextField())
        return qs.annotate(**annotation_kwargs)

    def filter_queryset(self, qs: QuerySet[Sample]) -> QuerySet[Sample]:
        filters = []
        ontology_filters = []
        if ontology_terms := self.get_query_param('ontology_term_id'):
            if q:= get_sample_ontology_q(ontology_terms):
                ontology_filters.append(q)

        if ontology_filters:
            q = reduce(operator.or_, ontology_filters)
            filters.append(q)

        if gene_symbol_str := self.get_query_param("gene_symbol"):
            if q := get_sample_qc_gene_list_gene_symbol_q(gene_symbol_str):
                filters.append(q)

        direct_sample_select_filters = []
        if project := self.get_query_param("project"):
            direct_sample_select_filters.append(Q(vcf__project=project))

        if vcf := self.get_query_param("vcf"):
            direct_sample_select_filters.append(Q(vcf=vcf))

        if sample_str := self.get_query_param("sample"):
            samples = sample_str.split(",")
            direct_sample_select_filters.append(Q(pk__in=samples))

        if direct_sample_select_filters:
            q = reduce(operator.or_, direct_sample_select_filters)
            filters.append(q)

        if filters:
            qs = qs.filter(*filters)
        return qs
