from django.conf import settings
from django.db.models.aggregates import Count
from django.db.models.query_utils import Q
from functools import reduce
from guardian.shortcuts import get_objects_for_user
import operator

from library.database_utils import get_queryset_column_names, \
    get_queryset_select_from_where_parts
from library.jqgrid_sql import JqGridSQL
from library.jqgrid_user_row_config import JqGridUserRowConfig
from library.utils import calculate_age
from snpdb.models import VCF, Cohort, Sample, ImportStatus, \
    GenomicIntervalsCollection, CustomColumnsCollection, Variant, Trio, UserGridConfig
from snpdb.tasks.soft_delete_tasks import soft_delete_vcfs, remove_soft_deleted_vcfs_task


class VCFListGrid(JqGridUserRowConfig):
    model = VCF
    caption = 'VCFs'
    fields = ["id", "name", "date", "import_status", "user__username", "source",
              "uploadedvcf__uploaded_file__import_source", "genotype_samples", "project__name", "cohort__import_status",
              "uploadedvcf__vcf_importer__name", 'uploadedvcf__vcf_importer__version']
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'width': 550,
                 'formatter': 'viewVCFLink',
                 'formatter_kwargs': {"icon_css_class": "vcf-icon",
                                      "url_name": "view_vcf",
                                      "url_object_column": "id"}},
        'import_status': {'formatter': 'viewImportStatus'},
        'user__username': {'label': 'Uploaded by', 'width': 60},
        'source': {'label': 'VCF source'},
        "project__name": {'label': "Project"},
        'cohort__import_status': {'hidden': True},
        'uploadedvcf__vcf_importer__name': {"label": 'VCF Importer', "hide_non_admin": True},
        'uploadedvcf__vcf_importer__version': {"label": 'VCF Importer Version', "hide_non_admin": True},
    }

    def __init__(self, user):
        super().__init__(user)
        user_grid_config = UserGridConfig.get(user, self.caption)
        queryset = VCF.filter_for_user(user, group_data=user_grid_config.show_group_data)
        self.queryset = queryset.order_by("-pk").values(*self.get_field_names())
        self.extra_config.update({'shrinkToFit': False,
                                  'sortname': 'id',
                                  'sortorder': 'desc'})

    def delete_row(self, vcf_id):
        """ Do async as it may be slow """
        soft_delete_vcfs(self.user, vcf_id)


# TODO: Merge this an cohort grid below into 1
class SamplesListGrid(JqGridUserRowConfig):
    model = Sample
    caption = 'Samples'
    fields = ["id", "name", "import_status", "variants_type",
              "vcf__name", "vcf__user__username", "vcf__uploadedvcf__uploaded_file__import_source",
              "mutationalsignature__id", "mutationalsignature__summary", "samplestats__variant_count",
              "patient__first_name", "patient__last_name", "patient__sex", "patient__date_of_birth", "patient__date_of_death",
              "specimen__reference_id", "specimen__tissue", "specimen__collection_date"]
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {"width": 400,
                 'formatter': 'viewSampleLink',
                 'formatter_kwargs': {"icon_css_class": "sample-icon",
                                      "url_name": "view_sample",
                                      "url_object_column": "id"}},
        'vcf__name': {'label': 'VCF Name'},
        'import_status': {'formatter': 'viewImportStatus'},
        'mutationalsignature__id': {'hidden': True},
        'mutationalsignature__summary': {'label': 'Mutational Signature',
                                         'formatter': 'viewMutationalSignature'},
        'patient__last_name': {'label': 'Last Name'},
        'patient__sex': {'label': 'Sex'},
        'patient__date_of_birth': {'label': 'D.O.B.'},
        'patient__date_of_death': {'hidden': True},
        'vcf__user__username': {'label': 'Owner', 'width': 50},
        'vcf__uploadedvcf__uploaded_file__import_source': {'width': 55},
        'samplestats__variant_count': {'label': 'Variant Count', 'width': 50},
        "specimen__reference_id": {'label': 'Specimen'}
    }

    def __init__(self, user):
        super().__init__(user)

        user_grid_config = UserGridConfig.get(user, self.caption)
        queryset = Sample.filter_for_user(user, group_data=user_grid_config.show_group_data)

        # If you don't have permission to view a patient - blank it out
        # If you have read only and
        # TODO: We need to pass whole row in - as we need date of death to display age
        if settings.PATIENTS_READ_ONLY_SHOW_AGE_NOT_DOB:
            dob_colmodel = self._overrides.get('patient__date_of_birth', {})
            dob_colmodel['label'] = "Age"
            dob_colmodel['server_side_formatter'] = lambda row, field: calculate_age(row[field])
            self._overrides['patient__date_of_birth'] = dob_colmodel

        # Only show mut sig column if we have any
        if not queryset.filter(mutationalsignature__isnull=False).exists():
            mut_sig_colmodel = self._overrides.get('mutationalsignature__summary', {})
            mut_sig_colmodel['hidden'] = True
            self._overrides['mutationalsignature__summary'] = mut_sig_colmodel

        self.queryset = queryset.order_by("-pk").values(*self.get_field_names())
        self.extra_config.update({'shrinkToFit': False,
                                  'sortname': 'id',
                                  'sortorder': 'desc'})

    def delete_row(self, sample_id):
        """ Do async as it may take a few secs to delete """

        sample = Sample.get_for_user(self.user, sample_id)
        sample.check_can_write(self.user)
        Sample.objects.filter(pk=sample.pk).update(import_status=ImportStatus.MARKED_FOR_DELETION)
        task = remove_soft_deleted_vcfs_task.si()  # @UndefinedVariable
        task.apply_async()


class CohortSampleListGrid(JqGridUserRowConfig):
    model = Sample
    caption = 'Cohort Samples'
    fields = ["id", "name", "vcf__name", "patient__family_code",
              "patient__first_name", "patient__first_name",
              "patient__sex", "patient__date_of_birth"]
    colmodel_overrides = {'id': {'width': 20, 'formatter': 'viewSampleLink'},
                          'vcf__name': {'label': 'VCF'},
                          'patient__family_code': {'label': 'Family Code'},
                          'patient__first_name': {'label': 'First Name'},
                          'patient__last_name': {'label': 'Last Name'},
                          'patient__sex': {'label': 'Sex'},
                          'patient__date_of_birth': {'label': 'D.O.B.'}}

    def __init__(self, user, cohort_id, extra_filters=None):
        super().__init__(user)

        if extra_filters is None:
            extra_filters = {}

        cohort = Cohort.get_for_user(user, cohort_id)
        sample_filters = [Q(vcf__genome_build=cohort.genome_build),
                          Q(import_status=ImportStatus.SUCCESS)]
        SHOW_COHORT = "show_cohort"
        EXCLUDE_COHORT = "exclude_cohort"
        cohort_op = extra_filters.get("cohort_op", EXCLUDE_COHORT)
        cohort_q = Q(cohortsample__cohort=cohort)
        if cohort_op == SHOW_COHORT:
            pass
        elif cohort_op == EXCLUDE_COHORT:
            cohort_q = ~cohort_q
        else:
            raise ValueError(f"Unknown cohort_op: '{cohort_op}'")

        sample_filters.append(cohort_q)
        q = reduce(operator.and_, sample_filters)
        queryset = Sample.filter_for_user(user).filter(q).order_by("-pk")
        self.queryset = queryset.values(*self.get_field_names())


class CohortListGrid(JqGridUserRowConfig):
    model = Cohort
    caption = 'Cohorts'
    fields = ["id", "name", "import_status", "modified", "sample_count"]
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'formatter': 'linkFormatter',
                 'formatter_kwargs': {"icon_css_class": "cohort-icon",
                                      "url_name": "view_cohort",
                                      "url_object_column": "id"}},
    }

    def __init__(self, user):
        super().__init__(user)
        queryset = Cohort.filter_for_user(user, success_status_only=False).order_by("-pk")
        queryset = queryset.filter(vcf__isnull=True)  # Don't show auto-cohorts

        self.queryset = queryset.values(*self.get_field_names())
        self.extra_config.update({'sortname': "modified",
                                  'sortorder': "desc"})

    def get_colmodels(self, *args, **kwargs):
        colmodels = super().get_colmodels(*args, **kwargs)
        extra = {'index': 'sample_count', 'name': 'sample_count', 'label': 'Sample Count', 'sorttype': 'int'}
        colmodels.append(extra)
        return colmodels


class TriosListGrid(JqGridUserRowConfig):
    model = Trio
    caption = 'Trios'
    fields = ["id", "name", "mother__sample__name", "mother_affected",
              "father__sample__name", "father_affected", "proband__sample__name"]
    colmodel_overrides = {
        'id': {'formatter': 'linkFormatter',
               'formatter_kwargs': {"icon_css_class": "trio-icon",
                                    "display_column": "name",
                                    "url_name": "view_trio"}},
        "name": {"hidden": True},
        "mother__sample__name": {"label": "Mother"},
        "father__sample__name": {"label": "Father"},
        "proband__sample__name": {"label": "Proband"}
    }

    def __init__(self, user):
        super().__init__(user)
        queryset = Trio.filter_for_user(user).order_by("-pk")
        self.queryset = queryset.values(*self.get_field_names())


class GenomicIntervalsListGrid(JqGridUserRowConfig):
    model = GenomicIntervalsCollection
    caption = 'Genomic Intervals'
    fields = ["id", "name", "user__username", "import_status"]
    colmodel_overrides = {
        'id': {"hidden": True},
        "name": {'formatter': 'linkFormatter',
                 'formatter_kwargs': {"icon_css_class": "bed-icon",
                                      "url_name": "view_genomic_intervals",
                                      "url_object_column": "id"}},
        'user__username': {'label': 'Uploaded by'}
    }

    def __init__(self, user):
        super().__init__(user)
        queryset = get_objects_for_user(user, 'snpdb.view_genomicintervalscollection', accept_global_perms=False)
        self.queryset = queryset.order_by("-pk").values(*self.get_field_names())


class CustomColumnsCollectionListGrid(JqGridUserRowConfig):
    model = CustomColumnsCollection
    caption = 'Custom Columns'
    fields = ["id", "name", "user__username"]
    colmodel_overrides = {'id': {'width': 60, 'formatter': 'viewColumnLink'},
                          'user__username': {'label': 'Owner'}}

    def __init__(self, user):
        super().__init__(user)
        queryset = self.model.filter_for_user(user).order_by("-pk")
        queryset = queryset.annotate(num_columns=Count("customcolumn"))
        field_names = self.get_field_names() + ["num_columns"]
        self.queryset = queryset.values(*field_names)

    def get_colmodels(self, *args, **kwargs):
        colmodels = super().get_colmodels(*args, **kwargs)
        extra = {'index': 'num_columns', 'name': 'num_columns', 'label': 'Number of columns', 'sorttype': 'int'}
        colmodels.append(extra)
        return colmodels


class AbstractVariantGrid(JqGridSQL):
    model = Variant
    INTERNAL_CLASSIFICATION_ALIASES_AND_SELECT = {
        "internally_classified": "select string_agg(classification_variantclassification.clinical_significance, '|') from classification_variantclassification where classification_variantclassification.variant_id = snpdb_variant.id",
        "max_internal_classification": "select max(classification_variantclassification.clinical_significance) from classification_variantclassification where classification_variantclassification.variant_id = snpdb_variant.id",
    }

    def column_in_queryset_fields(self, field):
        colmodel = self.get_override(field)
        return colmodel.get("queryset_field", True)

    def get_queryset_field_names(self):
        field_names = []
        for f in super().get_field_names():
            if self.column_in_queryset_fields(f):
                field_names.append(f)

        return field_names

    def get_sql_params_and_columns(self, request):
        queryset = self.filter_items(request, self.queryset)

        sidx = request.GET.get('sidx', 'id')
        if self.column_in_queryset_fields(sidx):
            queryset = self.sort_items(request, queryset)

        (select_part, from_part, where_part) = get_queryset_select_from_where_parts(queryset)

        extra_columns = []
        for (alias, select) in self.INTERNAL_CLASSIFICATION_ALIASES_AND_SELECT.items():
            extra_columns.append(f'({select}) as "{alias}"')
        select_part += ",\n" + ",\n".join(extra_columns)

        sql = '\n'.join([select_part, from_part, where_part])
        #logging.info(sql)

        column_names = get_queryset_column_names(queryset, list(self.INTERNAL_CLASSIFICATION_ALIASES_AND_SELECT.keys()))

        params = []
        return sql, params, column_names, True

    def get_count(self):
        return self.queryset.count()
