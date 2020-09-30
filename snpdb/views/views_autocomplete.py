from abc import ABC

from django.contrib.auth.models import User
from django.db.models.functions import Length
from django.db.models.query_utils import Q
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from library.constants import MINUTE_SECS
from library.django_utils.autocomplete_utils import AutocompleteView
from snpdb.models import VCF, Sample, Cohort, CustomColumnsCollection, CustomColumn, Tag, Trio, \
    Lab, GenomicIntervalsCollection, GenomeBuild, ImportStatus


class GenomeBuildAutocompleteView(AutocompleteView, ABC):

    def filter_to_genome_build(self, qs, path_to_genome_build):
        genome_build_id = self.forwarded.get('genome_build_id')
        if genome_build_id:
            genome_build = GenomeBuild.objects.get(pk=genome_build_id)
            qs = qs.filter(**{path_to_genome_build: genome_build})
        return qs


@method_decorator(cache_page(MINUTE_SECS), name='dispatch')
class UserAutocompleteView(AutocompleteView):
    fields = ['last_name', 'first_name', 'username']

    def get_user_queryset(self, user):
        return User.objects.all()


@method_decorator(cache_page(MINUTE_SECS), name='dispatch')
class UsernameAutocompleteView(AutocompleteView):
    """
    Needed to make separate from UserAutocompleteView for the sake of sort order
    """
    fields = ['username', 'first_name', 'last_name']

    def get_user_queryset(self, user):
        return User.objects.all()

    def get_result_label(self, obj):
        return obj.username


@method_decorator(cache_page(MINUTE_SECS), name='dispatch')
class LabAutocompleteView(AutocompleteView):
    fields = ['organization__name', 'name']

    def get_user_queryset(self, user):
        return Lab.objects.filter(organization__active=True)

    def get_result_label(self, obj):
        return f'{obj.organization.name} - {obj.name}'


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class VCFAutocompleteView(GenomeBuildAutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        qs = VCF.filter_for_user(user, True).filter(import_status=ImportStatus.SUCCESS)
        return self.filter_to_genome_build(qs, "genome_build")


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class SampleAutocompleteView(GenomeBuildAutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        sample_qs = Sample.filter_for_user(user, True).filter(import_status=ImportStatus.SUCCESS)
        return self.filter_to_genome_build(sample_qs, "vcf__genome_build")


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class CohortAutocompleteView(GenomeBuildAutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        vcf_success_if_exists = Q(vcf__isnull=True) | Q(vcf__import_status=ImportStatus.SUCCESS)
        qs = Cohort.filter_for_user(user, success_status_only=True).filter(vcf_success_if_exists)
        return self.filter_to_genome_build(qs, "genome_build")


class CustomColumnAutocompleteView(AutocompleteView):
    fields = ['column__grid_column_name']

    def get_user_queryset(self, user):
        custom_columns_collections_qs = CustomColumnsCollection.filter_for_user(user)
        # Called different things in Analysis/UserSettings
        columns = self.forwarded.get('columns') or self.forwarded.get('custom_columns_collection')
        if columns:
            custom_columns_collections_qs = custom_columns_collections_qs.filter(pk=columns)

        return CustomColumn.objects.filter(custom_columns_collection__in=custom_columns_collections_qs)


class GenomicIntervalsCollectionAutocompleteView(GenomeBuildAutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        qs = GenomicIntervalsCollection.filter_for_user(user).filter(import_status=ImportStatus.SUCCESS)
        return self.filter_to_genome_build(qs, "genome_build")


@method_decorator(cache_page(5), name='dispatch')
class TagAutocompleteView(AutocompleteView):
    fields = ['id']

    def get_user_queryset(self, user):
        return Tag.objects.all().order_by(Length("id").asc())


@method_decorator([cache_page(MINUTE_SECS), vary_on_cookie], name='dispatch')
class TrioAutocompleteView(GenomeBuildAutocompleteView):
    fields = ['name']

    def get_user_queryset(self, user):
        qs = Trio.filter_for_user(user, success_status_only=True)
        return self.filter_to_genome_build(qs, "cohort__genome_build")
