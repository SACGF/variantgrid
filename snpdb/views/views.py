import itertools
import json
import logging
from collections import OrderedDict, defaultdict
from typing import Iterable, Tuple

import pandas as pd
from celery.result import AsyncResult
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import User, Group
from django.core.exceptions import PermissionDenied, ImproperlyConfigured, ObjectDoesNotExist
from django.db.utils import IntegrityError
from django.forms.models import inlineformset_factory, ALL_FIELDS
from django.forms.widgets import TextInput
from django.http.response import HttpResponse, HttpResponseRedirect, HttpResponseServerError, JsonResponse
from django.shortcuts import get_object_or_404, render, redirect
from django.urls.base import reverse
from django.views.decorators.cache import cache_page
from django.views.decorators.http import require_POST
from django.views.decorators.vary import vary_on_cookie
from django_messages.models import Message
from global_login_required import login_not_required
from guardian.shortcuts import get_objects_for_group, get_objects_for_user
from termsandconditions.decorators import terms_required

from analysis.analysis_templates import get_sample_analysis
from analysis.forms import AnalysisOutputNodeChoiceForm
from analysis.models import AnalysisTemplate
from annotation.forms import GeneCountTypeChoiceForm
from annotation.manual_variant_entry import create_manual_variants, can_create_variants
from annotation.models import AnnotationVersion, SampleVariantAnnotationStats, SampleGeneAnnotationStats, \
    SampleClinVarAnnotationStats, SampleVariantAnnotationStatsPassingFilter, SampleGeneAnnotationStatsPassingFilter, \
    SampleClinVarAnnotationStatsPassingFilter
from annotation.models.models import ManualVariantEntryCollection, VariantAnnotationVersion
from annotation.models.models_gene_counts import GeneValueCountCollection, \
    GeneCountType, SampleAnnotationVersionVariantSource, CohortGeneCounts
from annotation.serializers import ManualVariantEntryCollectionSerializer
from classification.classification_stats import get_grouped_classification_counts
from classification.models.clinvar_export_sync import clinvar_export_sync
from classification.views.classification_accumulation_graph import get_accumulation_graph_data, \
    AccumulationReportMode
from classification.views.classification_datatables import ClassificationColumns
from genes.custom_text_gene_list import create_custom_text_gene_list
from genes.forms import CustomGeneListForm, UserGeneListForm, GeneAndTranscriptForm
from genes.models import GeneListCategory, CustomTextGeneList, GeneList
from library.constants import WEEK_SECS, HOUR_SECS
from library.django_utils import add_save_message, get_model_fields, set_form_read_only
from library.guardian_utils import DjangoPermission
from library.keycloak import Keycloak
from library.utils import full_class_name, import_class, rgb_invert
from ontology.models import OntologyTerm
from patients.forms import PatientForm
from patients.models import Patient, Clinician
from patients.views import get_patient_upload_csv
from snpdb import forms
from snpdb.sample_file_path import get_example_replacements
from snpdb.forms import SampleChoiceForm, VCFChoiceForm, \
    UserSettingsOverrideForm, UserForm, UserContactForm, SampleForm, TagForm, SettingsInitialGroupPermissionForm, \
    OrganizationForm, LabForm, LabUserSettingsOverrideForm, OrganizationUserSettingsOverrideForm
from snpdb.graphs import graphcache
from snpdb.graphs.allele_frequency_graph import AlleleFrequencyHistogramGraph
from snpdb.graphs.chromosome_density_graph import SampleChromosomeDensityGraph
from snpdb.graphs.chromosome_intervals_graph import ChromosomeIntervalsGraph
from snpdb.graphs.homozygosity_percent_graph import HomozygosityPercentGraph
from snpdb.import_status import set_vcf_and_samples_import_status
from snpdb.models import CachedGeneratedFile, VariantGridColumn, UserSettings, \
    VCF, UserTagColors, CustomColumnsCollection, CustomColumn, Cohort, \
    CohortSample, GenomicIntervalsCollection, Sample, UserDataPrefix, UserGridConfig, \
    get_igv_data, SampleLocusCount, UserContact, Tag, Wiki, Organization, GenomeBuild, \
    Trio, AbstractNodeCountSettings, CohortGenotypeCollection, UserSettingsOverride, NodeCountSettingsCollection, Lab, \
    LabUserSettingsOverride, OrganizationUserSettingsOverride, LabHead, SomalierRelatePairs, \
    VariantZygosityCountCollection, VariantZygosityCountForVCF, ClinVarKey, AvatarDetails, State, SampleStats, \
    SampleStatsPassingFilter
from snpdb.models.models_enums import ProcessingStatus, ImportStatus, BuiltInFilters
from snpdb.tasks.soft_delete_tasks import soft_delete_vcfs
from snpdb.utils import LabNotificationBuilder
from upload.models import UploadedVCF
from upload.uploaded_file_type import retry_upload_pipeline


@terms_required
def index(request):
    if Clinician.user_is_clinician(request.user):
        return redirect('clinician_login')

    return render(request, 'index.html')


def data(request):
    return render(request, 'snpdb/data/data.html')


def maps(request):
    return render(request, 'maps.html')


def get_writable_class_object(user, class_name, primary_key):
    klass = import_class(class_name)
    name = klass.__name__
    obj = klass.objects.get(pk=primary_key)

    if not obj.can_write(user):
        write_perm = DjangoPermission.perm(obj, DjangoPermission.WRITE)
        msg = f"You do not have permission {write_perm} needed to modify {name}"
        raise PermissionDenied(msg)

    return obj, name


def get_writable_class_objects(user, class_name):
    klass = import_class(class_name)
    name = klass.__name__
    write_perm = DjangoPermission.perm(klass, DjangoPermission.WRITE)
    qs = get_objects_for_user(user, write_perm, klass=klass, accept_global_perms=False)
    return qs, name


def group_permissions(request, class_name, primary_key):
    obj, name = get_writable_class_object(request.user, class_name, primary_key)

    try:
        # If object has "get_permission_object" it can delegate it.
        permission_obj = obj.get_permission_object()
        perm_obj_name = permission_obj.__class__.__name__
    except AttributeError:
        # Default is use itself
        permission_obj = obj
        perm_obj_name = name

    permission_forms = get_group_permission_forms(request, permission_obj)

    if request.method == 'POST':
        valid = all([pf.is_valid() for pf in permission_forms])
        if valid:
            for pf in permission_forms:
                pf.save()
        add_save_message(request, valid, f"{perm_obj_name} group permissions")

    get_listing_url = getattr(obj, "get_listing_url", None)
    if get_listing_url:
        delete_redirect_url = get_listing_url()
    else:
        delete_redirect_url = "/"

    context = {'permission_forms': permission_forms,
               'class_name': class_name,
               'name': name,
               'perm_obj_name': perm_obj_name,
               'permission_obj': permission_obj,
               'instance': obj,
               'delete_redirect_url': delete_redirect_url}
    return render(request, 'snpdb/data/group_permissions.html', context)


@require_POST
def group_permissions_object_delete(request, class_name, primary_key):
    if class_name == 'snpdb.models.VCF':  # TODO: Hack? Make some class object?
        soft_delete_vcfs(request.user, primary_key)
    else:
        obj, _ = get_writable_class_object(request.user, class_name, primary_key)
        try:
            obj.delete()
        except IntegrityError as ie:
            pks = ", ".join(str(o.pk) for o in ie.args[1])
            error_message = f"{ie.args[0]}: {pks}"
            return HttpResponseServerError(content=error_message)

    return HttpResponse()


def bulk_group_permissions(request, class_name):
    qs, name = get_writable_class_objects(request.user, class_name)
    groups = list(request.user.groups.all().order_by("name"))

    objects_and_forms = []
    for obj in qs:
        permission_forms = get_group_permission_forms(request, obj, groups=groups)
        objects_and_forms.append((obj, permission_forms))

    if request.method == 'POST':
        all_forms = []
        for _, permission_forms in objects_and_forms:
            all_forms.extend(permission_forms)
        valid = all([pf.is_valid() for pf in all_forms])
        if valid:
            for pf in all_forms:
                pf.save()
        add_save_message(request, valid, f"{name} group permissions")

    context = {"name": name,
               "groups": groups,
               "objects_and_forms": objects_and_forms}
    return render(request, 'snpdb/data/bulk_group_permissions.html', context)


def _get_vcf_sample_stats(vcf, klass):
    """ Count is het + hom """
    ss_fields = ("sample_id", "sample__name", "variant_count", "ref_count", "het_count", "hom_count", "unk_count")
    ss_values_qs = klass.objects.filter(sample__vcf=vcf).order_by("sample").values(*ss_fields)

    sample_stats_het_hom_count = {}
    sample_names = []
    sample_zygosities = defaultdict(list)
    for value_dict in ss_values_qs:
        sample_id = value_dict.pop("sample_id")
        sample_names.append(value_dict.pop("sample__name"))
        value_dict.pop("variant_count")
        sample_stats_het_hom_count[sample_id] = value_dict["het_count"] + value_dict["hom_count"]
        for k, v in value_dict.items():
            sample_zygosities[k].append(v)

    return sample_stats_het_hom_count, sample_names, tuple(sample_zygosities.items())


def view_vcf(request, vcf_id):
    vcf = VCF.get_for_user(request.user, vcf_id)
    # I couldn't get prefetch_related_objects([vcf], "sample_set__samplestats") to work - so storing in a dict

    sample_stats_het_hom_count, sample_names, sample_zygosities = _get_vcf_sample_stats(vcf, SampleStats)
    sample_stats_pass_het_hom_count, _, sample_zygosities_pass = _get_vcf_sample_stats(vcf, SampleStatsPassingFilter)

    VCFSampleFormSet = inlineformset_factory(VCF, Sample, extra=0, can_delete=False,
                                             fields=["vcf_sample_name", "name", "patient", "specimen"],
                                             widgets=SampleForm.Meta.widgets)

    post = request.POST or None
    vcf_form = forms.VCFForm(post, instance=vcf)
    samples_form = VCFSampleFormSet(post, instance=vcf)
    for form in samples_form.forms:
        form.fields["vcf_sample_name"].disabled = True

    requires_user_input = vcf.import_status == ImportStatus.REQUIRES_USER_INPUT
    reload_vcf = False
    if request.method == 'POST':
        valid = all(f.is_valid() for f in [vcf_form, samples_form])
        if valid:
            vcf = vcf_form.save()
            reload_vcf = requires_user_input and vcf.genome_build

            samples_form.save()

        add_save_message(request, valid, "VCF")

    try:
        # Some legacy data was too hard to fix and relies on being re-imported
        _ = vcf.cohort
        _ = vcf.cohort.cohort_genotype_collection
    except (Cohort.DoesNotExist, CohortGenotypeCollection.DoesNotExist):
        messages.add_message(request, messages.ERROR, "This legacy VCF is missing data and needs to be reloaded.")

    if reload_vcf:
        set_vcf_and_samples_import_status(vcf, ImportStatus.IMPORTING)
        retry_upload_pipeline(vcf.uploadedvcf.uploaded_file.uploadpipeline)
        vcf_form = forms.VCFForm(post, instance=vcf)  # Reload as import status has changed
        messages.add_message(request, messages.INFO, "Reloading VCF")

    for warning, _ in vcf.get_warnings():
        messages.add_message(request, messages.WARNING, warning, extra_tags='import-message')

    has_write_permission = vcf.can_write(request.user)
    if not has_write_permission:
        messages.add_message(request, messages.WARNING, "You can view but not modify this data.")

    variant_zygosity_count_collections = {}
    for vzcc in VariantZygosityCountCollection.objects.all():
        vzc_vcf = VariantZygosityCountForVCF.objects.filter(vcf=vcf, collection=vzcc).first()
        variant_zygosity_count_collections[vzcc] = vzc_vcf

    try:
        can_view_upload_pipeline = vcf.uploadedvcf.uploaded_file.can_view(request.user)
    except UploadedVCF.DoesNotExist:
        can_view_upload_pipeline = False

    context = {
        'vcf': vcf,
        'sample_stats_het_hom_count': sample_stats_het_hom_count,
        'sample_stats_pass_het_hom_count': sample_stats_pass_het_hom_count,
        'sample_names': sample_names,
        'sample_zygosities': sample_zygosities,
        'vcf_form': vcf_form,
        'samples_form': samples_form,
        'patient_form': PatientForm(user=request.user),  # blank
        'has_write_permission': has_write_permission,
        'can_download_vcf': (not settings.VCF_DOWNLOAD_ADMIN_ONLY) or request.user.is_superuser,
        'can_view_upload_pipeline': can_view_upload_pipeline,
        "variant_zygosity_count_collections": variant_zygosity_count_collections,
    }
    return render(request, 'snpdb/data/view_vcf.html', context)


def get_patient_upload_csv_for_vcf(request, pk):
    vcf = VCF.get_for_user(request.user, pk)
    sample_qs = vcf.sample_set.all()
    filename = f"vcf_{pk}_patient_upload"
    return get_patient_upload_csv(filename, sample_qs)


def _sample_stats(sample) -> Tuple[pd.DataFrame, pd.DataFrame]:
    annotation_version = AnnotationVersion.latest(sample.genome_build)
    STATS = {
        "Total": (SampleStats, SampleStatsPassingFilter, set()),
        "dbSNP": (SampleVariantAnnotationStats, SampleVariantAnnotationStatsPassingFilter, {"variant_annotation_version"}),
        "OMIM pheno": (SampleGeneAnnotationStats, SampleGeneAnnotationStatsPassingFilter, {"gene_annotation_version"}),
        "ClinVar LP/P": (SampleClinVarAnnotationStats, SampleClinVarAnnotationStatsPassingFilter, {"clinvar_version"})
    }

    VARIANT_CLASS = ["variant", "snp", "insertions", "deletions"]
    ZYGOSITY = ["ref", "het", "hom", "unk"]

    variant_class_data = {}
    zygosity_data = {}
    for name, (stats_klass, stats_passing_filters_klass, shared_fields) in STATS.items():
        kwargs = {"sample": sample}
        for sf in shared_fields:
            kwargs[sf] = getattr(annotation_version, sf)

        objs = {}
        try:
            objs[name] = stats_klass.objects.get(**kwargs)
        except ObjectDoesNotExist:
            pass

        try:
            objs[f"{name} PASS filters"] = stats_passing_filters_klass.objects.get(**kwargs)
        except ObjectDoesNotExist:
            pass

        for n, o in objs.items():
            obj_variant_class_data = {}
            for field in get_model_fields(o):
                for k in VARIANT_CLASS:
                    if field.startswith(k):
                        obj_variant_class_data[k] = getattr(o, field)
                        break
            if obj_variant_class_data:
                variant_class_data[n] = obj_variant_class_data

            obj_zygosity_data = {}
            for field in get_model_fields(o):
                for k in ZYGOSITY:
                    if field.startswith(k):
                        obj_zygosity_data[k] = getattr(o, field)
                        break
            if obj_zygosity_data:
                zygosity_data[n] = obj_zygosity_data

    sample_stats_variant_class_df = pd.DataFrame.from_dict(variant_class_data).reindex(VARIANT_CLASS)
    if "Total" in sample_stats_variant_class_df.columns:
        total = sample_stats_variant_class_df["Total"]
        sample_stats_variant_class_df["Total %"] = 100 * total / total["variant"]

    sample_stats_zygosity_df = pd.DataFrame.from_dict(zygosity_data).reindex(ZYGOSITY)
    return sample_stats_variant_class_df, sample_stats_zygosity_df


def view_sample(request, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)
    has_write_permission = sample.can_write(request.user)

    form = forms.SampleForm(request.POST or None, instance=sample)
    if not has_write_permission:
        set_form_read_only(form)
        messages.add_message(request, messages.WARNING, "You can view but not modify this data.")

    if request.method == 'POST':
        if not has_write_permission:
            raise PermissionDenied("Can't modify public data")

        valid = form.is_valid()
        if valid:
            form.save()
        add_save_message(request, valid, "Sample")

    sample_locus_count = list(SampleLocusCount.objects.filter(sample=sample).order_by("locus_count"))
    igv_data = get_igv_data(request.user, genome_build=sample.genome_build)
    patient_form = PatientForm(user=request.user)  # blank
    related_samples = None
    if settings.SOMALIER.get("enabled"):
        related_samples = SomalierRelatePairs.get_for_sample(sample).order_by("relate")

    sample_stats_variant_class_df, sample_stats_zygosity_df = _sample_stats(sample)
    context = {
        'sample': sample,
        'samples': [sample],
        'sample_locus_count': sample_locus_count,
        'form': form,
        'patient_form': patient_form,
        'cohorts': cohorts,
        'has_write_permission': has_write_permission,
        'igv_data': igv_data,
        "bam_list": sample.get_bam_files(),
        "sample_stats_variant_class_df": sample_stats_variant_class_df,
        "sample_stats_zygosity_df": sample_stats_zygosity_df,
        "related_samples": related_samples
    }
    return render(request, 'snpdb/data/view_sample.html', context)


def sample_files_tab(request, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)
    if request.method == "POST":
        sample_files_formset = forms.SampleFilesFormSet(request.POST, instance=sample)
        valid = sample_files_formset.is_valid()
        if valid:
            sample_files_formset.save()
        add_save_message(request, valid, "Sample Files")

    # We shouldn't re-use after POST - so generate fresh
    sample_files_formset = forms.SampleFilesFormSet(None, instance=sample)

    context = {
        "sample": sample,
        "sample_files_formset": sample_files_formset,
        'has_write_permission': sample.can_write(request.user),
    }
    return render(request, 'snpdb/data/sample_files_tab.html', context)


def sample_variants_tab(request, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)
    analysis = None
    error_message = None
    if settings.ANALYSIS_TEMPLATES_AUTO_SAMPLE:
        try:
            analysis_template = AnalysisTemplate.objects.get(name=settings.ANALYSIS_TEMPLATES_AUTO_SAMPLE)
            analysis = get_sample_analysis(sample, analysis_template)
        except AnalysisTemplate.DoesNotExist:
            error_message = f"Analysis Template '{settings.ANALYSIS_TEMPLATES_AUTO_SAMPLE}' does not exist!"
    else:
        error_message = "settings.ANALYSIS_TEMPLATES_AUTO_SAMPLE not set. Talk to your administrator"

    if error_message:
        messages.add_message(request, messages.ERROR, error_message)

    context = {
        'sample': sample,
        "analysis": analysis,
        'output_node_form': AnalysisOutputNodeChoiceForm(analysis=analysis)
    }
    return render(request, 'snpdb/data/sample_variants_tab.html', context)


def sample_variants_gene_detail(request, sample_id, gene_symbol):
    sample = Sample.get_for_user(request.user, sample_id)

    context = {'sample': sample,
               'sample_ids': [sample.pk],
               'gene_symbol': gene_symbol,
               "datatable_config": ClassificationColumns(request)}
    return render(request, 'snpdb/data/sample_variants_gene_detail.html', context)


def sample_graphs_tab(request, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)

    context = {'sample': sample}
    return render(request, 'snpdb/data/sample_graphs_tab.html', context)


def get_group_permission_forms(request, obj, groups=None):
    if groups is None:
        groups = request.user.groups.all().order_by("name")
    return [forms.GroupPermissionForm(request.POST or None, obj=obj, group=group) for group in groups]


def sample_permissions_tab(request, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)

    context = {'sample': sample,
               'class_name': full_class_name(Sample)}
    return render(request, 'snpdb/data/sample_permissions_tab.html', context)


def view_genomic_intervals(request, genomic_intervals_collection_id):
    gic = get_object_or_404(GenomicIntervalsCollection, pk=genomic_intervals_collection_id)
    if not request.user.has_perm('view_genomicintervalscollection', gic):
        raise PermissionDenied()

    form = forms.GenomicIntervalsCollectionForm(request.POST or None, instance=gic)
    if request.method == "POST":
        valid = form.is_valid()
        if valid:
            gic = form.save()
        add_save_message(request, valid, "Genomic Intervals")

    if gic.genome_build is None:
        msg = "Unable to automatically set build, please select manually."
        messages.add_message(request, messages.WARNING, msg, extra_tags='import-message')

    context = {'gic': gic,
               'form': form,
               "has_write_permission": gic.can_write(request.user)}
    return render(request, 'snpdb/data/view_genomic_intervals.html', context)


@require_POST
def cached_generated_file_delete(request):
    cgf_id = request.POST["cgf_id"]
    cgf = get_object_or_404(CachedGeneratedFile, pk=cgf_id)
    cgf.delete()
    return HttpResponse()


def vcfs(request):
    context = {
        "form": VCFChoiceForm(),
    }
    return render(request, 'snpdb/data/vcfs.html', context=context)


def samples(request):
    groups = request.user.groups.values_list("name", flat=True)
    groups_str = ', '.join(groups)
    num_groups = len(groups)
    if num_groups > 1:
        group_info = f"(or owned by one of your groups: {groups_str})"
    elif num_groups:
        group_info = f"(or owned by your group: {groups_str})"
    else:
        group_info = ''

    context = {
        "form": SampleChoiceForm(),
        "group_info": group_info,
    }
    return render(request, 'snpdb/data/samples.html', context=context)


def bed_files(request):
    return render(request, 'snpdb/data/bed_files.html')


@require_POST
def messages_bulk_delete(request):
    messages_str = request.POST['message_ids']
    message_ids = json.loads(messages_str)
    user_messages_qs = Message.objects.filter(recipient=request.user)
    user_messages_qs.filter(pk__in=message_ids).delete()
    return HttpResponse()


def manual_variant_entry(request):

    if can_create_variants(request.user):
        form = forms.ManualVariantEntryForm(request.POST or None, user=request.user)
        if request.method == 'POST':
            valid = form.is_valid()
            if valid:
                variants_text = form.cleaned_data['variants_text']
                genome_build_pk = form.cleaned_data['genome_build']
                genome_build = GenomeBuild.objects.get(pk=genome_build_pk)
                create_manual_variants(request.user, genome_build, variants_text)
                form = forms.ManualVariantEntryForm(None, user=request.user)  # Reset form

            add_save_message(request, valid, "Manually entered variants")
    else:
        form = None
        messages.add_message(request, messages.INFO, "Manual variant entry has been disabled by an admin.")

    mvec_qs = ManualVariantEntryCollection.objects.order_by("-id")
    context = {"form": form,
               "mvec_qs": mvec_qs}
    return render(request, 'snpdb/data/manual_variant_entry.html', context=context)


def watch_manual_variant_entry(request, pk):
    mvec = ManualVariantEntryCollection.get_for_user(request.user, pk)
    # TODO: Quick redirect to variant if it's already ready

    mvec_data = ManualVariantEntryCollectionSerializer(mvec).data
    context = {"mvec": mvec,
               "initial_json": json.dumps(mvec_data)}
    return render(request, 'snpdb/data/watch_manual_variant_entry.html', context=context)


@require_POST
def set_user_row_config(request):
    """ This is set from jqgrid.html setRowChangeCallbacks when changing grid rows """

    grid_name = request.POST["grid_name"]
    grid_rows = int(request.POST["grid_rows"])

    UserGridConfig.objects.update_or_create(user=request.user, grid_name=grid_name, defaults={"rows": grid_rows})
    return HttpResponse()


@require_POST
def set_user_data_grid_config(request):
    """ This is set from user_data_grid_filter.html, should contain either filter_level+checked or filter_name """

    grid_name = request.POST["grid_name"]
    user_grid_config = UserGridConfig.get(request.user, grid_name)
    filter_level = request.POST.get("filter_level")
    if filter_level:
        checked = json.loads(request.POST["checked"])

        if filter_level == 'groups':
            user_grid_config.show_group_data = checked
        elif filter_level == 'incomplete':
            user_grid_config.show_incomplete_data = checked
        elif filter_level == 'hidden':
            user_grid_config.show_hidden_data = checked
        else:
            msg = f"Unknown value for filter_level: '{filter_level}'"
            raise ValueError(msg)
    else:
        user_grid_config.filter_name = request.POST["filter_name"]

    user_grid_config.save()
    return HttpResponse()


def view_user_settings(request):
    user = request.user
    user_contact = UserContact.get_for_user(user)

    action = request.POST.get('action') if request.POST else None
    post = request.POST or None if not action else None
    user_form = UserForm(post, instance=user)
    user_contact_form = UserContactForm(post, instance=user_contact)
    user_settings = UserSettings.get_for_user(user)
    override_source, override_values = user_settings.get_override_source_and_values_before_user()
    user_settings_override = UserSettingsOverride.objects.get(user=user)
    user_settings_override_form = UserSettingsOverrideForm(post, instance=user_settings_override)
    labs_by_group_name = {l.group_name: l for l in Lab.valid_labs_qs(user)}

    group_initial_perm_forms = {}
    if settings.USER_SETTINGS_SHOW_GROUPS:
        read_groups, write_groups = user_settings.initial_perm_read_and_write_groups
        for group in user.groups.all().order_by('name'):
            initial = {"read": group in read_groups, "write": group in write_groups}
            group_initial_perm_forms[group] = SettingsInitialGroupPermissionForm(request.POST or None, initial=initial,
                                                                                 settings_override=user_settings_override,
                                                                                 group=group)
    if request.method == "POST":
        all_valid = True

        action = request.POST.get('action')
        if action == 'password-reset':
            keycloak = Keycloak()
            keycloak.change_password(user)
            messages.add_message(request, level=messages.INFO, message='Password reset email sent',
                                 extra_tags='save-message')
        else:
            if not settings.USE_OIDC:
                if user_form.is_valid():
                    user = user_form.save()
                else:
                    all_valid = False

            for form in itertools.chain([user_contact_form, user_settings_override_form],
                                        group_initial_perm_forms.values()):
                if form.is_valid():
                    form.save()
                else:
                    all_valid = False
            add_save_message(request, all_valid, "User Settings")

    context = {
        'user': user,
        'user_form': user_form,
        'user_contact_form': user_contact_form,
        'user_settings_form': user_settings_override_form,
        'group_initial_perm_forms': group_initial_perm_forms,
        'accounts_email': settings.ACCOUNTS_EMAIL,
        'account_manage_url': settings.OIDC_USER_SERVICES,
        'override_source': override_source,
        'override_values': override_values,
        'labs_by_group_name': labs_by_group_name,
        'avatar_details': AvatarDetails.avatar_for(user)
    }
    return render(request, 'snpdb/settings/view_user_settings.html', context)


def user_settings_node_counts_tab(request):
    user_settings_override = UserSettingsOverride.objects.get_or_create(user=request.user)[0]
    return _settings_override_node_counts_tab(request, user_settings_override)


def lab_settings_node_counts_tab(request, pk):
    lab = get_object_or_404(Lab, pk=pk)
    has_write_permission = lab.can_write(request.user)
    if has_write_permission is False:
        _add_read_only_settings_message(request, [lab])
    lab_settings_override = LabUserSettingsOverride.objects.get_or_create(lab=lab)[0]
    return _settings_override_node_counts_tab(request, lab_settings_override, has_write_permission=has_write_permission)


def organization_settings_node_counts_tab(request, pk):
    organization = get_object_or_404(Organization, pk=pk)
    has_write_permission = organization.can_write(request.user)
    if has_write_permission is False:
        _add_read_only_settings_message(request, organization.lab_set.all())

    org_settings_override = OrganizationUserSettingsOverride.objects.get_or_create(organization=organization)[0]
    return _settings_override_node_counts_tab(request, org_settings_override, has_write_permission=has_write_permission)


def _settings_override_node_counts_tab(request, settings_override, has_write_permission=True):
    # This calls _analysis_settings_node_counts_tab with a FakeAnalysis object that
    # handles loading/saving a global one against User settings objects instead of analysis
    class FakeAnalysis:

        def set_node_count_types(self, node_counts_array):
            collection, _ = NodeCountSettingsCollection.objects.get_or_create(settings=settings_override)
            AbstractNodeCountSettings.save_count_configs_from_array(collection.nodecountsettings_set, node_counts_array)

        def get_node_count_types(self):
            try:
                node_count_config = settings_override.nodecountsettingscollection
                node_count_filters = node_count_config.get_node_count_filters()
            except:
                node_count_filters = BuiltInFilters.DEFAULT_NODE_COUNT_FILTERS

            return AbstractNodeCountSettings.get_types_from_labels(node_count_filters)

    fake_analysis = FakeAnalysis()
    from analysis.views.views import _analysis_settings_node_counts_tab  # Circular import
    return _analysis_settings_node_counts_tab(request, fake_analysis,
                                              pass_analysis_settings=False, has_write_permission=has_write_permission)


def view_user(request, pk):
    user = get_object_or_404(User, pk=pk)
    user_contact = UserContact.get_for_user(user)

    context = {"user": user,
               'user_contact': user_contact}
    return render(request, 'snpdb/settings/view_user.html', context)


def _add_read_only_settings_message(request, lab_list: Iterable[Lab]):
    """ lab_list: labs where lab heads can modify settings """
    lab_heads_qs = LabHead.objects.filter(lab__in=lab_list).distinct()
    lab_head_names = ", ".join([str(lh.user) for lh in lab_heads_qs])
    if lab_head_names:
        lab_head_msg = f" or lab heads: {lab_head_names}"
    else:
        lab_head_msg = ""
    read_only_message = f"Only administrators{lab_head_msg} can modify these settings"
    messages.add_message(request, messages.INFO, read_only_message)


def view_lab(request, lab_id: int):
    lab = get_object_or_404(Lab, pk=lab_id)

    lab_form = LabForm(request.POST or None, instance=lab)
    lab_settings_override = LabUserSettingsOverride.objects.get_or_create(lab=lab)[0]

    override_fields = set(get_model_fields(LabUserSettingsOverride)) - {"id", "settingsoverride_ptr", "lab"}
    parent_overrides = UserSettings.get_settings_overrides(organization=lab.organization)
    override_source, override_values = UserSettings.get_override_source_and_values(override_fields, parent_overrides)
    settings_overrides = parent_overrides + [lab_settings_override]
    read_groups, write_groups = UserSettings.get_initial_perm_read_and_write_groups([lab.group], settings_overrides)

    initial = {"read": lab.group in read_groups, "write": lab.group in write_groups}
    group_initial_perm_form = None
    if settings.USER_SETTINGS_SHOW_GROUPS:
        group_initial_perm_form = SettingsInitialGroupPermissionForm(request.POST or None, initial=initial,
                                                                 settings_override=lab_settings_override,
                                                                 group=lab.group)
    lab_settings_override_form = LabUserSettingsOverrideForm(request.POST or None, instance=lab_settings_override)

    has_write_permission = lab.can_write(request.user)
    all_forms = [form for form in [lab_form, group_initial_perm_form, lab_settings_override_form] if form]

    if request.method == "POST":
        lab.check_can_write(request.user)

        if debug_method := request.POST.get("debug_method"):
            if "Test Slack" == debug_method:
                if not lab.slack_webhook:
                    messages.add_message(request, messages.ERROR, "Slack URL not configured correctly")
                else:
                    #try:
                    notification_builder = LabNotificationBuilder(lab=lab, message="Testing Slack Integration", notification_type=LabNotificationBuilder.NotificationType.SLACK_ONLY)
                    notification_builder.add_header(f"{settings.SITE_NAME} -> Slack Integration Test")
                    notification_builder.add_markdown("If you can see this, then integration has worked! :smile:")
                    notification_builder.send()
                    messages.add_message(request, messages.SUCCESS, "Message sent, check your Slack to confirm")
                    #except:
                    #    report_exc_info()
                    #    messages.add_message(request, messages.ERROR, "Unable to send test notification")
                return redirect(reverse('view_lab', kwargs={"lab_id": lab_id}))
            else:
                raise ValueError(f"Un-supported debug method {debug_method}")
        else:
            all_valid = True
            for form in all_forms:
                if form.is_valid():
                    form.save()
                else:
                    all_valid = False
            add_save_message(request, all_valid, "Lab Settings")

    if has_write_permission is False:
        for form in all_forms:
            set_form_read_only(form)
        # we just hide the form now
        # _add_read_only_settings_message(request, [lab])

    if settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED:
        visibility = "Shared"
    else:
        visibility = f"Created"

    context = {
        "lab": lab,
        "visibility": visibility,
        "is_member": lab.is_member(request.user) or request.user.is_superuser,

        "lab_form": lab_form,
        'settings_override_form': lab_settings_override_form,
        'group_initial_perm_form': group_initial_perm_form,
        'override_source': override_source,
        'override_values': override_values,
        'has_write_permission': has_write_permission,
        'clinvar_export_enabled': clinvar_export_sync.is_enabled
    }
    return render(request, 'snpdb/settings/view_lab.html', context)


def view_clinvar_key(request, pk: str):
    clinvar_key = get_object_or_404(ClinVarKey, pk=pk)
    clinvar_key.check_user_can_access(request.user)

    return render(request, 'snpdb/settings/clinvar_key.html', {
        'clinvar_key': clinvar_key,
        'labs': Lab.objects.filter(clinvar_key=clinvar_key).order_by('name')
    })


def view_organization(request, organization_id: int):
    organization = get_object_or_404(Organization, pk=organization_id)
    organization_form = OrganizationForm(request.POST or None, instance=organization)
    org_settings_override = OrganizationUserSettingsOverride.objects.get_or_create(organization=organization)[0]
    override_fields = set(get_model_fields(OrganizationUserSettingsOverride)) - {"id", "settingsoverride_ptr", "organization"}
    parent_overrides = UserSettings.get_settings_overrides()
    override_source, override_values = UserSettings.get_override_source_and_values(override_fields, parent_overrides)
    org_settings_override_form = OrganizationUserSettingsOverrideForm(request.POST or None, instance=org_settings_override)
    all_forms = [organization_form, org_settings_override_form]

    if request.method == "POST":
        organization.check_can_write(request.user)
        all_valid = True
        for form in all_forms:
            if form.is_valid():
                form.save()
            else:
                all_valid = False
        add_save_message(request, all_valid, "Organization Settings")

    has_write_permission = organization.can_write(request.user)
    if has_write_permission is False:
        for form in all_forms:
            set_form_read_only(form)
        # put on individual tabs now
        # _add_read_only_settings_message(request, organization.lab_set.all())

    context = {
        "organization": organization,
        "is_member": organization.is_member(request.user) or request.user.is_superuser,
        "organization_form": organization_form,
        'settings_override_form': org_settings_override_form,
        'override_source': override_source,
        'override_values': override_values,
        'has_write_permission': has_write_permission,
    }
    return render(request, 'snpdb/settings/view_organization.html', context)


def custom_columns(request):
    context = {}
    form = forms.CustomColumnsCollectionForm(request.POST or None, user=request.user)

    if request.method == "POST":
        if form.is_valid():
            ccc = form.save()
            return HttpResponseRedirect(reverse("view_custom_columns", kwargs={"custom_columns_collection_id": ccc.pk}))

        add_save_message(request, False, "Columns", created=True)

    context["form"] = form
    return render(request, 'snpdb/settings/custom_columns.html', context)


# Based on code from http://j-syk.com/weblog/2012/10/18/jquery-sortables-ajax-django/

def view_custom_columns(request, custom_columns_collection_id):
    ccc = CustomColumnsCollection.get_for_user(request.user, custom_columns_collection_id)

    custom_columns_qs = VariantGridColumn.objects.filter(customcolumn__custom_columns_collection=ccc)
    my_columns = list(custom_columns_qs.order_by("customcolumn__sort_order"))
    available_columns = list(VariantGridColumn.objects.exclude(grid_column_name__in=my_columns))
    variant_grid_columns = {}
    for vgc in VariantGridColumn.objects.all():
        variant_grid_columns[vgc.pk] = vgc

    has_write_permission = ccc.can_write(request.user)
    if not has_write_permission:
        msg = "You do not have permission to edit these columns. " \
              "If you wish to customise them, click 'clone' and modify the copy"
        messages.add_message(request, messages.WARNING, msg)

    if request.method == "POST":
        ccc.check_can_write(request.user)
        if name := request.POST.get("name"):
            ccc.name = name
            ccc.save()
        elif my_columns_str := request.POST.get("columns"):
            def update_user_columns(id_list, active):
                for i, col in enumerate(id_list):
                    column = variant_grid_columns[col]
                    CustomColumn.objects.update_or_create(custom_columns_collection=ccc, column=column,
                                                          defaults={"sort_order": i})
                # Delete any not in id_list
                CustomColumn.objects.filter(custom_columns_collection=ccc).exclude(column__in=id_list).delete()

            my_columns_list = my_columns_str.split(',') if my_columns_str else []
            active = 'my_columns' in request.POST
            update_user_columns(my_columns_list, active)
        return HttpResponse()  # Nobody ever looks at this

    context_dict = {
        'available_columns_list': available_columns,
        'my_columns_list': my_columns,
        'custom_columns': ccc,
        'has_write_permission': has_write_permission,
    }
    return render(request, 'snpdb/settings/view_custom_columns.html', context_dict)


def tag_settings(request):
    form = forms.CreateTagForm(request.POST or None)
    if request.method == "POST":
        valid = form.is_valid()
        if valid:
            tag_name = form.cleaned_data['tag']
            name = f"Tag {tag_name}"
            try:
                Tag.objects.create(pk=tag_name)
            except:
                valid = False
        else:
            name = "Tag"
        add_save_message(request, valid, name, created=True)

    user_tag_styles, user_tag_colors = UserTagColors.get_tag_styles_and_colors(request.user)
    context_dict = {'form': form,
                    'user_tag_styles': user_tag_styles,
                    'user_tag_colors': user_tag_colors}
    return render(request, 'snpdb/settings/tag_settings.html', context_dict)


@require_POST
def set_user_tag_color(request):
    tag = request.POST['tag']
    rgb = request.POST['rgb']
    (utc, _) = UserTagColors.objects.get_or_create(user=request.user, tag_id=tag)
    utc.rgb = rgb
    utc.save()
    logging.info("saved %s", utc)
    return HttpResponse()


def igv_integration(request):
    widgets = {"prefix": TextInput(attrs={'placeholder': 'from...'}),
               "replacement": TextInput(attrs={'placeholder': 'to...'})}

    UserDataPrefixFormSet = inlineformset_factory(User,
                                                  UserDataPrefix,
                                                  can_delete=True,
                                                  fields=ALL_FIELDS,
                                                  widgets=widgets,
                                                  max_num=10,
                                                  extra=3)
    if request.method == "POST":
        formset = UserDataPrefixFormSet(request.POST, instance=request.user)
        valid = formset.is_valid()
        if valid:
            formset.save()
        add_save_message(request, valid, "IGV Integration")

    formset = UserDataPrefixFormSet(instance=request.user)
    context_dict = {'user': request.user,
                    'formset': formset,
                    'example_replacements': get_example_replacements(request.user)}
    return render(request, 'snpdb/settings/igv_integration.html', context_dict)


def cohorts(request):
    user_settings = UserSettings.get_for_user(request.user)
    initial = {'user': request.user, 'genome_build': user_settings.default_genome_build}
    form = forms.CreateCohortForm(request.POST or None, initial=initial)
    if request.method == "POST":
        valid = form.is_valid()
        if valid:
            cohort = form.save()
            return HttpResponseRedirect(reverse('view_cohort', kwargs={'cohort_id': cohort.pk}))
        else:
            add_save_message(request, valid, "Cohort", created=True)

    context = {"form": form}
    return render(request, 'snpdb/patients/cohorts.html', context)


def view_cohort_details_tab(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    context = {"cohort": cohort,
               "has_write_permission": cohort.can_write(request.user)}
    return render(request, 'snpdb/patients/view_cohort_details_tab.html', context)


def view_cohort(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if cohort.vcf:
        return redirect('view_vcf', vcf_id=cohort.vcf.pk)

    try:
        cohort_genotype_collection = cohort.cohort_genotype_collection
    except CohortGenotypeCollection.DoesNotExist:
        cohort_genotype_collection = None

    form = forms.CohortForm(request.POST or None, instance=cohort)
    if request.method == "POST":
        if valid := form.is_valid():
            cohort = form.save()
        add_save_message(request, valid, "Cohort")

    sample_form = SampleChoiceForm(genome_build=cohort.genome_build)
    sample_form.fields['sample'].required = False

    context = {"form": form,
               "sample_form": sample_form,
               "cohort": cohort,
               "cohort_genotype_collection": cohort_genotype_collection,
               "has_write_permission": cohort.can_write(request.user)}
    return render(request, 'snpdb/patients/view_cohort.html', context)


def cohort_sample_edit(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)

    if request.method == "POST":
        cohort_op = request.POST['cohort_op']
        sample_ids_str = request.POST['sample_ids']
        sample_ids = json.loads(sample_ids_str)
        if cohort_op == 'add':
            for sample_id in sample_ids:
                cohort.add_sample(sample_id)
        elif cohort_op == 'remove':
            for sample_id in sample_ids:
                try:
                    cohort_sample = CohortSample.objects.get(cohort=cohort, sample_id=sample_id)
                    cohort_sample.delete()
                    logging.info("Removed: %s", sample_id)
                except CohortSample.DoesNotExist:
                    pass
        else:
            raise ValueError(f"Unknown cohort_op '{cohort_op}'")

    return HttpResponse()


def cohort_hotspot(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    form = GeneAndTranscriptForm(genome_build=cohort.genome_build)

    try:
        cohort_genotype_collection = cohort.cohort_genotype_collection
    except Exception as e:
        cohort_genotype_collection = None
        logging.error(e)

    context = {"cohort": cohort,
               "cohort_genotype_collection": cohort_genotype_collection,
               "form": form}
    return render(request, 'snpdb/patients/cohort_hotspot.html', context)


def cohort_gene_counts(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)

    COHORT_CUSTOM_GENE_LIST = f"__QC_COVERAGE_CUSTOM_GENE_LIST__{request.user}"

    # We only want to keep 1 per user
    custom_text_gene_list, _ = CustomTextGeneList.objects.get_or_create(name=COHORT_CUSTOM_GENE_LIST)

    custom_gene_list_form = CustomGeneListForm(request.POST or None,
                                               initial={"custom_gene_list_text": custom_text_gene_list.text})
    if custom_gene_list_form.is_valid():
        custom_text_gene_list.text = custom_gene_list_form.cleaned_data['custom_gene_list_text']
        custom_text_gene_list.save()
        create_custom_text_gene_list(custom_text_gene_list, request.user, GeneListCategory.QC_COVERAGE_CUSTOM_TEXT,
                                     hidden=True)
        gene_list_id = custom_text_gene_list.gene_list.pk
    else:
        gene_list_id = None

    context = {"cohort": cohort,
               'gene_list_id': gene_list_id,
               'gene_list_form': UserGeneListForm(),
               'custom_gene_list_form': custom_gene_list_form,
               'gene_count_type_choice_form': GeneCountTypeChoiceForm()}
    return render(request, 'snpdb/patients/cohort_gene_counts.html', context)


def cohort_gene_counts_matrix(request, cohort_id, gene_count_type_id, gene_list_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    gene_count_type = GeneCountType.objects.get(pk=gene_count_type_id)
    gene_list = GeneList.get_for_user(request.user, gene_list_id)
    samples = list(cohort.get_samples())

    annotation_version = AnnotationVersion.latest(cohort.genome_build)
    variant_annotation_version = annotation_version.variant_annotation_version
    cgc, created = CohortGeneCounts.objects.get_or_create(variant_annotation_version=variant_annotation_version,
                                                          gene_count_type=gene_count_type,
                                                          cohort=cohort,
                                                          cohort_version=cohort.version)

    graph_kwargs = {"cohort_id": cohort_id,
                    "gene_count_type_id": gene_count_type_id,
                    "gene_list_id": gene_list_id}
    redirect_url = reverse("cohort_gene_counts_matrix", kwargs=graph_kwargs)
    if created or (cgc.processing_status not in ProcessingStatus.FINISHED_STATES):
        celery_task = cgc.launch_task()
        wait_for_task_kwargs = {"celery_task": celery_task, "sleep_ms": 2000, "redirect_url": redirect_url}
        wait_url = reverse("wait_for_task", kwargs=wait_for_task_kwargs)
        return HttpResponseRedirect(wait_url)
    else:
        if cgc.processing_status == ProcessingStatus.SUCCESS:
            return sample_gene_matrix(request, variant_annotation_version, samples, gene_list, gene_count_type)
        else:
            raise ValueError(f"{cgc} had ProcessingStatus: {cgc.processing_status}")


def trios(request):
    context = {}
    return render(request, 'snpdb/patients/trios.html', context)


def view_trio(request, pk):
    trio = Trio.get_for_user(request.user, pk)
    context = {"trio": trio,
               "has_write_permission": trio.cohort.can_write(request.user)}
    return render(request, 'snpdb/patients/view_trio.html', context)


@login_not_required
def sample_gene_matrix(request, variant_annotation_version, samples, gene_list,
                       gene_count_type, highlight_gene_symbols=None):
    """ highlight_gene_symbols - put these genes 1st """
    # 19/07/18 - Plotly can't display a categorical color map. See: https://github.com/plotly/plotly.js/issues/1747
    # So just doing as HTML table

    if gene_list:
        genes = gene_list.get_genes(variant_annotation_version.gene_annotation_release)
        gene_symbols = set(gene_list.get_gene_names())
    else:
        # This was originally designed around a gene list, but now we need to support no gene list (only when uses
        # variant classifications)
        genes = []
        gene_symbols = []
        qs = gene_count_type.get_variant_queryset(variant_annotation_version)
        GS_PATH = "variantannotation__transcript_version__gene_version__gene_symbol"
        qs = qs.filter(**{GS_PATH + "__isnull": False})
        for gene, gene_symbol in qs.values_list("variantannotation__gene", GS_PATH).distinct():
            genes.append(gene)
            gene_symbols.append(gene_symbol)

    gene_values = list(gene_count_type.genevalue_set.all().order_by("id"))
    default_color = "#d9d9d9"
    default_text = ""
    empty_gene_value = list(filter(lambda x: x.use_as_empty_value, gene_values))
    if len(empty_gene_value) == 1:
        default_color = empty_gene_value[0].rgb

    phenotypes = ["Age", "HPO", "OMIM"]
    highlight_gene_labels = []
    other_gene_labels = []

    gene_links_lookup = OrderedDict()
    for gene_symbol in sorted(gene_symbols):
        gene_classes_list = ["gene-label", gene_symbol]

        highlight = highlight_gene_symbols and gene_symbol in highlight_gene_symbols
        if highlight:
            gene_classes_list.append("highlight-gene")
        gene_classes = ' '.join(gene_classes_list)

        if request.user.is_authenticated:  # Only display links to logged in users
            url = reverse('view_gene_symbol', kwargs={"gene_symbol": gene_symbol})
            gene_symbol_text = f'<a class="{gene_classes}" href="{url}">{gene_symbol}</a>'
        else:
            gene_symbol_text = f"<span class='{gene_classes}'>{gene_symbol}</span>"

        if highlight:
            highlight_gene_labels.append(gene_symbol_text)
        else:
            other_gene_labels.append(gene_symbol_text)
        gene_links_lookup[gene_symbol] = gene_symbol_text

    matrix_rows = phenotypes + highlight_gene_labels + other_gene_labels
    color_df = pd.DataFrame(index=matrix_rows, dtype='O')
    text_df = pd.DataFrame(index=matrix_rows)

    sample_names = []
    used_sample_names = set()

    for i, sample in enumerate(samples):
        try:
            can_access = False
            if request.user.is_authenticated:  # Only display links to logged in users
                try:
                    Sample.get_for_user(request.user, sample.pk)  # Throws exception
                    can_access = True
                except (Sample.DoesNotExist, PermissionDenied):
                    pass

            source = SampleAnnotationVersionVariantSource.objects.get(sample=sample,
                                                                      variant_annotation_version=variant_annotation_version)

            gvcc = GeneValueCountCollection.objects.get(source=source,
                                                        gene_count_type=gene_count_type)
            gvc_qs = gvcc.genevaluecount_set.filter(gene__in=genes)

            sample_code = "%03d" % i
            if can_access:
                view_sample_url = reverse('view_sample', kwargs={'sample_id': sample.pk})

                sample_link = f'<a href="{view_sample_url}">{sample.name}</a>'
                if sample_link in used_sample_names:
                    uniq_sample_name = sample.name + "_" + sample_code
                    sample_link = f'<a href="{view_sample_url}">{uniq_sample_name}</a>'

                sample_name = sample_link
            else:
                sample_name = "S" + sample_code

            sample_names.append(sample_name)
            used_sample_names.add(sample_name)

            color_df[sample_name] = default_color
            color_df.loc["Age", sample_name] = '#FFFFFF'
            color_df.loc["HPO", sample_name] = '#FFFFFF'
            color_df.loc["OMIM", sample_name] = '#FFFFFF'

            text_df[sample_name] = default_text

            if sample.patient:
                try:
                    try:
                        age = sample.specimen.age_at_collection_date
                    except:
                        age = None
                    text_df.loc["Age", sample_name] = age or ''

                    # Check you have Patient permissions
                    patient = Patient.get_for_user(request.user, sample.patient.pk)
                    terms_dict = OntologyTerm.split_hpo_omim_mondo_as_dict(patient.get_ontology_term_ids())

                    def format_ontology(ontology_term):
                        return f"<div title='{ontology_term}'>{ontology_term.name}</div>"

                    for ontology_name, terms_qs in terms_dict.items():
                        ontology_text = " ".join(map(format_ontology, terms_qs))
                        text_df.loc[ontology_name, sample_name] = ontology_text
                except PermissionDenied:
                    pass
                except Patient.DoesNotExist:
                    pass

            FIELDS = ["gene__geneversion__gene_symbol", "value__rgb", "value__show_counts", "count"]
            for gene_symbol, rgb, show_counts, count in gvc_qs.values_list(*FIELDS):
                gene_link = gene_links_lookup[gene_symbol]
                color_df.loc[gene_link, sample_name] = rgb
                if show_counts:
                    text_df.loc[gene_link, sample_name] = count
        except (SampleAnnotationVersionVariantSource.DoesNotExist, GeneValueCountCollection.DoesNotExist):
            pass

    def set_style(s):
        color_series = color_df[s.name]
        styles = []
        for color in color_series:
            styles.append(f"color: {rgb_invert(color)}; background-color: {color};")

        return styles

    style = text_df.style.apply(set_style)
    style = style.set_table_attributes('class="sample-gene-matrix"')
    text_table_html = style.render()

    context = {"text_table_html": text_table_html,
               "gene_values": gene_values}
    return render(request, 'snpdb/patients/cohort_gene_counts_matrix.html', context)


def cohort_sort(request, cohort_id):
    cohort = Cohort.get_for_user(request.user, cohort_id)
    if request.method == "POST":
        cohort_samples_str = request.POST.get("cohort_samples")
        cohort_samples_ids = cohort_samples_str.split(',') if cohort_samples_str else []
        cohort_samples = []
        for i, cs_id in enumerate(cohort_samples_ids):
            cohort_sample = CohortSample.objects.get(pk=cs_id, cohort=cohort)
            cohort_sample.sort_order = i
            cohort_sample.save()
            cohort_samples.append(cohort_sample)
    else:
        cohort_samples = cohort.get_cohort_samples()

    context = {'cohort': cohort,
               'cohort_samples': cohort_samples}
    return render(request, 'snpdb/patients/cohort_sort.html', context)


def help_static_page(request, page_name):
    """ This embeds static pages in a help template """
    context = {"page_name": page_name}
    return render(request, 'snpdb/help/help_static_page.html', context)


def ajax_hello_world(request, data:str):
    return render(request, 'snpdb/ajax_hello_world.html', {'data': data})


def staff_only(request):
    return render(request, 'snpdb/staff_only.html')


@cache_page(WEEK_SECS)
def tag_autocomplete_form(request):
    """ This is an absolutely minimal HTML to create a Tag autocomplete form (used for load()) """

    context = {"tag_form": TagForm()}
    return render(request, 'snpdb/tag_autocomplete_form.html', context)


def wait_for_task(request, celery_task, sleep_ms, redirect_url):
    async_result = AsyncResult(celery_task)
    if async_result.successful():
        return HttpResponseRedirect(redirect_url)

    kwargs = {"celery_task": celery_task, "sleep_ms": sleep_ms, "redirect_url": redirect_url}
    url = reverse("wait_for_task", kwargs=kwargs)

    context = {"url": url,
               "sleep_ms": sleep_ms,
               "async_result": async_result}
    return render(request, 'snpdb/wait_for_task.html', context)


@require_POST
def wiki_save(request, class_name, unique_keyword, unique_value):
    wiki = Wiki.get_or_create(class_name, unique_keyword, unique_value)
    markdown = request.POST["markdown"]
    wiki.check_user_edit_permission(request.user)  # Throws 403

    wiki.markdown = markdown
    wiki.last_edited_by = request.user
    wiki.save()
    return JsonResponse({})


def labs(request):
    # Use short names is available
    show_unclassified = not settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED

    """
    vc_state_data_json = get_grouped_classification_counts(
        user=request.user,
        field=state_field,
        max_groups=15,
        show_unclassified=show_unclassified)
    """
    active_organizations = Organization.objects.filter(active=True).order_by('name')
    organization_labs = {}
    for org in active_organizations:
        if settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED:
            org_labs = org.sharing_labs
        else:
            org_labs = org.classifying_labs
        if org_labs:
            organization_labs[org] = list(org_labs)
    lab_list = [l for ll in organization_labs.values() for l in ll]
    context = {
        "organization_labs": organization_labs,
        "labs": lab_list,
        "shared_classifications": settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED,
        # "vc_state_data": vc_state_data_json,
        "show_unclassified": show_unclassified,
    }

    return render(request, "snpdb/labs.html", context)


def labs_graph_detail(request):
    short_names_qs = Organization.objects.filter(short_name__isnull=False)
    name_to_short_name = dict(short_names_qs.values_list("name", "short_name"))
    state_field = "classification__lab__state"
    show_unclassified = not settings.VARIANT_CLASSIFICATION_STATS_USE_SHARED
    org_field = "classification__lab__organization__name"

    vc_org_data_json = get_grouped_classification_counts(
        user=request.user,
        field=org_field,
        max_groups=15,
        field_labels=name_to_short_name,
        show_unclassified=show_unclassified,
        allele_level=True)

    context = dict()
    context["vc_org_data"] = vc_org_data_json

    graph_data = get_accumulation_graph_data(mode=AccumulationReportMode.Allele)
    context["accumulation_by_status"] = graph_data["status"]
    if request.user.is_superuser:
        # TODO, do we really need to hide this graph away?
        context["accumulation_by_lab"] = graph_data["lab"]

        state_pop_multiplier = {}
        for state in State.objects.filter(population__gt=0):
            state_pop_multiplier[state.name] = 100_000 / state.population

        vc_normalized_state_data_json = get_grouped_classification_counts(
            user=request.user,
            field=state_field,
            max_groups=15,
            show_unclassified=show_unclassified,
            norm_factor=state_pop_multiplier,
            allele_level=True)

        context["vc_normalized_state_data_json"] = vc_normalized_state_data_json
    return render(request, "snpdb/labs_graph_detail.html", context)


@login_not_required
@cache_page(4 * HOUR_SECS)
def public_global_sample_gene_matrix(request):
    # No auth required - rendered w/o links etc
    return global_sample_gene_matrix(request)


@cache_page(HOUR_SECS)
@vary_on_cookie
def user_global_sample_gene_matrix(request):
    # global_sample_gene_matrix is rendered differently for external/logged in users
    # So keep as separate views so we can cache them
    return global_sample_gene_matrix(request)


def global_sample_gene_matrix(request):
    gene_count_type = GeneCountType.objects.get(pk=settings.PUBLIC_SAMPLE_GENE_MATRIX_TYPE)
    gene_list_id = settings.PUBLIC_SAMPLE_GENE_MATRIX_GENE_LIST_ID
    if gene_list_id:
        gene_list = GeneList.objects.get(pk=gene_list_id)
    else:
        gene_list = None
        if gene_count_type.uses_classifications is False:
            raise PermissionDenied("settings.PUBLIC_SAMPLE_GENE_MATRIX_GENE_LIST_ID must be set "
                                   "if GeneCountType.uses_classifications is False")

    if settings.PUBLIC_SAMPLE_GENE_MATRIX_SHOW_PRIVATE_SAMPLES:
        sample_qs = Sample.objects.filter(import_status=ImportStatus.SUCCESS)
    else:
        public = Group.objects.get(name='public')
        read_perm = DjangoPermission.perm(Sample, DjangoPermission.READ)
        sample_qs = get_objects_for_group(public, read_perm, Sample)

    if gene_count_type.uses_classifications:
        vc_qs = gene_count_type.get_classification_qs()
        sample_qs = sample_qs.filter(classification__in=vc_qs)

    genome_build_name = settings.PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD
    if genome_build_name is None:
        try:
            genome_build = GenomeBuild.builds_with_annotation().get()
        except GenomeBuild.MultipleObjectsReturned:
            msg = f"settings.PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD must be set when there are multiple genome builds"
            raise ImproperlyConfigured(msg)
    else:
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

    samples_list = list(sample_qs.filter(vcf__genome_build=genome_build).order_by("name").distinct())
    variant_annotation_version = VariantAnnotationVersion.latest(genome_build)
    highlight_gene_symbols = settings.PUBLIC_SAMPLE_GENE_MATRIX_HIGHLIGHT_GENE_SYMBOLS

    return sample_gene_matrix(request, variant_annotation_version, samples_list, gene_list, gene_count_type,
                              highlight_gene_symbols=highlight_gene_symbols)


def genomic_intervals_graph(request, genomic_intervals_collection_id):
    graph_class_name = full_class_name(ChromosomeIntervalsGraph)
    cached_graph = graphcache.async_graph(graph_class_name, genomic_intervals_collection_id)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


def chrom_density_graph(request, sample_id, cmap):
    graph_class_name = full_class_name(SampleChromosomeDensityGraph)

    cached_graph = graphcache.async_graph(graph_class_name, cmap, sample_id)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


def homozygosity_graph(request, sample_id, cmap):
    graph_class_name = full_class_name(HomozygosityPercentGraph)
    cached_graph = graphcache.async_graph(graph_class_name, cmap, sample_id)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))


def sample_allele_frequency_histogram_graph(request, sample_id, min_read_depth):
    graph_class_name = full_class_name(AlleleFrequencyHistogramGraph)
    cached_graph = graphcache.async_graph(graph_class_name, sample_id, min_read_depth)
    return HttpResponseRedirect(reverse("cached_generated_file_check", kwargs={"cgf_id": cached_graph.id}))
