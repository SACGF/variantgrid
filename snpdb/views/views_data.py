import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd
from django.conf import settings
from django.contrib import messages
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.forms.models import inlineformset_factory
from django.http import HttpRequest
from django.http.response import (
    HttpResponseRedirect,
)
from django.shortcuts import get_object_or_404, render
from django.template.defaultfilters import pluralize
from django.utils.html import format_html

from analysis.analysis_templates import get_sample_analysis
from analysis.forms import AnalysisOutputNodeChoiceForm
from analysis.models import AnalysisTemplate
from analysis.tasks.analysis_grid_export_tasks import get_annotated_download_files_cgf
from annotation.manual_variant_entry import can_create_variants, create_manual_variants
from annotation.models import (
    AnnotationVersion,
    CohortGenotypeClinVarAnnotationStats,
    CohortGenotypeGeneAnnotationStats,
    CohortGenotypeVariantAnnotationStats,
    VCFAnnotationStats,
)
from annotation.models.models import ManualVariantEntryCollection
from annotation.serializers import ManualVariantEntryCollectionSerializer
from annotation.tasks.calculate_sample_stats import enqueue_cohort_stats_recompute
from classification.views.classification_datatables import ClassificationColumns
from library.django_utils import (
    add_save_message,
    get_model_fields,
    set_form_read_only,
)
from library.utils import full_class_name
from patients.forms import PatientForm
from patients.views import get_patient_upload_csv
from snpdb import forms
from snpdb.archive import (
    ArchivePreconditionError,
    DataArchivedError,
    check_vcf_archive_precondition,
    mark_vcf_archive_started,
)
from snpdb.forms import (
    SampleChoiceForm,
    SampleForm,
    VCFChoiceForm,
)
from snpdb.import_status import set_vcf_and_samples_import_status
from snpdb.models import (
    VCF,
    Cohort,
    CohortGenotypeCollection,
    CohortGenotypeStats,
    GenomeBuild,
    GenomicIntervalsCollection,
    Sample,
    SampleLocusCount,
    SomalierRelatePairs,
    VariantZygosityCountCollection,
    VariantZygosityCountForVCF,
    VCFLengthStatsCollection,
    get_igv_data,
)
from snpdb.models.models_enums import (
    ImportStatus,
)
from snpdb.tasks.vcf_archive_tasks import archive_vcf_task
from upload.models import UploadedVCF
from upload.uploaded_file_type import retry_upload_pipeline
from upload.views.views import get_remaining_annotation_runs


def data(request):
    return render(request, 'snpdb/data/data.html')


def _get_vcf_sample_stats(vcf, passing_filter: bool):
    """ Count is het + hom. Reads CohortGenotypeStats per-sample rows
        (sample IS NOT NULL, filter_key NULL) for the VCF's cohort. """
    try:
        cgc = vcf.cohort.cohort_genotype_collection
    except (Cohort.DoesNotExist, CohortGenotypeCollection.DoesNotExist, DataArchivedError):
        return {}, [], ()

    ss_fields = ("sample_id", "sample__name", "variant_count", "ref_count",
                 "het_count", "hom_count", "unk_count")
    ss_values_qs = (CohortGenotypeStats.objects
                    .filter(cohort_genotype_collection=cgc, sample__vcf=vcf,
                            filter_key__isnull=True, passing_filter=passing_filter)
                    .order_by("sample")
                    .values(*ss_fields))

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


def _get_vcf_length_stats(vcf: VCF) -> dict:
    vcf_length_stats = {}
    try:
        for vcf_ls in vcf.vcflengthstatscollection.vcflengthstats_set.all():
            variant_class = vcf_ls.get_variant_class_display()
            bin_edges = np.array(vcf_ls.histogram_bin_edges)
            bin_centers = [edge / 2 for edge in (bin_edges[:-1] + bin_edges[1:])]
            vcf_length_stats[variant_class] = {
                "is_log": vcf_ls.is_log,
                "bin_centers": bin_centers,
                "counts": vcf_ls.histogram_counts,
                # "width": np.diff(bin_edges),
            }
    except VCFLengthStatsCollection.DoesNotExist:
        pass
    return vcf_length_stats


def _get_vcf_infos(vcf) -> list[str]:
    infos = []
    if vcf.import_status == ImportStatus.IMPORTING:
        if num_annotation_runs := get_remaining_annotation_runs(vcf.uploadedvcf, vcf.genome_build):
            infos.append(f"Variants in VCF are being annotated: {num_annotation_runs} annotation runs remaining")
    return infos


def _get_vcf_skipped_annotation_count(vcf) -> int:
    """ Number of variants VEP was unable to annotate for a VCF (latest annotation version) - see issue #1409 """
    if annotation_version := AnnotationVersion.latest_or_none(vcf.genome_build, validate=False):
        if vas := VCFAnnotationStats.objects.filter(
                vcf=vcf, variant_annotation_version=annotation_version.variant_annotation_version).first():
            return vas.vep_skipped_count
    return 0


def _skipped_annotation_message(count: int, tab_anchor: str) -> str:
    """ HTML warning message linking to the skipped-annotation tab - see issue #1409 """
    plural = pluralize(count)
    return format_html('Variant Effect Predictor (VEP) was unable to annotate {} variant{}. '
                       '<a class="activate-tab" href="#{}">View skipped variants</a>',
                       count, plural, tab_anchor)


def view_vcf(request, vcf_id):
    vcf = VCF.get_for_user(request.user, vcf_id)
    # I couldn't get prefetch_related_objects([vcf], "sample_set__samplestats") to work - so storing in a dict

    sample_stats_het_hom_count, sample_names, sample_zygosities = _get_vcf_sample_stats(vcf, passing_filter=False)
    sample_stats_pass_het_hom_count, _, sample_zygosities_pass = _get_vcf_sample_stats(vcf, passing_filter=True)

    VCFSampleFormSet = inlineformset_factory(VCF, Sample, extra=0, can_delete=False,
                                             fields=["vcf_sample_name", "name", "patient", "extraction"],
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

    cohort_id = None
    try:
        # Some legacy data was too hard to fix and relies on being re-imported
        _ = vcf.cohort
        cohort_id = vcf.cohort.pk
        _ = vcf.cohort.cohort_genotype_collection
    except (Cohort.DoesNotExist, CohortGenotypeCollection.DoesNotExist):
        messages.add_message(request, messages.ERROR, "This legacy VCF is missing data and needs to be reloaded.")
    except DataArchivedError:
        # Banner from _data_archived_banner.html shows the message; cohort_id stays set if available.
        pass

    if reload_vcf:
        set_vcf_and_samples_import_status(vcf, ImportStatus.IMPORTING)
        retry_upload_pipeline(vcf.uploadedvcf.file_upload.uploadpipeline)
        vcf_form = forms.VCFForm(post, instance=vcf)  # Reload as import status has changed
        messages.add_message(request, messages.INFO, "Reloading VCF")

    for warning, _ in vcf.get_warnings():
        messages.add_message(request, messages.WARNING, warning, extra_tags='import-message')

    if not reload_vcf:
        for info in _get_vcf_infos(vcf):
            messages.add_message(request, messages.INFO, info, extra_tags='import-message')

    has_write_permission = vcf.can_write(request.user)
    if not has_write_permission:
        set_form_read_only(vcf_form)
        for form in samples_form.forms:
            set_form_read_only(form)
        messages.add_message(request, messages.WARNING, "You can view but not modify this data.")

    variant_zygosity_count_collections = {}
    for vzcc in VariantZygosityCountCollection.objects.all():
        vzc_vcf = VariantZygosityCountForVCF.objects.filter(vcf=vcf, collection=vzcc).first()
        variant_zygosity_count_collections[vzcc] = vzc_vcf

    try:
        can_view_upload_pipeline = vcf.uploadedvcf.file_upload.can_view(request.user)
    except UploadedVCF.DoesNotExist:
        can_view_upload_pipeline = False

    annotated_download_files = {}
    if not settings.VCF_DOWNLOAD_ADMIN_ONLY or request.user.is_superuser:
        if vcf.import_status == ImportStatus.SUCCESS and cohort_id:
            annotated_download_files = get_annotated_download_files_cgf("export_cohort_to_downloadable_file", cohort_id)

    vcf_length_stats = _get_vcf_length_stats(vcf)

    # VEP-skipped variants for the latest annotation version (VG only - see issue #1409)
    skipped_annotation_count = _get_vcf_skipped_annotation_count(vcf)
    if skipped_annotation_count:
        messages.add_message(request, messages.WARNING,
                             _skipped_annotation_message(skipped_annotation_count, "vcf-skipped-annotation"),
                             extra_tags='html import-message')

    from snpdb.archive import vcf_can_be_archived
    can_archive = has_write_permission and vcf_can_be_archived(vcf)
    restore_source_exists = bool(vcf.data_restorable_from) and os.path.exists(vcf.data_restorable_from)
    restore_source_kind = "uploaded"
    if vcf.data_restorable_from and vcf.data_restorable_from.startswith(settings.PARTITION_ARCHIVE_DIR):
        restore_source_kind = "backend"

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
        'can_view_upload_pipeline': can_view_upload_pipeline,
        'annotated_download_files': annotated_download_files,
        "variant_zygosity_count_collections": variant_zygosity_count_collections,
        "vcf_length_stats": vcf_length_stats,
        "can_archive": can_archive,
        "restore_source_exists": restore_source_exists,
        "restore_source_kind": restore_source_kind,
        "skipped_annotation_count": skipped_annotation_count,
    }
    return render(request, 'snpdb/data/view_vcf.html', context)


def archive_vcf_view(request, vcf_id):
    vcf = VCF.get_for_user(request.user, vcf_id)
    if not vcf.can_write(request.user):
        raise PermissionDenied("You do not have write permission for this VCF")
    if request.method != "POST":
        return HttpResponseRedirect(vcf.get_absolute_url())
    reason = request.POST.get("reason", "")
    force = request.POST.get("force") == "1"
    if vcf.data_archived:
        messages.add_message(request, messages.INFO, "VCF data is already archived.")
        return HttpResponseRedirect(vcf.get_absolute_url())
    if vcf.data_archive_in_progress:
        messages.add_message(request, messages.INFO, "VCF data archiving is already in progress.")
        return HttpResponseRedirect(vcf.get_absolute_url())
    try:
        # Validate synchronously so the user gets immediate feedback, then queue the
        # slow work (zygosity walk + dropping partition data) to avoid request timeout.
        check_vcf_archive_precondition(vcf, force=force)
    except ArchivePreconditionError as e:
        messages.add_message(request, messages.ERROR, str(e))
        return HttpResponseRedirect(vcf.get_absolute_url())

    # Mark in-progress (committed now) so a reload won't offer Archive again, then queue.
    mark_vcf_archive_started(vcf)
    archive_vcf_task.delay(vcf.pk, request.user.pk, reason=reason, force=force)
    if force:
        messages.add_message(request, messages.SUCCESS,
                             "VCF data archiving has been queued (non-recoverable — no source file to restore from).")
    else:
        messages.add_message(request, messages.SUCCESS, "VCF data archiving has been queued.")
    return HttpResponseRedirect(vcf.get_absolute_url())


def restore_vcf_view(request, vcf_id):
    from snpdb.archive import restore_vcf
    vcf = VCF.get_for_user(request.user, vcf_id)
    if not vcf.can_write(request.user):
        raise PermissionDenied("You do not have write permission for this VCF")
    if request.method != "POST":
        return HttpResponseRedirect(vcf.get_absolute_url())
    try:
        restore_vcf(vcf, request.user)
        messages.add_message(request, messages.SUCCESS, "VCF restore started — re-importing data.")
    except ValueError as e:
        messages.add_message(request, messages.ERROR, str(e))
    return HttpResponseRedirect(vcf.get_absolute_url())


def get_patient_upload_csv_for_vcf(request, pk):
    vcf = VCF.get_for_user(request.user, pk)
    sample_qs = vcf.sample_set.all()
    filename = f"vcf_{pk}_patient_upload"
    return get_patient_upload_csv(filename, sample_qs)


def _sample_stats(sample) -> tuple[pd.DataFrame, pd.DataFrame]:
    annotation_version = AnnotationVersion.latest(sample.genome_build)
    try:
        cohort = sample.vcf.cohort
        cgc = cohort.cohort_genotype_collection
    except (Cohort.DoesNotExist, CohortGenotypeCollection.DoesNotExist, DataArchivedError):
        return pd.DataFrame(), pd.DataFrame()

    STATS = {
        "Total": (CohortGenotypeStats, set()),
        "dbSNP": (CohortGenotypeVariantAnnotationStats, {"variant_annotation_version"}),
        "OMIM pheno": (CohortGenotypeGeneAnnotationStats, {"gene_annotation_version"}),
        "ClinVar LP/P": (CohortGenotypeClinVarAnnotationStats, {"clinvar_version"}),
    }

    VARIANT_CLASS = ["variant", "snp", "insertions", "deletions"]
    ZYGOSITY = ["ref", "het", "hom", "unk"]

    variant_class_data = {}
    zygosity_data = {}
    missing_stats = False
    for name, (stats_klass, shared_fields) in STATS.items():
        base_kwargs = {
            "cohort_genotype_collection": cgc,
            "sample": sample,
            "filter_key__isnull": True,
        }
        for sf in shared_fields:
            base_kwargs[sf] = getattr(annotation_version, sf)

        objs = {}
        try:
            objs[name] = stats_klass.objects.get(passing_filter=False, **base_kwargs)
        except ObjectDoesNotExist:
            # Stats absent for the latest annotation version (eg it was bumped since import)
            missing_stats = True

        try:
            objs[f"{name} PASS filters"] = stats_klass.objects.get(passing_filter=True, **base_kwargs)
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

    if missing_stats:
        enqueue_cohort_stats_recompute(cohort, annotation_version)

    return sample_stats_variant_class_df, sample_stats_zygosity_df


def view_sample(request, sample_id):
    sample = Sample.get_for_user(request.user, sample_id)
    has_write_permission = sample.can_write(request.user)

    form = forms.SampleForm(request.POST or None, user=request.user, instance=sample)
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

    for info in _get_vcf_infos(sample.vcf):
        messages.add_message(request, messages.INFO, info, extra_tags='import-message')

    sample_locus_count = list(SampleLocusCount.objects.filter(sample=sample).order_by("locus_count"))
    igv_data = get_igv_data(request.user, genome_build=sample.genome_build)
    patient_form = PatientForm(user=request.user)  # blank
    related_samples = None
    if settings.SOMALIER.get("enabled"):
        related_samples = SomalierRelatePairs.get_for_sample(sample).order_by("relate")

    sample_stats_variant_class_df, sample_stats_zygosity_df = _sample_stats(sample)
    sample_genotype_stats = sample.get_genotype_stats()

    # VEP-skipped variants for the latest annotation version (VG only - see issue #1409)
    skipped_annotation_count = _get_vcf_skipped_annotation_count(sample.vcf)
    if skipped_annotation_count:
        messages.add_message(request, messages.WARNING,
                             _skipped_annotation_message(skipped_annotation_count, "sample-skipped-annotation"),
                             extra_tags='html import-message')

    annotated_download_files = {}
    if not settings.VCF_DOWNLOAD_ADMIN_ONLY or request.user.is_superuser:
        if sample.import_status == ImportStatus.SUCCESS:
            annotated_download_files = get_annotated_download_files_cgf("export_sample_to_downloadable_file", sample.pk)

    context = {
        'annotated_download_files': annotated_download_files,
        'sample': sample,
        'samples': [sample],
        'sample_locus_count': sample_locus_count,
        'form': form,
        'patient_form': patient_form,
        'has_write_permission': has_write_permission,
        'igv_data': igv_data,
        "bam_list": sample.get_bam_files(),
        "sample_stats_variant_class_df": sample_stats_variant_class_df,
        "sample_stats_zygosity_df": sample_stats_zygosity_df,
        "sample_genotype_stats": sample_genotype_stats,
        "related_samples": related_samples,
        "skipped_annotation_count": skipped_annotation_count,
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
    try:
        analysis_template = AnalysisTemplate.get_template_from_setting("ANALYSIS_TEMPLATES_AUTO_SAMPLE")
        analysis = get_sample_analysis(sample, analysis_template)
    except ValueError as e:
        error_message = str(e)

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


def genomic_intervals_graphs_tab(request, genomic_intervals_collection_id):
    gic = get_object_or_404(GenomicIntervalsCollection, pk=genomic_intervals_collection_id)
    if not request.user.has_perm('view_genomicintervalscollection', gic):
        raise PermissionDenied()
    context = {'gic': gic}
    return render(request, 'snpdb/data/genomic_intervals_graphs_tab.html', context)


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


def manual_variant_entry(request):
    if can_create_variants(request.user):
        form = forms.ManualVariantEntryForm(request.POST or None, user=request.user)
        if request.method == 'POST':
            valid = form.is_valid()
            if valid:
                variants_text = form.cleaned_data['variants_text']
                genome_build_pk = form.cleaned_data['genome_build']
                genome_build = GenomeBuild.objects.get(pk=genome_build_pk)
                try:
                    create_manual_variants(request.user, genome_build, variants_text)
                except ValueError as ve:
                    messages.add_message(request, messages.ERROR, ve)
                    valid = False
                form = forms.ManualVariantEntryForm(None, user=request.user)  # Reset form

            add_save_message(request, valid, "Manually entered variants")
    else:
        form = None
        messages.add_message(request, messages.INFO, "Manual variant entry has been disabled by an admin.")

    context = {"form": form}
    return render(request, 'snpdb/data/manual_variant_entry.html', context=context)


def manual_variant_entry_collection_detail(request: HttpRequest, pk: int):
    mvec = ManualVariantEntryCollection.get_for_user(request.user, pk)
    return render(request, 'snpdb/data/manual_variant_entry_collection_detail.html', context={'mvec': mvec})


def watch_manual_variant_entry(request, pk):
    mvec = ManualVariantEntryCollection.get_for_user(request.user, pk)
    # TODO: Quick redirect to variant if it's already ready

    mvec_data = ManualVariantEntryCollectionSerializer(mvec).data
    context = {"mvec": mvec,
               "initial_json": json.dumps(mvec_data)}
    return render(request, 'snpdb/data/watch_manual_variant_entry.html', context=context)
