import logging
from collections import defaultdict, Counter

import cdot
from django.conf import settings
from django.contrib import messages
from django.http.response import HttpResponse, HttpResponseRedirect, Http404, \
    JsonResponse
from django.shortcuts import get_object_or_404, render, redirect
from django.template.loader import render_to_string
from django.urls.base import reverse
from django.views.decorators.cache import cache_page
from django.views.decorators.http import require_POST
from django.views.decorators.vary import vary_on_cookie
from htmlmin.decorators import not_minified_response

from annotation.annotation_versions import vav_diff_vs_kwargs
from annotation.clinvar_fetch_request import ClinVarFetchRequest
from annotation.manual_variant_entry import create_manual_variants
from annotation.models import AnnotationVersion, AnnotationRun, VariantAnnotationVersion, \
    VariantAnnotationVersionDiff, Citation
from annotation.models.models import CachedWebResource, HumanProteinAtlasAnnotationVersion, \
    HumanProteinAtlasAnnotation, ColumnVEPField, DBNSFPGeneAnnotationVersion
from annotation.models.models_citations import CitationFetchRequest
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.models.models_version_diff import VersionDiff
from annotation.tasks.annotate_variants import annotation_run_retry
from annotation.tasks.annotation_scheduler_task import annotation_scheduler
from annotation.vep_annotation import get_vep_command, get_vep_variant_annotation_version_kwargs
from genes.models import GeneListCategory, GeneAnnotationImport, GeneVersion, TranscriptVersion, GeneSymbolAlias
from genes.models_enums import AnnotationConsortium, GeneSymbolAliasSource
from library.constants import WEEK_SECS
from library.django_utils import require_superuser, get_field_counts
from library.utils import first
from ontology.models import OntologyTerm, OntologyService, OntologyImport, OntologyVersion
from snpdb.models import VariantGridColumn, SomalierConfig, GenomeBuild, VCF, UserSettings, ColumnAnnotationLevel
from variantgrid.celery import app
from variantgrid.deployment_validation.somalier_check import verify_somalier_config


def get_build_contigs():
    # TODO: We know what build we have when we call this - so could only retrieve for that
    build_contigs = {}
    for genome_build in GenomeBuild.objects.all():
        contigs_list = list(genome_build.contigs.values_list("name", flat=True))
        # Only report build if contigs inserted
        if contigs_list:
            contigs = ', '.join(contigs_list)
            data = {"num_contigs": len(contigs_list),
                    "contigs": contigs}
            build_contigs[genome_build.name] = data
    return build_contigs


def _get_gene_and_transcript_stats(genome_build: GenomeBuild, annotation_consortium):
    genes_and_transcripts = {}
    gene_counts = GeneVersion.objects.filter(gene__annotation_consortium=annotation_consortium,
                                             genome_build=genome_build).count()
    transcripts_qs = TranscriptVersion.objects.filter(transcript__annotation_consortium=annotation_consortium,
                                                      genome_build=genome_build)
    transcript_counts = transcripts_qs.count()
    if gene_counts and transcript_counts:
        genes_and_transcripts = {"genes": gene_counts, "transcripts": transcript_counts}
        field_counts = get_field_counts(transcripts_qs, "import_source__url")
        import_sources = {}
        for import_source in GeneAnnotationImport.objects.filter(url__in=field_counts).order_by("created"):
            field = import_source.url
            num_transcripts = field_counts[field]
            import_sources[field] = {"transcripts": num_transcripts, "created": import_source.created}
        genes_and_transcripts["import_sources"] = import_sources
        if last_tv := transcripts_qs.order_by("pk").last():
            genes_and_transcripts["latest_cdot"] = last_tv.data.get("cdot")

    return genes_and_transcripts


@not_minified_response
def annotation_build_detail(request, genome_build_name):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    build_contigs = get_build_contigs()

    annotation_details = {
        "contigs": build_contigs.get(genome_build.name),
        "annotation_consortium": genome_build.settings["annotation_consortium"],
    }

    expected_somalier_sites_name = None
    reference_ok = False
    try:
        annotation_details["reference_fasta"] = genome_build.reference_fasta
        reference_ok = True
    except Exception as e:
        annotation_details["reference_fasta_error"] = str(e)

    if av := AnnotationVersion.latest(genome_build, active=False, validate=False):
        annotation_details["latest"] = av

        sync_with_current_vep = False
        out_of_sync_with_current_vep = None
        if av.variant_annotation_version is None:
            out_of_sync_with_current_vep = "Not created"
        else:
            vav = av.variant_annotation_version
            vep_vav_kwargs = get_vep_variant_annotation_version_kwargs(vav.genome_build)
            if diff := vav_diff_vs_kwargs(vav, vep_vav_kwargs):
                out_of_sync_with_current_vep = f"Latest is out of sync with current VEP: {diff}"
            else:
                sync_with_current_vep = "In sync with current VEP"

        annotation_details["sync_with_current_vep"] = sync_with_current_vep
        annotation_details["out_of_sync_with_current_vep"] = out_of_sync_with_current_vep
        genes_and_transcripts = None
        try:
            genes_and_transcripts = _get_gene_and_transcript_stats(genome_build, genome_build.annotation_consortium)
            annotation_details["genes_and_transcripts"] = genes_and_transcripts

            annotation_consortia = dict(AnnotationConsortium.choices)
            other_consortia = set(annotation_consortia.keys()) - {genome_build.annotation_consortium}
            other_gene_annotation = {}
            for other_ac in other_consortia:
                annotation_consortium_display = annotation_consortia[other_ac]
                other_gene_annotation[annotation_consortium_display] = _get_gene_and_transcript_stats(genome_build,
                                                                                                      other_ac)
            annotation_details["other_consortia"] = other_gene_annotation
        except:
            pass

        gene_annotation_release = None
        if av.variant_annotation_version:
            if gene_annotation_release := av.variant_annotation_version.gene_annotation_release:
                annotation_details["gene_annotation_release"] = {
                    "name": str(gene_annotation_release),
                    "gene_annotation_import": str(gene_annotation_release.gene_annotation_import),
                }

        if gene_annotation_counts := av.get_gene_annotation().count():
            annotation_details["gene_level_annotation"] = f"{gene_annotation_counts} gene annotations."

        clinvar_counts = av.get_clinvar().count()
        if clinvar_counts:
            annotation_details["clinvar"] = f"{clinvar_counts} ClinVar records"

        clinvar_counts = av.get_clinvar().count()
        annotation_sub_components = [sync_with_current_vep, reference_ok,
                                     genes_and_transcripts, gene_annotation_release, clinvar_counts]
        if settings.ANNOTATION_GENE_ANNOTATION_VERSION_ENABLED:
            annotation_sub_components.append(gene_annotation_counts)

        if settings.SOMALIER.get("enabled"):
            somalier_cfg = SomalierConfig()
            try:
                expected_somalier_sites_name = somalier_cfg.get_sites_vcf_name(genome_build)
                vcf = somalier_cfg.get_sites_vcf(genome_build)
                annotation_details["somalier"] = f"Sites VCF: {vcf.name}"
                somalier = True
            except VCF.DoesNotExist:
                somalier = None
            annotation_sub_components.append(somalier)

        annotation_details["ok"] = all(annotation_sub_components)

    context = {
        "build_name": genome_build.name,
        "details": annotation_details,
        "cdot_version": cdot.__version__,
        "expected_somalier_sites_name": expected_somalier_sites_name,
    }
    return render(request, "annotation/annotation_build_detail.html", context)


def annotation(request):
    context = {
        "genome_build_list": GenomeBuild.builds_with_annotation(),
    }
    return render(request, "annotation/annotation.html", context)


@not_minified_response
def annotation_detail(request):
    # Set Variables to None for uninstalled components, the template will show installation instructions
    ensembl_biomart_transcript_genes = None
    diagnostic_gene_list = None

    gene_symbol_alias_counts = get_field_counts(GeneSymbolAlias.objects.all(), "source")
    if gene_symbol_alias_counts:
        gene_symbol_alias_counts = {GeneSymbolAliasSource(k).label: v for k, v in gene_symbol_alias_counts.items()}

    all_ontologies_accounted_for = True
    ontology_counts = []
    for service in [OntologyService.MONDO, OntologyService.OMIM, OntologyService.HPO, OntologyService.HGNC]:
        # don't report HGNC as it's just there as a stub for other items to relate to
        count = OntologyTerm.objects.filter(ontology_service=service).count()
        ontology_counts.append({"service": service, "count": count})

    ontology_relationship_counts = {}
    ontology_services = [OntologyService.MONDO, OntologyService.OMIM, OntologyService.HPO, OntologyService.HGNC]
    if ontology_version := OntologyVersion.latest(validate=False):
        otr_qs = ontology_version.get_ontology_term_relations()
        for first_index, first_service in enumerate(ontology_services):
            for second_service in ontology_services[first_index:]:
                join_count = otr_qs.filter(source_term__ontology_service=first_service,
                                           dest_term__ontology_service=second_service).count()
                if first_service != second_service:
                    reverse_count = otr_qs.filter(source_term__ontology_service=second_service,
                                                  dest_term__ontology_service=first_service).count()
                    join_count += reverse_count
                ontology_relationship_counts[f"{first_service}{second_service}"] = join_count
                ontology_relationship_counts[f"{second_service}{first_service}"] = join_count

    ontology_imports = []
    for context in ["mondo_file", "gencc_file", "hpo_file", "omim_file", "biomart_omim_aliases", "phenotype_to_genes"]:
        last_import = OntologyImport.objects.filter(context=context).order_by('-created').first()
        if not last_import and context != "omim_file":  # don't complain about omim_file not being imported as not available to environments without license
            all_ontologies_accounted_for = False
        ontology_imports.append({"context": context, "last_import": last_import})

    diagnostic = GeneListCategory.get_or_create_category(category_name='Diagnostic')
    diagnostic_gene_list_count = diagnostic.genelist_set.count()
    if diagnostic_gene_list_count:
        diagnostic_gene_list = f"{diagnostic_gene_list_count} diagnostic gene lists"

    hpa_version = HumanProteinAtlasAnnotationVersion.objects.order_by("-annotation_date").first()
    hpa_counts = HumanProteinAtlasAnnotation.objects.filter(version=hpa_version).count()

    if dbnsfp_gene_annotation := DBNSFPGeneAnnotationVersion.latest():
        num_records = dbnsfp_gene_annotation.dbnsfpgeneannotation_set.count()
        dbnsfp_gene_annotation = f"dbNSFP gene - version: {dbnsfp_gene_annotation.version}, {num_records} records"

    somalier = None
    if somalier_enabled := settings.SOMALIER.get("enabled"):
        somalier = verify_somalier_config()

    # These are empty/None if not set.
    annotations_ok = [all_ontologies_accounted_for,
                      hpa_counts > 0]
    if somalier_enabled:
        annotations_ok.append(somalier)
    annotations_all_imported = all(annotations_ok)  # Any unset will show instructions header

    cached_web_resources = []
    for cwr_name in settings.ANNOTATION_CACHED_WEB_RESOURCES:
        cwr, _ = CachedWebResource.objects.get_or_create(name=cwr_name)
        cached_web_resources.append(cwr)

    context = {
        "annotations_all_imported": annotations_all_imported,
        "all_ontologies_accounted_for": all_ontologies_accounted_for,
        "ensembl_biomart_transcript_genes": ensembl_biomart_transcript_genes,
        "ontology_services": ontology_services,
        "ontology_counts": ontology_counts,
        "ontology_relationship_counts": ontology_relationship_counts,
        "ontology_imports": ontology_imports,
        "gene_symbol_alias_counts": gene_symbol_alias_counts,
        "diagnostic_gene_list": diagnostic_gene_list,
        "hpa_version": hpa_version,
        "hpa_counts": hpa_counts,
        "dbnsfp_gene_annotation": dbnsfp_gene_annotation,
        "num_annotation_columns": VariantGridColumn.objects.count(),
        "cached_web_resources": cached_web_resources,
        "somalier": somalier,
    }
    return render(request, "annotation/annotation_detail.html", context)


@require_POST
@require_superuser
def load_cached_web_resource(request, pk):
    """ Deletes it and creates a new one """
    cwr = get_object_or_404(CachedWebResource, pk=pk)
    cwr.delete()

    # this will call post_save signal handlers - which populate things appropriately
    CachedWebResource.objects.create(name=pk)
    return HttpResponse()


def annotation_versions(request):
    anno_versions = {}
    for genome_build in GenomeBuild.builds_with_annotation():
        try:
            latest = AnnotationVersion.latest(genome_build, validate=False)
        except AnnotationVersion.DoesNotExist:
            latest = None
        qs = AnnotationVersion.objects.filter(genome_build=genome_build).order_by("-annotation_date")
        has_annotation = qs.exists()
        has_active = qs.filter(variant_annotation_version__active=True).exists()
        vep_commands = {}
        for pt in VariantAnnotationPipelineType:
            vep_command = get_vep_command("in.vcf", "out.vcf", genome_build,
                                          genome_build.annotation_consortium, pt)
            vep_command = " ".join(vep_command).replace(" -", "\n")
            vep_commands[pt.value] = vep_command

        anno_versions[genome_build.name] = (vep_commands, qs, latest, has_annotation, has_active)

    context = {
        "annotation_versions": anno_versions,
        "pipeline_types": VariantAnnotationPipelineType,
    }
    return render(request, "annotation/annotation_versions.html", context)


def version_diffs(request):
    def get_last_2(version_klass, version_diff_klass):
        last_2 = list(version_klass.objects.order_by("-pk")[:2])
        data = {"num": len(last_2)}

        try:
            version_from = last_2[1]
            version_to = last_2[0]
            version_diff = version_diff_klass.objects.get(version_from=version_from, version_to=version_to)
            data["diff"] = version_diff
        except:
            pass

        for k, v in zip(['version', 'previous'], last_2):
            data[k] = v

        return data

    VERSION_LEVELS = [('Variant', VariantAnnotationVersion, VariantAnnotationVersionDiff)]

    version_levels = []
    for name, version_klass, version_diff_klass in VERSION_LEVELS:
        version_levels.append((name, get_last_2(version_klass, version_diff_klass)))
    context = {"version_levels": version_levels}
    return render(request, "annotation/version_diffs.html", context)


def view_version_diff(request, version_diff_id):
    try:
        diff = VersionDiff.objects.get_subclass(pk=version_diff_id)
    except VersionDiff.DoesNotExist:
        raise Http404(f"No VersionDiff pk={version_diff_id}")

    context = {"diff": diff}
    context.update(diff.get_diff_results())
    return render(request, "annotation/view_version_diff.html", context)


@require_superuser
def variant_annotation_runs(request):
    as_display = dict(AnnotationStatus.choices)

    genome_build_field_counts = defaultdict(dict)
    genome_build_summary = defaultdict(dict)

    if request.method == "POST":
        if "annotation-scheduler" in request.POST:
            annotation_scheduler.si().apply_async()

        for genome_build in GenomeBuild.builds_with_annotation():
            for vav in VariantAnnotationVersion.objects.filter(genome_build=genome_build):
                annotation_runs = AnnotationRun.objects.filter(annotation_range_lock__version=vav)
                message = None
                if f"set-non-finished-to-error-{genome_build.name}-{vav.pk}" in request.POST:
                    num_errored = 0
                    non_finished_statuses = [AnnotationStatus.FINISHED, AnnotationStatus.ERROR]
                    for annotation_run in annotation_runs.exclude(status__in=non_finished_statuses):
                        if celery_task := annotation_run.task_id:
                            logging.info("Terminating celery job '%s'", celery_task)
                            app.control.revoke(celery_task, terminate=True)  # @UndefinedVariable
                        annotation_run.error_exception = "Manually failed"
                        annotation_run.save()
                        num_errored += 1
                    message = f"{genome_build} - set {num_errored} annotation runs to Error"
                elif f"retry-annotation-runs-{genome_build.name}-{vav.pk}" in request.POST:
                    num_retrying = 0
                    for annotation_run in annotation_runs.filter(status=AnnotationStatus.ERROR):
                        annotation_run_retry(annotation_run)
                        num_retrying += 1
                    message = f"{genome_build} - retrying {num_retrying} annotation runs."

                if message:
                    messages.add_message(request, messages.INFO, message)

    for genome_build in GenomeBuild.builds_with_annotation():
        for vav in VariantAnnotationVersion.objects.filter(genome_build=genome_build).order_by("-annotation_date"):
            qs = AnnotationRun.objects.filter(annotation_range_lock__version=vav)
            field_counts = get_field_counts(qs, "status")
            summary_data = Counter()
            for field, count in field_counts.items():
                summary = AnnotationStatus.get_summary_state(field)
                summary_data[summary] += count

            genome_build_summary[genome_build.pk][vav.pk] = summary_data
            genome_build_field_counts[genome_build.pk][vav] = {as_display[k]: v for k, v in field_counts.items()}

    current_variant_annotation_versions = VariantAnnotationVersion.latest_for_all_builds()
    current_pks = current_variant_annotation_versions.values_list("pk", flat=True)
    historical_variant_annotation_versions = VariantAnnotationVersion.objects.exclude(pk__in=current_pks)

    context = {
        "genome_build_summary": dict(genome_build_summary),
        "genome_build_field_counts": dict(genome_build_field_counts),
        "current_variant_annotation_versions": current_variant_annotation_versions,
        "historical_variant_annotation_versions": historical_variant_annotation_versions,
    }
    return render(request, "annotation/variant_annotation_runs.html", context)


def view_annotation_run(request, annotation_run_id):
    annotation_run = get_object_or_404(AnnotationRun, pk=annotation_run_id)

    # There may be other runs of different types (don't care about them)
    other_annotation_runs_qs = AnnotationRun.objects.filter(
        annotation_range_lock=annotation_run.annotation_range_lock,
        pipeline_type=annotation_run.pipeline_type).exclude(pk=annotation_run.pk)
    can_retry_annotation_run = not other_annotation_runs_qs.exists()
    can_retry_annotation_run_upload = can_retry_annotation_run and annotation_run.vcf_annotated_filename

    context = {"annotation_run": annotation_run,
               "can_retry_annotation_run": can_retry_annotation_run,
               "can_retry_annotation_run_upload": can_retry_annotation_run_upload}
    return render(request, "annotation/view_annotation_run.html", context)


@require_POST
def retry_annotation_run(request, annotation_run_id, upload_only=False):
    """ Deletes - then re-tries using annotation lock """
    annotation_run = get_object_or_404(AnnotationRun, pk=annotation_run_id)
    annotation_run = annotation_run_retry(annotation_run, upload_only=upload_only)

    msg = 'Re-trying annotation'
    if upload_only:
        msg += " (upload only)"
    status = messages.INFO
    messages.add_message(request, status, msg, extra_tags='import-message')

    return HttpResponseRedirect(reverse("view_annotation_run", kwargs={"annotation_run_id": annotation_run.pk}))


@require_POST
def retry_annotation_run_upload(request, annotation_run_id):
    return retry_annotation_run(request, annotation_run_id, upload_only=True)


@cache_page(WEEK_SECS)
@vary_on_cookie  # the information isn't actually different per user, but hack to avoid showing other user's email/notifications etc in the top right
def view_annotation_descriptions(request, genome_build_name=None):
    genome_build = UserSettings.get_genome_build_or_default(request.user, genome_build_name)
    variantgrid_columns_by_annotation_level = defaultdict(list)
    vep_annotation_levels = [ColumnAnnotationLevel.TRANSCRIPT_LEVEL, ColumnAnnotationLevel.VARIANT_LEVEL]
    columns_and_vep_by_annotation_level = {al.label: {} for al in vep_annotation_levels}

    vep_qs = ColumnVEPField.filter_for_build(genome_build)
    for vgc in VariantGridColumn.objects.all().order_by("grid_column_name"):
        if vgc.annotation_level in vep_annotation_levels:
            # For Transcript/Variant that use VEP - only show if visible in that build
            if vep := vep_qs.filter(variant_grid_column=vgc).first():
                columns_and_vep_by_annotation_level[vgc.get_annotation_level_display()][vgc] = vep
        else:
            variantgrid_columns_by_annotation_level[vgc.annotation_level].append(vgc)

    context = {
        "genome_build": genome_build,
        "variantgrid_columns_by_annotation_level": variantgrid_columns_by_annotation_level,
        "columns_and_vep_by_annotation_level": columns_and_vep_by_annotation_level
    }
    return render(request, "annotation/view_annotation_descriptions.html", context)


@cache_page(WEEK_SECS)
@vary_on_cookie  # the information isn't actually different per user, but hack to avoid showing other user's email/notifications etc in the top right
def about_new_vep_columns(request):
    # Keep this around for another upgrade cycle
    return render(request, "annotation/about_new_vep_columns.html")


def view_annotation_version_details(request, annotation_version_id):
    annotation_version = get_object_or_404(AnnotationVersion, pk=annotation_version_id)

    context = {"annotation_version": annotation_version}
    return render(request, "annotation/view_annotation_version_details.html", context)


def create_manual_variant_entry_from_text(request, genome_build_name, variants_text):
    genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
    mvec = create_manual_variants(request.user, genome_build, variants_text)
    return redirect('watch_manual_variant_entry', pk=mvec.pk)


def view_citation(request, citation_id: str):
    """
    A dedicated page for a single citation, handy for testing.
    :param request: The Request
    :param citation_id: The full ID of the citation, e.g. PMID:4353345
    :return: A stand alone page with the loaded citation details
    """
    # important, we don't actually retrieve the citation here
    return render(request, "annotation/citation.html", {"citation": first(CitationFetchRequest.get_unfetched_citations([citation_id]))})


def view_citation_detail(request, citation_id: str):
    """
    A dedicated page for a single citation, handy for testing
    :param request: The Request
    :param citation_id: The full ID of the citation, e.g. PMID:4353345
    :return: A stand alone page with the loaded citation details
    """
    return render(request, "annotation/citation_detail.html", {"citation": CitationFetchRequest.fetch_all_now([citation_id]).first_citation})


def simple_citation_html(cd: Citation) -> str:
    return render_to_string('annotation/citation_simple.html', {"citation": cd}).replace('\n', '').strip()


def citations_json(request, citations_ids_list):
    """
    Request JSON of citations, accepts either variant grid internal citation ids, or PubMed:123456
    """
    citation_ids = citations_ids_list.split("/")
    return JsonResponse({"citations": CitationFetchRequest.fetch_all_now(citation_ids).to_json()})


def view_clinvar_records_detail(request, clinvar_variation_id: int, min_stars: int):
    clinvar_record_collection = ClinVarFetchRequest(
        clinvar_variation_id=clinvar_variation_id,
    ).fetch()
    return render(request, "annotation/clinvar_records_detail.html", {
        "clinvar_record_collection": clinvar_record_collection,
        "records": clinvar_record_collection.records_with_min_stars(min_stars)
    })
