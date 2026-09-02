from collections import OrderedDict

import pandas as pd
from django.conf import settings
from django.contrib.auth.models import Group
from django.core.exceptions import ImproperlyConfigured, PermissionDenied
from django.shortcuts import render
from django.urls.base import reverse
from django.utils.html import escape
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie
from global_login_required import login_not_required
from guardian.shortcuts import get_objects_for_group

from annotation.models.models import VariantAnnotationVersion
from annotation.models.models_gene_counts import (
    GeneCountType,
    GeneValueCountCollection,
    SampleAnnotationVersionVariantSource,
)
from genes.models import GeneList
from library.constants import HOUR_SECS
from library.guardian_utils import DjangoPermission
from library.utils import rgb_invert
from ontology.models import OntologyTerm
from patients.models import Patient
from snpdb.models import (
    GenomeBuild,
    Sample,
)
from snpdb.models.models_enums import (
    ImportStatus,
)


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
                safe_name = escape(sample.name)
                sample_link = f'<a href="{view_sample_url}">{safe_name}</a>'
                if sample_link in used_sample_names:
                    safe_uniq_name = escape(sample.name + "_" + sample_code)
                    sample_link = f'<a href="{view_sample_url}">{safe_uniq_name}</a>'

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
    text_table_html = style.to_html()

    context = {"text_table_html": text_table_html,
               "gene_values": gene_values}
    return render(request, 'snpdb/patients/cohort_gene_counts_matrix.html', context)


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
            msg = "settings.PUBLIC_SAMPLE_GENE_MATRIX_GENOME_BUILD must be set when there are multiple genome builds"
            raise ImproperlyConfigured(msg)
    else:
        genome_build = GenomeBuild.get_name_or_alias(genome_build_name)

    samples_list = list(sample_qs.filter(vcf__genome_build=genome_build).order_by("name").distinct())
    variant_annotation_version = VariantAnnotationVersion.latest(genome_build)
    highlight_gene_symbols = settings.PUBLIC_SAMPLE_GENE_MATRIX_HIGHLIGHT_GENE_SYMBOLS

    return sample_gene_matrix(request, variant_annotation_version, samples_list, gene_list, gene_count_type,
                              highlight_gene_symbols=highlight_gene_symbols)
