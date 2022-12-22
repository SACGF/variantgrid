import csv
from collections import defaultdict, OrderedDict

from django.core.exceptions import PermissionDenied
from django.http.response import StreamingHttpResponse
from django.shortcuts import get_object_or_404, render, redirect
from django.urls.base import reverse
from django.views.decorators.cache import cache_page
from django.views.decorators.vary import vary_on_cookie

from analysis.forms import KaryomappingGeneForm, UserTrioForm
from analysis.models.models_karyomapping import KaryomappingAnalysis, KaryotypeBins, \
    KaryomappingGene
from library.constants import DAY_SECS
from library.django_utils import add_save_message
from library.jqgrid.jqgrid_export import StashFile
from patients.models_enums import Zygosity
from snpdb.models import Trio
from snpdb.models.models_variant import Variant


def get_karyomapping_analysis_permission_check(request, pk):
    ka = get_object_or_404(KaryomappingAnalysis, pk=pk)
    if not request.user.has_perm(KaryomappingAnalysis.get_read_perm(), ka):
        msg = f"{request.user} does not have permission to access {ka}"
        raise PermissionDenied(msg)
    return ka


def karyomapping_analyses(request):
    context = {"trio_form": UserTrioForm()}
    return render(request, 'analysis/karyomapping/karyomapping_analyses.html', context)


def create_and_view_karyomapping_analysis_for_trio(trio, user):
    karyomapping = KaryomappingAnalysis.objects.create(user=user,
                                                       name=trio.name,
                                                       trio=trio)
    url = reverse("view_karyomapping_analysis", kwargs={"pk": karyomapping.pk})
    return redirect(url)


def create_karyomapping_analysis_for_trio_id(request, trio_id):
    trio = Trio.get_for_user(request.user, trio_id)
    return create_and_view_karyomapping_analysis_for_trio(trio, request.user)


def view_karyomapping_analysis(request, pk):
    karyomapping_analysis = get_karyomapping_analysis_permission_check(request, pk)

    gene_form = KaryomappingGeneForm(request.POST or None,
                                     karyomapping_analysis=karyomapping_analysis,
                                     initial={"upstream_kb": 2000,
                                              "downstream_kb": 2000})
    created_karyomapping_gene = None
    if request.method == "POST":
        valid = gene_form.is_valid()
        if valid:
            created_karyomapping_gene = gene_form.save()

        add_save_message(request, valid, "KaryomappingGene")

    context = {"karyomapping_analysis": karyomapping_analysis,
               "gene_form": gene_form,
               "has_write_permission": karyomapping_analysis.can_write(request.user),
               "created_karyomapping_gene": created_karyomapping_gene}
    return render(request, 'analysis/karyomapping/view_karyomapping_analysis.html', context)


def get_variant_lookup_and_scatter_data(karyomapping_bins):
    """ Dumped to JS to be used by Plotly scatterplot
        karyomapping_bins : Have separate entries for ALT/REF, we merge these for output """

    variant_id_lookup = {}
    data = defaultdict(lambda: defaultdict(list))
    for karyotype_code, variant_data in karyomapping_bins.items():
        x = []
        text = []
        for variant_id, chrom, position, ref, alt in variant_data:

            variant_string = Variant.format_tuple(chrom, position, ref, alt)
            variant_id_lookup[variant_string] = variant_id
            x.append(position)
            text.append(variant_string)

        collapsed_code = KaryotypeBins.COLLAPSED_BINS[karyotype_code]
        data[collapsed_code]["x"].extend(x)
        data[collapsed_code]["text"].extend(text)

    karyotype_bin_counts = OrderedDict()
    for kc in KaryotypeBins.KARYOTYPE_LABEL_ORDER:
        karyotype_bin_counts[kc] = len(data[kc]["x"])

    return variant_id_lookup, data, karyotype_bin_counts


@cache_page(DAY_SECS)  # Only caching this for a day due to high amount of development
@vary_on_cookie
def view_karyomapping_gene(request, pk):
    karyomapping_gene = get_object_or_404(KaryomappingGene, pk=pk)
    # Permission check on parent karyomapping_analysis
    get_karyomapping_analysis_permission_check(request, karyomapping_gene.karyomapping_analysis.pk)

    iv, strand = karyomapping_gene.get_genomic_interval_and_strand()
    karyomapping_bins = karyomapping_gene.get_karyomapping_bins()
    variant_id_lookup, karyotype_bin_scatter_data, karyotype_bin_counts = get_variant_lookup_and_scatter_data(karyomapping_bins)

    context = {"kag": karyomapping_gene,
               "iv": iv,
               "strand": strand,
               "karyotype_bin_labels": KaryotypeBins.KARYOTYPE_LABEL_ORDER,
               "karyotype_bin_labels": KaryotypeBins.KARYOTYPE_LABEL_ORDER,  # Want order
               "variant_id_lookup": variant_id_lookup,
               "karyotype_bin_scatter_data": karyotype_bin_scatter_data,
               "karyotype_bin_counts": karyotype_bin_counts}
    return render(request, 'analysis/karyomapping/view_karyomapping_gene.html', context)


def download_karyomapping_gene_csv(request, pk):
    karyomapping_gene = get_object_or_404(KaryomappingGene, pk=pk)
    # Permission check on parent karyomapping_analysis
    get_karyomapping_analysis_permission_check(request, karyomapping_gene.karyomapping_analysis.pk)

    variant_and_genotypes = karyomapping_gene.get_variant_and_genotypes()

    filename = f"karyomapping_gene_{karyomapping_gene.pk}_{karyomapping_gene}.csv"
    # TODO: merge code w/library.jqgrid_export.grid_export_csv

    karotype_bin_lookup = KaryotypeBins.get_karotype_bin_lookup()

    header = ['chrom', 'position', 'ref', 'alt', 'proband_gt', 'father_gt', 'mother_gt', 'karyotype_bin']
    pseudo_buffer = StashFile()
    writer = csv.DictWriter(pseudo_buffer, header, dialect='excel')

    def iter_row_writer():
        writer.writeheader()
        yield pseudo_buffer.value
        for variant_data, genotype_tuple in variant_and_genotypes:
            _, chrom, position, ref, alt = variant_data
            proband_gt, father_gt, mother_gt = genotype_tuple
            try:
                karotype_bin = karotype_bin_lookup[proband_gt][father_gt][mother_gt]
            except:
                karotype_bin = ''

            row = {'chrom': chrom,
                   'position': position,
                   'ref': ref,
                   'alt': alt,
                   'proband_gt': Zygosity.get_genotype(proband_gt),
                   'father_gt': Zygosity.get_genotype(father_gt),
                   'mother_gt': Zygosity.get_genotype(mother_gt),
                   'karyotype_bin': karotype_bin}
            writer.writerow(row)
            yield pseudo_buffer.value

    response = StreamingHttpResponse(iter_row_writer(), content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response
