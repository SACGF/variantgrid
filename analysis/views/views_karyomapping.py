"""
Per-gene karyomapping for a trio (#1754): pick a gene and flanks, see the proband's phasing variants by bin as a
Plotly scatter and download them as CSV. The genome-wide QC counts are on the trio page (karyomapping_tags).
"""
import csv
from collections import defaultdict
from typing import Optional

from django.core.exceptions import PermissionDenied
from django.http import Http404
from django.http.response import StreamingHttpResponse
from django.shortcuts import render

from analysis.forms import TrioKaryomappingGeneForm
from analysis.models.models_karyomapping import KaryotypeBins, TrioKaryomappingGeneRegion
from library.utils import StashFile
from patients.models_enums import Zygosity
from snpdb.models import Trio
from snpdb.models.models_variant import Variant


def _get_trio(request, trio_id) -> Trio:
    trio = Trio.get_for_user(request.user, trio_id)
    if trio.cohort.data_archived:
        raise PermissionDenied("Underlying VCF data is archived")
    return trio


def _get_gene_region(trio: Trio, form: TrioKaryomappingGeneForm) -> Optional[TrioKaryomappingGeneRegion]:
    if not form.is_valid():
        return None
    gene_region = TrioKaryomappingGeneRegion(trio=trio, **form.cleaned_data)
    if gene_region.gene_version is None:
        form.add_error("gene_symbol", f"{gene_region.gene_symbol} has no gene in {trio.genome_build}")
        return None
    return gene_region


def get_variant_lookup_and_scatter_data(karyomapping_bins):
    """ Dumped to JS to be used by Plotly scatterplot
        karyomapping_bins : Have separate entries for ALT/REF, we merge these for output """

    variant_id_lookup = {}
    data = defaultdict(lambda: defaultdict(list))
    for karyotype_code, variant_data in karyomapping_bins.items():
        collapsed_code = KaryotypeBins.COLLAPSED_BINS[karyotype_code]
        for variant_id, _contig_id, chrom, position, ref, alt, svlen in variant_data:
            variant_string = Variant.format_tuple(chrom, position, ref, alt, svlen)
            variant_id_lookup[variant_string] = variant_id
            data[collapsed_code]["x"].append(position)
            data[collapsed_code]["text"].append(variant_string)

    karyotype_bin_counts = {kc: len(data[kc]["x"]) for kc in KaryotypeBins.KARYOTYPE_LABEL_ORDER}
    return variant_id_lookup, data, karyotype_bin_counts


def trio_karyomapping(request, trio_id):
    trio = _get_trio(request, trio_id)
    form = TrioKaryomappingGeneForm(request.GET or None)
    context = {"trio": trio, "form": form}
    if gene_region := _get_gene_region(trio, form):
        variant_id_lookup, scatter_data, bin_counts = get_variant_lookup_and_scatter_data(gene_region.get_karyomapping_bins())
        gene_version = gene_region.gene_version
        context.update({
            "gene_region": gene_region,
            "gene_version": gene_version,
            "karyotype_bin_labels": KaryotypeBins.KARYOTYPE_LABEL_ORDER,  # Want order
            "variant_id_lookup": variant_id_lookup,
            "karyotype_bin_scatter_data": scatter_data,
            "karyotype_bin_counts": bin_counts,
        })
    return render(request, 'analysis/karyomapping/trio_karyomapping.html', context)


def trio_karyomapping_gene_csv(request, trio_id):
    trio = _get_trio(request, trio_id)
    form = TrioKaryomappingGeneForm(request.GET)
    if not (gene_region := _get_gene_region(trio, form)):
        raise Http404(form.errors.as_text())

    variant_and_genotypes = gene_region.get_variant_and_genotypes()
    filename = (f"karyomapping_trio_{trio.pk}_{gene_region.gene_symbol}_"
                f"up_{gene_region.upstream_kb}kb_down_{gene_region.downstream_kb}kb.csv")

    header = ['chrom', 'position', 'ref', 'alt', 'svlen', 'proband_gt', 'father_gt', 'mother_gt', 'karyotype_bin']
    pseudo_buffer = StashFile()
    writer = csv.DictWriter(pseudo_buffer, header, dialect='excel')

    def iter_row_writer():
        writer.writeheader()
        yield pseudo_buffer.value
        for variant_data, genotype_tuple in variant_and_genotypes:
            _, _contig_id, chrom, position, ref, alt, svlen = variant_data
            proband_gt, father_gt, mother_gt = genotype_tuple
            row = {'chrom': chrom,
                   'position': position,
                   'ref': ref,
                   'alt': alt,
                   'svlen': svlen,
                   'proband_gt': Zygosity.get_genotype(proband_gt),
                   'father_gt': Zygosity.get_genotype(father_gt),
                   'mother_gt': Zygosity.get_genotype(mother_gt),
                   'karyotype_bin': KaryotypeBins.get_bin(genotype_tuple) or ''}
            writer.writerow(row)
            yield pseudo_buffer.value

    response = StreamingHttpResponse(iter_row_writer(), content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="{filename}"'
    return response
