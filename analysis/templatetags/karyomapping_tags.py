from django.core.exceptions import ObjectDoesNotExist
from django.template import Library

from analysis.models.models_karyomapping import KaryotypeBins
from snpdb.models import Trio

register = Library()


@register.inclusion_tag("analysis/tags/trio_karyomapping_qc.html")
def trio_karyomapping_qc(trio: Trio):
    """ Relatedness summary and per-chromosome phase graph from the trio's stored karyomapping counts - a parent
        well off 50% genome-wide suggests a sample swap, on one chromosome UPD or a deletion (#1754) """
    try:
        genome_karyomapping_counts = trio.genomekaryomappingcounts
    except ObjectDoesNotExist:
        genome_karyomapping_counts = None

    contig_qc = []
    if genome_karyomapping_counts:
        for contig_counts in genome_karyomapping_counts.get_contig_counts():
            collapsed_counts = contig_counts.get_collapsed_counts()
            father_percent, mother_percent = contig_counts.get_phase_percents()
            contig_qc.append({
                "contig": contig_counts.contig.name,
                "father_percent": father_percent,
                "mother_percent": mother_percent,
                "father_in_phase": collapsed_counts[KaryotypeBins.FATHER_1],
                "father_out_of_phase": collapsed_counts[KaryotypeBins.FATHER_2],
                "mother_in_phase": collapsed_counts[KaryotypeBins.MOTHER_1],
                "mother_out_of_phase": collapsed_counts[KaryotypeBins.MOTHER_2],
            })
    return {"trio": trio,
            "genome_karyomapping_counts": genome_karyomapping_counts,
            "collapsed_counts": genome_karyomapping_counts.get_collapsed_counts() if genome_karyomapping_counts else None,
            "contig_qc": contig_qc}
