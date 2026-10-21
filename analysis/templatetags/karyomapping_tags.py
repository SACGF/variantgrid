from django.template import Library

from analysis.models.models_karyomapping import KaryotypeBins
from analysis.tasks.karyomapping_tasks import get_or_queue_genome_karyomapping_counts
from snpdb.enums import ImportStatus
from snpdb.models import Trio

register = Library()


@register.inclusion_tag("analysis/tags/trio_karyomapping_qc.html")
def trio_karyomapping_qc(trio: Trio):
    """ Relatedness summary and per-chromosome phase graph from the trio's stored karyomapping counts - a parent
        well off 50% genome-wide suggests a sample swap, on one chromosome UPD or a deletion (#1754) """
    genome_karyomapping_counts = get_or_queue_genome_karyomapping_counts(trio)
    calculated = genome_karyomapping_counts.import_status == ImportStatus.SUCCESS

    contig_qc = []
    if calculated:
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
            "calculated": calculated,
            "collapsed_counts": genome_karyomapping_counts.get_collapsed_counts() if calculated else None,
            "contig_qc": contig_qc}
