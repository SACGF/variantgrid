from collections import Counter, defaultdict

import celery

from analysis.models.models_karyomapping import (
    ContigKaryomappingCounts,
    GenomeKaryomappingCounts,
    KaryotypeBins,
)
from snpdb.models import Trio
from snpdb.models.models_enums import ImportStatus


@celery.shared_task
def create_genome_karyomapping_for_trio(trio_id):
    trio = Trio.objects.get(pk=trio_id)

    GenomeKaryomappingCounts.objects.filter(trio=trio).delete()
    genome_karyomapping_counts = GenomeKaryomappingCounts.objects.create(trio=trio)

    contig_code_count = defaultdict(Counter)
    for variant_data, genotype_tuple in KaryotypeBins.iter_variant_and_genotypes(trio):
        if karyotype_bin := KaryotypeBins.get_bin(genotype_tuple):
            contig_id = variant_data[1]
            contig_code_count[contig_id][karyotype_bin] += 1

    genome_counts = Counter()
    contig_karyomapping_counts_list = []
    for contig_id, code_count in contig_code_count.items():
        genome_counts.update(code_count)

        kwargs = {"genome_karyomapping_counts": genome_karyomapping_counts,
                  "contig_id": contig_id}
        for code, count in code_count.items():
            kwargs[code.lower()] = count

        contig_karyomapping_counts_list.append(ContigKaryomappingCounts(**kwargs))

    ContigKaryomappingCounts.objects.bulk_create(contig_karyomapping_counts_list)

    for code, count in genome_counts.items():
        setattr(genome_karyomapping_counts, code.lower(), count)
    genome_karyomapping_counts.import_status = ImportStatus.SUCCESS
    genome_karyomapping_counts.save()
