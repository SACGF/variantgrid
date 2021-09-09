from collections import Counter, defaultdict

import celery

from analysis.models.models_karyomapping import KaryotypeBins, GenomeKaryomappingCounts, ContigKaryomappingCounts
from snpdb.models import Trio
from snpdb.models.models_enums import ImportStatus


@celery.shared_task
def create_genome_karyomapping_for_trio(trio_id):
    trio = Trio.objects.get(pk=trio_id)

    GenomeKaryomappingCounts.objects.filter(trio=trio).delete()
    genome_karyomapping_counts = GenomeKaryomappingCounts.objects.create(trio=trio)

    karyotype_bin_lookup = KaryotypeBins.get_karotype_bin_lookup()
    variant_and_genotypes = KaryotypeBins.get_variant_and_genotypes(trio)

    contig_code_count = defaultdict(Counter)
    for variant_data, genotype_tuple in variant_and_genotypes:
        contig_id = variant_data[1]
        (proband_gt, father_gt, mother_gt) = genotype_tuple

        try:
            karyotype_bin = karyotype_bin_lookup[proband_gt][father_gt][mother_gt]
            contig_code_count[contig_id][karyotype_bin] += 1
        except KeyError:
            pass

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
