import logging
from collections import defaultdict

from django.core.management.base import BaseCommand
from django.db.models import Q
from snpdb.models import Variant, Sequence, GenomeBuild, Locus


class Command(BaseCommand):
    """
        Dragen/CNV kit import variants with ref=N

        We now change this using norm, but some existing ones are the DB

        To do this, we can:

        * If a locus has a N, and it's the only one, switch it out
        * Otherwise, switch all variants linking to it to the other locus.
        * If this fails, then
        * Then for Variants - switch to the switch locus. If this fails, it's because there is an existing variant with the new locus
        *

    """
    def handle(self, *args, **options):
        seq_n = Sequence.objects.filter(seq='N').first()
        single_base_seq = {seq.seq: seq for seq in Sequence.objects.filter(seq__in='GATC')}

        for genome_build in GenomeBuild.builds_with_annotation():
            # If loci with correct base exists, we want to swap all variants to that one
            loci_old_new = {}
            loci_update_ref = []  # For loci we can just switch out
            num_ref_actually_n = 0
            q_contig = Q(contig__genomebuildcontig__genome_build=genome_build)
            for locus in Locus.objects.filter(q_contig, ref=seq_n):
                contig_sequence = genome_build.genome_fasta.fasta[locus.chrom]
                # reference sequence is 0-based, we only want 1 base
                position = locus.position - 1
                ref_sequence = contig_sequence[position:position+1].upper()
                if ref_sequence != 'N':
                    if existing_locus := Locus.objects.filter(contig=locus.contig, position=locus.position, ref__seq=ref_sequence).first():
                        loci_old_new[locus] = existing_locus
                    else:
                        # Update this locus to be the new one
                        locus.ref = single_base_seq[ref_sequence]
                        loci_update_ref.append(locus)
                else:
                    num_ref_actually_n += 1

            logging.info("%s: search found num_ref_actually_n=%d, loci to update ref=%d, switch=%d",
                         genome_build, num_ref_actually_n, len(loci_update_ref), len(loci_old_new))
            if loci_update_ref:
                logging.info("%s: Updating %d ref sequences", genome_build, len(loci_update_ref))
                Locus.objects.bulk_update(loci_update_ref, ['ref'], batch_size=2000)

            if loci_old_new:
                # loci_old_new

                # Now we need to update the variants that use these old loci to the new ones
                # Some would fail as there are already existing variants with the new loci breaking constraint
                # There is no ignore conflict on bulk update so we will need to get these out
                new_loci_variants = defaultdict(set)
                for v in Variant.objects.filter(locus__in=loci_old_new.values()):
                    new_loci_variants[v.locus].add(v)

                variant_not_updated = []
                variant_update_locus = []
                for variant in Variant.objects.filter(locus__in=loci_old_new.keys()):
                    new_locus = loci_old_new[variant.locus]
                    conflicted_variant = None
                    for potential_conflict in new_loci_variants[new_locus]:
                        if variant.alt == potential_conflict.alt and variant.svlen == potential_conflict.svlen:
                            conflicted_variant = potential_conflict
                            break

                    if conflicted_variant:
                        variant_not_updated.append(variant)
                    else:
                        variant.locus = new_locus
                        variant_update_locus.append(variant)

                if variant_update_locus:
                    Variant.objects.bulk_update(variant_update_locus, ['locus'], batch_size=256)

                    # TODO: Delete unused loci now (new QS so it isn't cached)
                    logging.info("Deleting unused loci")
                    Locus.objects.filter(q_contig, ref=seq_n, variant__isnull=True).delete()

                if variant_not_updated:
                    logging.info("%s: Variants that could not be updated due to dupes: %d", genome_build, len(variant_not_updated))
                    logging.info("%s: It may be easiest to re-import variants that caused these.", genome_build)
                    if len(variant_not_updated) < 50:
                        variants_str = ",".join([str(v.pk) for v in variant_not_updated])
                        logging.info("%s: %s",genome_build, variants_str)
