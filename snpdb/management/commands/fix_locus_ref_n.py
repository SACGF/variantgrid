#         Dragen/CNV kit import variants with ref=N
#
#         We now change this using norm, but some existing ones are the DB
import logging
import os
import subprocess
from collections import defaultdict, Counter

import cyvcf2
from django.conf import settings
from django.core.management.base import BaseCommand
from django.db.models.query_utils import Q

from library.genomics.vcf_utils import get_contigs_header_lines, write_vcf_from_variant_coordinates
from library.utils import mk_path
from snpdb.models import GenomeBuild, Sequence, Locus, Variant, VariantCoordinate
from upload.models import ModifiedImportedVariant


class Command(BaseCommand):
    """
        Dragen/CNV kit import variants with ref=N

        We now change this using norm, but some existing ones are the DB

        We can't just look up the reference as it may need to be normalized

    """
    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **options):
        dry_run = options["dry_run"]

        seq_n = Sequence.objects.filter(seq='N').first()
        single_base_seq = {seq.seq: seq for seq in Sequence.objects.filter(seq__in='GATC')}

        for genome_build in GenomeBuild.builds_with_annotation():
            chrom_contig_id_mappings = genome_build.get_chrom_contig_id_mappings()

            # Collect all the coordinates, we are going to write them out and normalize them
            variant_and_coordinate_by_id = {}
            vcf_ids = []
            variant_coordinates = []

            logging.info("Finding variants with reference N")
            q_contig = Variant.get_contigs_q(genome_build)
            q_non_ref = Variant.get_no_reference_q()
            for variant in Variant.objects.filter(q_contig & q_non_ref, locus__ref=seq_n):
                vc = variant.coordinate.as_contig_accession(genome_build)
                vcf_ids.append(variant.pk)
                variant_coordinates.append(vc)
                variant_and_coordinate_by_id[variant.pk] = variant, vc

            logging.info("Found %d variants with ref N", len(vcf_ids))

            # Write out
            change_count = Counter()
            if vcf_ids:
                processing_dir = os.path.join(settings.PRIVATE_DATA_ROOT, 'fix_variant_ref_n')
                mk_path(processing_dir)
                vcf_input_filename = os.path.join(processing_dir, f"old_variants_ref_n_{genome_build.name}.vcf")
                used_chroms = set((vc.chrom for vc in variant_coordinates))
                header_lines = get_contigs_header_lines(genome_build, use_accession=True,
                                                        contig_allow_list=used_chroms)
                logging.info("Writing old variants ref n VCF: %s", vcf_input_filename)
                write_vcf_from_variant_coordinates(vcf_input_filename, variant_coordinates=variant_coordinates,
                                                   vcf_ids=vcf_ids, header_lines=header_lines)

                vcf_output_filename = os.path.join(processing_dir, f"normalized_ref_n_{genome_build.name}.vcf")
                bcftools_normalize = [
                    settings.BCFTOOLS_COMMAND, "norm",
                    "--check-ref=s",  # Set ref (ie replace N with actual ref base)
                    f"--old-rec-tag={ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG}",
                    f"--fasta-ref={genome_build.reference_fasta}",
                    "--output", vcf_output_filename,
                    vcf_input_filename,
                ]

                logging.info("Running bcftools normalize")
                subprocess.check_call(bcftools_normalize)
                logging.info("Wrote normalized VCF: %s", vcf_output_filename)

                reader = cyvcf2.Reader(vcf_output_filename)
                modified_imported_variants = []
                loci_update_ref = []
                loci_old_new = {}

                for record in reader:
                    logging.info("alt: %s", record.ALT)
                    variant, vc = variant_and_coordinate_by_id[int(record.ID)]
                    svlen = record.INFO.get("SVLEN")
                    normalized_vc = VariantCoordinate(chrom=record.CHROM, position=record.POS,
                                                      ref=record.REF, alt=record.ALT[0],
                                                      svlen=svlen)

                    # Verify these haven't changed
                    for field in ["chrom", "alt", "svlen"]:
                        old_val = getattr(vc, field)
                        new_val = getattr(vc, field)
                        if old_val != new_val:
                            raise ValueError(f"{record.ID=}: {field} changed after normalization from {old_val=} to {new_val=}")

                    logging.info("old: %s", vc)
                    logging.info("new: %s", normalized_vc)
                    locus = variant.locus
                    ref = single_base_seq[normalized_vc.ref.upper()]
                    contig_id = chrom_contig_id_mappings[normalized_vc.chrom]

                    # if contig/position is the same - can update locus
                    if (vc.chrom, vc.position) == (normalized_vc.chrom, normalized_vc.position):
                        if existing_locus := Locus.objects.filter(contig_id=contig_id,
                                                                  position=normalized_vc.position,
                                                                  ref=ref).first():
                            change_count["change to existing locus"] += 1
                            loci_old_new[locus] = existing_locus
                        else:
                            change_count["replace locus ref"] += 1
                            # Update this locus to be the new one
                            locus.ref = ref
                            loci_update_ref.append(locus)
                    else:
                        bcftools_old_variant = variant.INFO[ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG]
                        for ov in ModifiedImportedVariant.bcftools_format_old_variant(bcftools_old_variant, svlen, genome_build):
                            miv = ModifiedImportedVariant(variant=variant,
                                                          old_variant=bcftools_old_variant,
                                                          old_variant_formatted=ov)
                            modified_imported_variants.append(miv)

                        # Shifted position
                        new_locus, existing = Locus.objects.get_or_create(contig_id=contig_id, position=normalized_vc.position, ref=ref)
                        loci_old_new[locus] = new_locus
                        if existing:
                            change_count["shift to existing locus"] += 1
                        else:
                            change_count["shift to new locus"] += 1

                logging.info("Updates: %s", change_count)
                if loci_update_ref:
                    logging.info("%s: Updating %d ref sequences", genome_build, len(loci_update_ref))
                    if not dry_run:
                        Locus.objects.bulk_update(loci_update_ref, ['ref'], batch_size=2000)

                if loci_old_new:
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
                        if not dry_run:
                            Variant.objects.bulk_update(variant_update_locus, ['locus'], batch_size=256)

                        if modified_imported_variants:
                            logging.info("Inserting MIV for normalized variants: %d", len(modified_imported_variants))
                            ModifiedImportedVariant.objects.bulk_create(modified_imported_variants)

                        if not dry_run:
                            logging.info("Deleting unused loci")
                            q_locus_contig = Q(contig__genomebuildcontig__genome_build=genome_build)
                            Locus.objects.filter(q_locus_contig, ref=seq_n, variant__isnull=True).delete()



                    if variant_not_updated:
                        logging.info("%s: Variants that could not be updated due to dupes: %d", genome_build, len(variant_not_updated))
                        logging.info("%s: It may be easiest to re-import variants that caused these.", genome_build)
                        if len(variant_not_updated) < 50:
                            variants_str = ",".join([str(v.pk) for v in variant_not_updated])
                            logging.info("%s: %s",genome_build, variants_str)



