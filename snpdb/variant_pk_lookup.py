"""
Variant table is massive, so we need a way to quickly look up a variant PK for insertion
"""

import logging
import os
from time import time
from typing import List

from django.conf import settings
from django.db.models import Q
from django.db.models.aggregates import Max

from library.django_utils import get_redis
from library.utils import md5sum_str
from snpdb.models import Variant, Locus, Sequence, GenomeBuild, ProcessingStatus, VariantCoordinate
from upload.models import UploadStep
from upload.vcf import sql_copy_files


class VariantPKLookup:
    MAX_LOCI_ID = 'MAX_LOCI_ID'
    MAX_VARIANT_ID = 'MAX_VARIANT_ID'

    def __init__(self, genome_build: GenomeBuild = None, redis_check=True, working_dir=None):
        self.genome_build = genome_build
        self._working_dir = working_dir
        self._file_batch_id = 0
        self.redis = get_redis()
        self.variant_hashes: List[str] = []  # Variant hash and coordinate lists are kept in sync
        self.variant_coordinates: List[VariantCoordinate] = []
        self.variant_pk_by_hash = {}
        self.unknown_variant_coordinates: List[VariantCoordinate] = []

        if genome_build:
            self.chrom_contig_id_mappings = genome_build.get_chrom_contig_id_mappings()
        else:
            self.chrom_contig_id_mappings = None
        self.sequence_pk_by_seq = Sequence.get_pk_by_seq()
        if redis_check:
            self.redis_check = self.redis_sanity_check()

    def redis_sanity_check(self) -> bool:
        """ Check that the first and last Variant is stored in Redis, throws exception on error
            Returns whether it was able to perform a check or not """

        running_insert_jobs = UploadStep.objects.filter(name=UploadStep.CREATE_UNKNOWN_LOCI_AND_VARIANTS_TASK_NAME,
                                                        status=ProcessingStatus.PROCESSING)
        if running_insert_jobs.exists():
            logging.warning("Unknown variant insert job running - unable to perform Redis sanity check")
            return False

        qs = Variant.objects.order_by("pk")
        first_variant = qs.first()
        last_variant = qs.last()

        for v in [first_variant, last_variant]:
            if v is None:
                continue
            v_hash = self.get_variant_object_hash(v)
            variant_id = self.redis.get(v_hash)
            if variant_id is None:
                missing_key_msg = f"Redis sanity check failed - variant: {v.pk} ({v_hash}) missing!"
                logging.error(missing_key_msg)
                raise ValueError(missing_key_msg)
            if v.pk != int(variant_id):
                diff_key_msg = f"Redis sanity check failed - hash: {v_hash} value from variant {v.pk} != '{variant_id}'"
                logging.error(diff_key_msg)
                raise ValueError(diff_key_msg)
        return True  # Checked ok

    @staticmethod
    def _get_locus_hash(contig_id, position, ref_id):
        return '_'.join((str(contig_id), str(position), str(ref_id)))

    @staticmethod
    def _get_variant_hash(contig_id, position, ref_id, alt_id) -> str:
        locus_hash = VariantPKLookup._get_locus_hash(contig_id, position, ref_id)
        return f"{locus_hash}_{alt_id}"

    @staticmethod
    def get_variant_object_hash(v: Variant):
        return VariantPKLookup._get_variant_hash(v.locus.contig_id, v.locus.position, v.locus.ref_id, v.alt_id)

    def get_variant_coordinate_hash(self, chrom, position, ref, alt):
        """ For VCF records (needs GenomeBuild supplied) """
        if self.chrom_contig_id_mappings is None:
            raise ValueError("Need to initialise w/GenomeBuild to call get_variant_coordinate_hash")
        contig_id = self.chrom_contig_id_mappings[chrom]
        ref_id = self.sequence_pk_by_seq[ref]
        alt_id = self.sequence_pk_by_seq[alt]
        return self._get_variant_hash(contig_id, position, ref_id, alt_id)

    def add(self, chrom, position, ref, alt):
        variant_coordinate = (chrom, position, ref, alt)
        # If sequence isn't known, variant is definitely unknown
        if ref in self.sequence_pk_by_seq and alt in self.sequence_pk_by_seq:
            # Maybe unknown, check redis
            variant_hash = self.get_variant_coordinate_hash(*variant_coordinate)
            self.variant_hashes.append(variant_hash)
            self.variant_coordinates.append(variant_coordinate)
        else:
            variant_hash = None
            self.unknown_variant_coordinates.append(variant_coordinate)
        return variant_hash

    def get_variant_ids(self, variant_hashes: List[str]):
        return self.redis.mget(variant_hashes)

    def batch_check(self, minimum_redis_pipeline_size=0, insert_unknown=False):
        # Check if known/unknown via Redis
        # Don't do this in a separate thread, as we need to ensure this gets back to populate unknown_variants
        # during final finish call of batch_process_check(insert_all=True)
        num_variant_hashes = len(self.variant_hashes)
        if num_variant_hashes and num_variant_hashes >= minimum_redis_pipeline_size:
            variant_ids = self.get_variant_ids(self.variant_hashes)
            for variant_hash, variant_id, variant_coordinate in zip(self.variant_hashes, variant_ids, self.variant_coordinates):
                if variant_id:
                    self.variant_pk_by_hash[variant_hash] = variant_id
                else:
                    self.unknown_variant_coordinates.append(variant_coordinate)

            self.variant_hashes = []
            self.variant_coordinates = []

        num_unknown = len(self.unknown_variant_coordinates)
        if insert_unknown and num_unknown >= minimum_redis_pipeline_size:
            self.insert_unknown()

    def clear(self):
        self.variant_pk_by_hash.clear()
        self.unknown_variant_coordinates = []

    # Note we use setnx - to not overwrite any variants/loci
    # This means that if there's ever a duplicated locus/variant from some race condition error,
    # it will never be used (lower primary key will be kept as insert in pk order).
    # See https://redis.io/commands/setnx
    def update_redis_hashes(self):
        """ This updates both loci and variants  """
        print("update_redis_hashes....")
        self._update_loci_hash()
        self._update_variants_hash()

    def _update_variants_hash(self):
        max_variant_id = int(self.redis.get(self.MAX_VARIANT_ID) or 0)

        old_max_variant_id = max_variant_id
        max_dict = Variant.objects.all().aggregate(Max("id"))
        max_variant_id = max_dict["id__max"] or 0
        num_potential_new_variants = max_variant_id - old_max_variant_id

        if num_potential_new_variants:
            start = time()
            logging.debug("Updating approx ~%d variants.", num_potential_new_variants)

            variants_qs = Variant.objects.all()

            if old_max_variant_id:
                variants_qs = variants_qs.filter(id__gt=old_max_variant_id)

            values_qs = variants_qs.values_list('id', 'locus__id', 'locus__contig', 'locus__position',
                                                'locus__ref', 'alt')

            i = 0
            redis_pipeline = self.redis.pipeline()
            inserted = 0

            # Count takes ~1% of the time and nice to give a running progress
            num_to_insert = values_qs.count()
            for variant_id, locus_id, contig_id, position, ref_id, alt_id in values_qs.order_by("pk").iterator():
                locus_hash = self._get_locus_hash(contig_id, position, ref_id)

                redis_pipeline.setnx(locus_hash, locus_id)
                i += 1

                if variant_id:
                    variant_hash = self._get_variant_hash(contig_id, position, ref_id, alt_id)
                    redis_pipeline.setnx(variant_hash, variant_id)

                if i >= settings.REDIS_PIPELINE_SIZE:
                    redis_pipeline.execute()
                    redis_pipeline = self.redis.pipeline()
                    inserted += i
                    i = 0

                    percent_done = 100.0 * inserted / num_to_insert
                    logging.debug("Inserted %d of %d (%.2f%%)", inserted, num_to_insert, percent_done)

            if i:
                logging.debug("Inserting remaining %d", i)
                redis_pipeline.execute()
                inserted += i

            end = time()
            time_taken = end - start
            logging.debug("Inserted %d of %d - Update took %.3f seconds", inserted, num_to_insert, time_taken)

        self.redis.set(self.MAX_VARIANT_ID, max_variant_id)

    def _update_loci_hash(self):
        """ Updates Just the loci (for when you have inserted loci but not variants yet) """

        max_loci_id = int(self.redis.get(self.MAX_LOCI_ID) or 0)
        old_max_loci_id = max_loci_id
        max_dict = Locus.objects.all().aggregate(Max("id"))
        max_loci_id = max_dict["id__max"] or 0
        num_new_loci = max_loci_id - old_max_loci_id

        if not num_new_loci:  # Look and double check
            if last_locus := Locus.objects.order_by("pk").last():
                locus_hash = self._get_locus_hash(last_locus.contig_id, last_locus.position, last_locus.ref_id)
                locus_pk = self.redis.get(locus_hash)
                if locus_pk is None:
                    num_new_loci = "unknown...."

        if num_new_loci:
            start = time()
            i = 0
            redis_pipeline = self.redis.pipeline()
            logging.debug("Updating %d loci.", num_new_loci)

            loci_values_qs = Locus.objects.values_list('id', 'contig', 'position', 'ref')
            if old_max_loci_id:
                loci_values_qs = loci_values_qs.filter(id__gt=old_max_loci_id)

            for locus_id, contig_id, position, ref_id in loci_values_qs.order_by("pk").iterator():
                locus_hash = self._get_locus_hash(contig_id, position, ref_id)

                redis_pipeline.setnx(locus_hash, locus_id)

                i += 1
                if i >= settings.REDIS_PIPELINE_SIZE:
                    logging.info("Pipeline execute (%d)", settings.REDIS_PIPELINE_SIZE)
                    redis_pipeline.execute()
                    redis_pipeline = self.redis.pipeline()
                    i = 0

            if i:
                redis_pipeline.execute()

            end = time()
            time_taken = end - start
            logging.debug("Update took %.1f seconds", time_taken)

        self.redis.set(self.MAX_LOCI_ID, max_loci_id)

    def _insert_unknown_sequences(self):
        """ Inserts unknown from coordinates, updates sequence_pk_by_seq """
        unknown_sequences = set()
        for chrom, position, ref, alt in self.unknown_variant_coordinates:
            if ref not in self.sequence_pk_by_seq:
                unknown_sequences.add(ref)
            if alt not in self.sequence_pk_by_seq:
                unknown_sequences.add(alt)

        if unknown_sequences:
            # bulk create not returning pks, so have to retrieve any new ones ourselves
            max_dict = Sequence.objects.all().aggregate(highest_pk=Max("pk"))
            old_highest_pk = max_dict["highest_pk"] or 0
            sequences = [Sequence(seq=seq, seq_md5_hash=md5sum_str(seq), length=len(seq)) for seq in unknown_sequences]
            Sequence.objects.bulk_create(sequences, ignore_conflicts=True)
            self.sequence_pk_by_seq.update(Sequence.get_pk_by_seq(Q(pk__gt=old_highest_pk)))

    def _get_csv_filename(self, prefix) -> str:
        return os.path.join(self._working_dir, f"{prefix}_{self._file_batch_id}.csv")

    def insert_unknown(self):
        self._insert_unknown_sequences()  # All ref/alt need to be in sequence_pk_by_seq to be able to get hash

        loci_parts_by_hash = {}
        locus_hash_and_alt_id_by_variant_hash = {}
        for chrom, position, ref, alt in self.unknown_variant_coordinates:
            contig_id = self.chrom_contig_id_mappings[chrom]
            ref_id = self.sequence_pk_by_seq[ref]
            alt_id = self.sequence_pk_by_seq[alt]
            loci_parts = (contig_id, position, ref_id)
            locus_hash = self._get_locus_hash(*loci_parts)
            loci_parts_by_hash[locus_hash] = loci_parts
            variant_parts = (contig_id, position, ref_id, alt_id)
            variant_hash = self._get_variant_hash(*variant_parts)
            locus_hash_and_alt_id_by_variant_hash[variant_hash] = (locus_hash, alt_id)

        self.unknown_variant_coordinates.clear()
        loci_ids = self.redis.mget(loci_parts_by_hash)  # unknowns will be None
        unknown_loci_parts = []
        for locus_hash, locus_pk in zip(loci_parts_by_hash, loci_ids):
            if locus_pk is None:
                unknown_loci_parts.append(loci_parts_by_hash[locus_hash])

        if unknown_loci_parts:
            logging.info("1st unknown_loci_parts: %r", str(unknown_loci_parts[0]))
            unknown_loci_filename = self._get_csv_filename("unknown_loci")
            sql_copy_files.write_sql_copy_csv(unknown_loci_parts, unknown_loci_filename)
            logging.info("wrote loci file: %s", unknown_loci_filename)
            sql_copy_files.loci_sql_copy_csv(unknown_loci_filename)
            self._update_loci_hash()  # Read in new loci

        loci_ids = self.redis.mget(loci_parts_by_hash)  # should all be there now
        locus_pk_by_hash = dict(zip(loci_parts_by_hash, loci_ids))
        variant_ids = self.redis.mget(locus_hash_and_alt_id_by_variant_hash)  # unknowns will be None
        unknown_variants = []
        for variant_hash, variant_pk in zip(locus_hash_and_alt_id_by_variant_hash, variant_ids):
            if variant_pk is None:
                (locus_hash, alt_id) = locus_hash_and_alt_id_by_variant_hash[variant_hash]
                locus_pk = locus_pk_by_hash[locus_hash]
                unknown_variants.append((locus_pk, alt_id))
        logging.debug("loci_hashes = %d, unknown_variants = %d", len(loci_parts_by_hash), len(unknown_variants))

        if unknown_variants:
            logging.info("1st unknown_variant_locus_id_and_alt: %s", str(unknown_variants[0]))
            unknown_variants_filename = self._get_csv_filename("unknown_variants")
            sql_copy_files.write_sql_copy_csv(unknown_variants, unknown_variants_filename)
            logging.info("wrote variants file: %s", unknown_variants_filename)
            sql_copy_files.variants_sql_copy_csv(unknown_variants_filename)
            self._update_variants_hash()  # Read in new variants

        self._file_batch_id += 1
