"""
Variant table is massive, so we need a way to quickly look up a variant PK for insertion

The original implementation RedisVariantPKLookup used Redis to store a hash, was deleted 2021-08-05
"""
import logging
import os
from collections import defaultdict
from typing import Iterable, TypeAlias, Tuple, Callable, Collection, Any

from django.db.models import Q, Value, TextField, QuerySet
from django.db.models.aggregates import Max
from django.db.models.functions import Concat

from library.utils import sha256sum_str
from snpdb.models import Variant, Locus, Sequence, GenomeBuild, VariantCoordinate
from upload.vcf import sql_copy_files

VariantHash: TypeAlias = Tuple[int, int, int, int, int | str]
LocusHash: TypeAlias = Tuple[int, int, int]
LociHash: TypeAlias = Any
AnyHash: TypeAlias = VariantHash | LocusHash | LociHash


class VariantPKLookup:
    """ variant hash is something we can compute via SQL from the Variant """

    def __init__(self, genome_build: GenomeBuild = None, working_dir=None):
        self.genome_build: GenomeBuild = genome_build
        self._working_dir: str = working_dir
        self._file_batch_id = 0
        self.variant_hashes: list[str] = []  # Variant hash and coordinate lists are kept in sync
        self.variant_coordinates: list[VariantCoordinate] = []
        self.variant_pk_by_hash: dict[VariantHash, int] = {}
        self.unknown_variant_coordinates: list[VariantCoordinate] = []

        if genome_build:
            self.chrom_contig_id_mappings = genome_build.get_chrom_contig_id_mappings()
        else:
            self.chrom_contig_id_mappings = None
        self.sequence_pk_by_seq = Sequence.get_pk_by_seq()
        defaults = {
            "seq_sha256_hash": sha256sum_str(Variant.REFERENCE_ALT),
        }
        self.reference_seq_id = Sequence.objects.get_or_create(seq=Variant.REFERENCE_ALT,
                                                               defaults=defaults)[0].pk

    def _get_locus_hash(self, contig_id, position, ref_id):
        return contig_id, position, ref_id

    def _get_variant_hash(self, contig_id: int, position: int, ref_id: int, alt_id: int, svlen: int) -> VariantHash:
        return contig_id, position, ref_id, alt_id, svlen or ''

    def _get_ids_for_hashes(self, hashes: Iterable[AnyHash], get_queryset: Callable[[int, Collection[int]], QuerySet]) -> list[int]:
        """ hashes tuples 1st 2 elements should be contig_id, position """
        contig_positions = defaultdict(set)
        for contig_id, position, *_ in hashes:
            contig_positions[contig_id].add(position)

        contig_hashes: dict[int, dict] = defaultdict(dict)
        for contig_id, positions in contig_positions.items():
            qs = get_queryset(contig_id, positions)
            contig_hashes[contig_id] = dict(qs.values_list("hash", "id"))

        return [contig_hashes[h[0]].get("_".join([str(s) for s in h[1:]])) for h in hashes]

    def _get_variant_ids(self, variant_hashes: Iterable[VariantHash]) -> list[int]:
        annotate_kwargs = {
            "hash": Concat("locus__position", Value("_"), "locus__ref_id", Value("_"), "alt_id",
                           Value("_"), "svlen",
                           output_field=TextField())
        }

        def get_queryset(contig_id, positions) -> QuerySet:
            qs = Variant.objects.filter(locus__contig_id=contig_id, locus__position__in=positions)
            return qs.annotate(**annotate_kwargs)

        return self._get_ids_for_hashes(variant_hashes, get_queryset)

    def get_variant_ids(self, variant_hashes: Iterable[VariantHash], validate_not_null=True) -> list[int]:
        variant_ids = self._get_variant_ids(variant_hashes)
        if validate_not_null:
            if not all(variant_ids):  # Quick test if ok
                # Slowly find the first bad record
                for variant_hash, pk in zip(variant_hashes, variant_ids):
                    if not pk:
                        raise ValueError(f"Variant hash {variant_hash} had no PK in DB!")
        return variant_ids

    def _get_loci_ids(self, loci_hashes: Iterable) -> list[int]:
        annotate_kwargs = {
            "hash": Concat("position", Value("_"), "ref_id",
                           output_field=TextField())
        }

        def get_queryset(contig_id, positions):
            qs = Locus.objects.filter(contig_id=contig_id, position__in=positions)
            return qs.annotate(**annotate_kwargs)

        return self._get_ids_for_hashes(loci_hashes, get_queryset)

    def get_loci_ids(self, loci_hashes: Iterable, validate_not_null=True) -> list[int]:
        loci_ids = self._get_loci_ids(loci_hashes)
        if validate_not_null:
            if not all(loci_ids):  # Quick test if ok
                # Slowly find the first bad record
                for loci_hash, pk in zip(loci_hashes, loci_ids):
                    if not pk:
                        raise ValueError(f"Loci hash {loci_hash} had no PK in DB!")
        return loci_ids

    def filter_non_reference(self, variant_hashes: Iterable[VariantHash], variant_ids: Iterable[int]) -> list[int]:
        return [vh_vi[1] for vh_vi in zip(variant_hashes, variant_ids) if vh_vi[0][3] != self.reference_seq_id]

    def get_variant_coordinate_hash(self, variant_coordinate: VariantCoordinate) -> VariantHash:
        """ For VCF records (needs GenomeBuild supplied) """

        variant_coordinate = variant_coordinate.as_internal_canonical_form(self.genome_build)
        if self.chrom_contig_id_mappings is None:
            raise ValueError("Need to initialise w/GenomeBuild to call get_variant_coordinate_hash")
        contig_id = self.chrom_contig_id_mappings[variant_coordinate.chrom]
        ref_id = self.sequence_pk_by_seq[variant_coordinate.ref]
        alt_id = self.sequence_pk_by_seq[variant_coordinate.alt]
        return self._get_variant_hash(contig_id, variant_coordinate.position, ref_id, alt_id, variant_coordinate.svlen)

    def add(self, variant_coordinate: VariantCoordinate) -> VariantHash:
        variant_coordinate = variant_coordinate.as_internal_canonical_form(self.genome_build)

        # If sequence isn't known, variant is definitely unknown
        if variant_coordinate.ref in self.sequence_pk_by_seq and variant_coordinate.alt in self.sequence_pk_by_seq:
            # Maybe unknown, need to check
            variant_hash = self.get_variant_coordinate_hash(variant_coordinate)
            self.variant_hashes.append(variant_hash)
            self.variant_coordinates.append(variant_coordinate)
        else:
            variant_hash = None
            self.unknown_variant_coordinates.append(variant_coordinate)
        return variant_hash

    def batch_check(self, minimum_pipeline_size=0, insert_unknown=False):
        # Check if known/unknown
        # Don't do this in a separate thread, as we need to ensure this gets back to populate unknown_variants
        # during final finish call of batch_process_check(insert_all=True)
        num_variant_hashes = len(self.variant_hashes)
        if num_variant_hashes and num_variant_hashes >= minimum_pipeline_size:
            variant_ids = self.get_variant_ids(self.variant_hashes, validate_not_null=False)
            for variant_hash, variant_id, variant_coordinate in zip(self.variant_hashes, variant_ids,
                                                                    self.variant_coordinates):
                if variant_id:
                    self.variant_pk_by_hash[variant_hash] = variant_id
                else:
                    self.unknown_variant_coordinates.append(variant_coordinate)

            self.variant_hashes = []
            self.variant_coordinates = []

        num_unknown = len(self.unknown_variant_coordinates)
        if insert_unknown and num_unknown >= minimum_pipeline_size:
            self._insert_unknown()

    def clear(self):
        self.variant_pk_by_hash.clear()
        self.unknown_variant_coordinates = []

    def _insert_unknown_sequences(self):
        """ Inserts unknown from coordinates, updates sequence_pk_by_seq """
        unknown_sequences = set()
        for variant_coordinate in self.unknown_variant_coordinates:
            if variant_coordinate.ref not in self.sequence_pk_by_seq:
                unknown_sequences.add(variant_coordinate.ref)
            if variant_coordinate.alt not in self.sequence_pk_by_seq:
                unknown_sequences.add(variant_coordinate.alt)

        if unknown_sequences:
            # bulk create not returning pks, so have to retrieve any new ones ourselves
            max_dict = Sequence.objects.all().aggregate(highest_pk=Max("pk"))
            old_highest_pk = max_dict["highest_pk"] or 0
            sequences = [Sequence(seq=seq, seq_sha256_hash=sha256sum_str(seq)) for seq in unknown_sequences]
            Sequence.objects.bulk_create(sequences, ignore_conflicts=True)
            self.sequence_pk_by_seq.update(Sequence.get_pk_by_seq(Q(pk__gt=old_highest_pk)))

    def _get_csv_filename(self, prefix) -> str:
        return os.path.join(self._working_dir, f"{prefix}_{self._file_batch_id}.csv")

    def _insert_new_loci(self, new_loci_rows: list[tuple]):
        logging.info("1st new_loci_rows: %r", str(new_loci_rows[0]))
        unknown_loci_filename = self._get_csv_filename("unknown_loci")
        sql_copy_files.write_sql_copy_csv(new_loci_rows, unknown_loci_filename)
        logging.info("wrote loci file: %s", unknown_loci_filename)
        sql_copy_files.loci_sql_copy_csv(unknown_loci_filename)

    def _insert_new_variants(self, new_variant_rows: list[tuple]):
        logging.info("1st unknown variant: locus_id, alt_id, end : %s", str(new_variant_rows[0]))
        unknown_variants_filename = self._get_csv_filename("unknown_variants")
        sql_copy_files.write_sql_copy_csv(new_variant_rows, unknown_variants_filename)
        logging.info("wrote variants file: %s", unknown_variants_filename)
        sql_copy_files.variants_sql_copy_csv(unknown_variants_filename)

    def _insert_unknown(self):
        self._insert_unknown_sequences()  # All ref/alt need to be in sequence_pk_by_seq to be able to get hash

        # We need these to be ordered, so will use dict keys
        loci_hashes = {}  # Uses None as values, just want orderable
        variant_coordinates_by_hash = {}

        for variant_coordinate in self.unknown_variant_coordinates:
            variant_hash = self.get_variant_coordinate_hash(variant_coordinate)
            variant_coordinates_by_hash[variant_hash] = variant_coordinate

            (contig_id, start, ref_id, _alt_id, _svlen) = variant_hash
            locus_hash = (contig_id, start, ref_id)
            loci_hashes[locus_hash] = None

        self.unknown_variant_coordinates.clear()
        loci_ids = self.get_loci_ids(loci_hashes, validate_not_null=False)  # Could have some unknowns as None
        unknown_loci_parts = []
        for locus_hash, locus_pk in zip(loci_hashes, loci_ids):
            if locus_pk is None:
                unknown_loci_parts.append(locus_hash)

        if unknown_loci_parts:
            self._insert_new_loci(unknown_loci_parts)

        loci_ids = self.get_loci_ids(loci_hashes)  # Validates all are not null
        locus_pk_by_hash = dict(zip(loci_hashes, loci_ids))
        variant_ids = self.get_variant_ids(variant_coordinates_by_hash, validate_not_null=False)
        unknown_variants = []
        for variant_hash, variant_pk in zip(variant_coordinates_by_hash, variant_ids):
            if variant_pk is None:
                (contig_id, start, ref_id, alt_id, svlen) = variant_hash
                locus_hash = (contig_id, start, ref_id)
                locus_pk = locus_pk_by_hash[locus_hash]
                vc = variant_coordinates_by_hash[variant_hash]
                unknown_variants.append((locus_pk, alt_id, vc.end, svlen))
        logging.debug("loci_hashes = %d, unknown_variants = %d", len(loci_hashes), len(unknown_variants))

        if unknown_variants:
            self._insert_new_variants(unknown_variants)

        self._file_batch_id += 1
