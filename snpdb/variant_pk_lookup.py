"""
Variant table is massive, so we need a way to quickly look up a variant PK for insertion

An explicit set of sequences (every GATCN sequence up to SHORT_SEQUENCE_MAX_LENGTH, the reference alt and the
supported symbolic alts) is looked up by hash up front, so a preloaded REF/ALT missing from the map means the variant
is new. That covers ~98% of variants in 1 indexed query, whatever the size of the Sequence table. It is an enumeration
rather than a length cutoff because the hash is the only index on Sequence - any other sequence stays as a string in
the hash, and is resolved to a pk in 1 query per batch (#2538)

The original implementation RedisVariantPKLookup used Redis to store a hash, was deleted 2021-08-05
"""
import itertools
import logging
import os
from collections import defaultdict
from collections.abc import Callable, Collection, Iterable
from functools import cache, cached_property
from typing import NamedTuple, Optional, TypeAlias

from django.conf import settings
from django.db.models import QuerySet, TextField, Value
from django.db.models.functions import Concat

from library.utils import sha256sum_str
from snpdb.models import GenomeBuild, Locus, Sequence, Variant, VariantCoordinate
from upload.vcf import sql_copy_files

# ref/alt are a Sequence pk for preloaded sequences, or the sequence itself for the rest
SequenceKey: TypeAlias = int | str


def _resolve_sequence_key(key: SequenceKey, sequence_pk_by_seq: dict[str, Optional[int]]) -> Optional[int]:
    return sequence_pk_by_seq[key] if isinstance(key, str) else key


class LocusHash(NamedTuple):
    contig_id: int
    position: int
    ref: SequenceKey

    @property
    def sql_key(self) -> str:
        """ Matches the Concat annotation in _get_loci_ids """
        return f"{self.position}_{self.ref}"

    @property
    def long_sequences(self) -> set[str]:
        return {s for s in (self.ref,) if isinstance(s, str)}

    def resolve(self, sequence_pk_by_seq: dict[str, Optional[int]]) -> Optional["LocusHash"]:
        """ Long sequence swapped for its pk. None = not in the DB, so no Locus can use it """
        ref = _resolve_sequence_key(self.ref, sequence_pk_by_seq)
        return None if ref is None else self._replace(ref=ref)


class VariantHash(NamedTuple):
    contig_id: int
    position: int
    ref: SequenceKey
    alt: SequenceKey
    svlen: int | str  # '' for null - 0 is a real value gene-level variants rely on, @see snpdb.gene_level_variants

    @property
    def locus_hash(self) -> LocusHash:
        return LocusHash(self.contig_id, self.position, self.ref)

    @property
    def sql_key(self) -> str:
        """ Matches the Concat annotation in _get_variant_ids """
        return f"{self.position}_{self.ref}_{self.alt}_{self.svlen}"

    @property
    def long_sequences(self) -> set[str]:
        return {s for s in (self.ref, self.alt) if isinstance(s, str)}

    def resolve(self, sequence_pk_by_seq: dict[str, Optional[int]]) -> Optional["VariantHash"]:
        """ Long sequences swapped for their pks. None = one isn't in the DB, so no Variant can use it """
        ref = _resolve_sequence_key(self.ref, sequence_pk_by_seq)
        alt = _resolve_sequence_key(self.alt, sequence_pk_by_seq)
        return None if ref is None or alt is None else self._replace(ref=ref, alt=alt)


AnyHash: TypeAlias = VariantHash | LocusHash

SHORT_SEQUENCE_MAX_LENGTH = 5


@cache
def _preloaded_sequences() -> frozenset[str]:
    """ Gene-level alts (<FUSION:HGNC:76>) are unique to their event, so are left to the per-batch query """
    sequences = {Variant.REFERENCE_ALT, *settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES}
    for length in range(1, SHORT_SEQUENCE_MAX_LENGTH + 1):
        sequences.update("".join(p) for p in itertools.product("GATCN", repeat=length))
    return frozenset(sequences)


@cache
def _preloaded_sequence_hashes() -> list[str]:
    return [sha256sum_str(seq) for seq in _preloaded_sequences()]


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
        preloaded_qs = Sequence.objects.filter(seq_sha256_hash__in=_preloaded_sequence_hashes())
        self.sequence_pk_by_seq: dict[str, int] = dict(preloaded_qs.values_list("seq", "pk"))
        self.long_sequence_pk_by_seq: dict[str, Optional[int]] = {}  # None = not in DB

    @cached_property
    def reference_seq_id(self) -> int:
        if reference_seq_id := self.sequence_pk_by_seq.get(Variant.REFERENCE_ALT):
            return reference_seq_id
        defaults = {
            "seq_sha256_hash": sha256sum_str(Variant.REFERENCE_ALT),
        }
        return Sequence.objects.get_or_create(seq=Variant.REFERENCE_ALT, defaults=defaults)[0].pk

    @staticmethod
    def _is_preloaded(seq: str) -> bool:
        return seq in _preloaded_sequences()

    def _get_sequence_key(self, seq: str) -> Optional[SequenceKey]:
        """ None = preloaded sequence not in DB, so no variant can use it """
        if self._is_preloaded(seq):
            return self.sequence_pk_by_seq.get(seq)
        return seq

    def _load_long_sequences(self, sequences: Iterable[str]):
        """ 1 query for every long sequence not looked up yet """
        new_sequence_by_hash = {sha256sum_str(seq): seq for seq in sequences
                                if seq not in self.long_sequence_pk_by_seq}
        if new_sequence_by_hash:
            self.long_sequence_pk_by_seq.update(dict.fromkeys(new_sequence_by_hash.values()))
            qs = Sequence.objects.filter(seq_sha256_hash__in=new_sequence_by_hash)
            for seq_sha256_hash, pk in qs.values_list("seq_sha256_hash", "pk"):
                self.long_sequence_pk_by_seq[new_sequence_by_hash[seq_sha256_hash]] = pk

    def _resolve_hashes(self, hashes: Iterable[Optional[AnyHash]]) -> list[Optional[AnyHash]]:
        """ Swaps long sequences for their pks. None where a hash is None or uses a sequence not in the DB """
        hashes = list(hashes)
        self._load_long_sequences(set().union(*(h.long_sequences for h in hashes if h is not None)))
        return [None if h is None else h.resolve(self.long_sequence_pk_by_seq) for h in hashes]

    def _get_variant_hash(self, contig_id: int, position: int, ref_key: SequenceKey, alt_key: SequenceKey,
                          svlen: Optional[int]) -> VariantHash:
        return VariantHash(contig_id, position, ref_key, alt_key, '' if svlen is None else svlen)

    def _get_ids_for_hashes(self, hashes: Iterable[Optional[AnyHash]],
                            get_queryset: Callable[[int, Collection[int]], QuerySet]) -> list[Optional[int]]:
        resolved_hashes = self._resolve_hashes(hashes)
        contig_positions = defaultdict(set)
        for h in resolved_hashes:
            if h is not None:
                contig_positions[h.contig_id].add(h.position)

        contig_hashes: dict[int, dict] = defaultdict(dict)
        for contig_id, positions in contig_positions.items():
            qs = get_queryset(contig_id, positions)
            contig_hashes[contig_id] = dict(qs.values_list("hash", "id"))

        return [None if h is None else contig_hashes[h.contig_id].get(h.sql_key) for h in resolved_hashes]

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

    def _get_loci_ids(self, loci_hashes: Iterable[LocusHash]) -> list[int]:
        annotate_kwargs = {
            "hash": Concat("position", Value("_"), "ref_id",
                           output_field=TextField())
        }

        def get_queryset(contig_id, positions):
            qs = Locus.objects.filter(contig_id=contig_id, position__in=positions)
            return qs.annotate(**annotate_kwargs)

        return self._get_ids_for_hashes(loci_hashes, get_queryset)

    def get_loci_ids(self, loci_hashes: Iterable[LocusHash], validate_not_null=True) -> list[int]:
        loci_ids = self._get_loci_ids(loci_hashes)
        if validate_not_null:
            if not all(loci_ids):  # Quick test if ok
                # Slowly find the first bad record
                for loci_hash, pk in zip(loci_hashes, loci_ids):
                    if not pk:
                        raise ValueError(f"Loci hash {loci_hash} had no PK in DB!")
        return loci_ids

    def filter_non_reference(self, variant_hashes: Iterable[VariantHash], variant_ids: Iterable[int]) -> list[int]:
        return [variant_id for variant_hash, variant_id in zip(variant_hashes, variant_ids)
                if variant_hash.alt != self.reference_seq_id]

    def get_variant_coordinate_hash(self, variant_coordinate: VariantCoordinate) -> Optional[VariantHash]:
        """ For VCF records (needs GenomeBuild supplied). None = uses a sequence not in the DB, so the variant
            can't be either. A hash isn't proof the variant exists - get_variant_ids returns None for those """

        variant_coordinate = variant_coordinate.as_internal_canonical_form(self.genome_build)
        if self.chrom_contig_id_mappings is None:
            raise ValueError("Need to initialise w/GenomeBuild to call get_variant_coordinate_hash")
        contig_id = self.chrom_contig_id_mappings[variant_coordinate.chrom]
        ref_key = self._get_sequence_key(variant_coordinate.ref)
        alt_key = self._get_sequence_key(variant_coordinate.alt)
        if ref_key is None or alt_key is None:
            return None
        return self._get_variant_hash(contig_id, variant_coordinate.position, ref_key, alt_key, variant_coordinate.svlen)

    def add(self, variant_coordinate: VariantCoordinate) -> Optional[VariantHash]:
        variant_coordinate = variant_coordinate.as_internal_canonical_form(self.genome_build)

        # If a preloaded sequence isn't known, variant is definitely unknown. Others are checked in batch_check
        if variant_hash := self.get_variant_coordinate_hash(variant_coordinate):
            # Maybe unknown, need to check
            self.variant_hashes.append(variant_hash)
            self.variant_coordinates.append(variant_coordinate)
        else:
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
        """ Inserts unknown from coordinates, updates sequence_pk_by_seq / long_sequence_pk_by_seq """
        sequences = set()
        for variant_coordinate in self.unknown_variant_coordinates:
            sequences.update((variant_coordinate.ref, variant_coordinate.alt))
        self._load_long_sequences(seq for seq in sequences if not self._is_preloaded(seq))

        unknown_sequence_by_hash = {}
        for seq in sequences:
            if self._is_preloaded(seq):
                known = seq in self.sequence_pk_by_seq
            else:
                known = self.long_sequence_pk_by_seq[seq] is not None
            if not known:
                unknown_sequence_by_hash[sha256sum_str(seq)] = seq

        if unknown_sequence_by_hash:
            # bulk create doesn't return pks with ignore_conflicts, and another import may have inserted some
            sequences = [Sequence(seq=seq, seq_sha256_hash=seq_sha256_hash)
                         for seq_sha256_hash, seq in unknown_sequence_by_hash.items()]
            Sequence.objects.bulk_create(sequences, ignore_conflicts=True)
            qs = Sequence.objects.filter(seq_sha256_hash__in=unknown_sequence_by_hash)
            for seq_sha256_hash, pk in qs.values_list("seq_sha256_hash", "pk"):
                seq = unknown_sequence_by_hash[seq_sha256_hash]
                if self._is_preloaded(seq):
                    self.sequence_pk_by_seq[seq] = pk
                else:
                    self.long_sequence_pk_by_seq[seq] = pk

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
        self._insert_unknown_sequences()  # All ref/alt need to be in the DB to resolve hashes to sequence pks

        # We need these to be ordered, so will use dict keys
        loci_hashes = {}  # Uses None as values, just want orderable
        variant_coordinates_by_hash = {}

        for variant_coordinate in self.unknown_variant_coordinates:
            variant_hash = self.get_variant_coordinate_hash(variant_coordinate)
            variant_coordinates_by_hash[variant_hash] = variant_coordinate
            loci_hashes[variant_hash.locus_hash] = None

        self.unknown_variant_coordinates.clear()
        loci_ids = self.get_loci_ids(loci_hashes, validate_not_null=False)  # Could have some unknowns as None
        unknown_loci_parts = []
        for resolved_locus_hash, locus_pk in zip(self._resolve_hashes(loci_hashes), loci_ids):
            if locus_pk is None:
                unknown_loci_parts.append(resolved_locus_hash)

        if unknown_loci_parts:
            self._insert_new_loci(unknown_loci_parts)

        loci_ids = self.get_loci_ids(loci_hashes)  # Validates all are not null
        locus_pk_by_hash = dict(zip(loci_hashes, loci_ids))
        variant_ids = self.get_variant_ids(variant_coordinates_by_hash, validate_not_null=False)
        unknown_variants = []
        resolved_variant_hashes = self._resolve_hashes(variant_coordinates_by_hash)
        for variant_hash, resolved_variant_hash, variant_pk in zip(variant_coordinates_by_hash, resolved_variant_hashes,
                                                                   variant_ids):
            if variant_pk is None:
                locus_pk = locus_pk_by_hash[variant_hash.locus_hash]
                vc = variant_coordinates_by_hash[variant_hash]
                unknown_variants.append((locus_pk, resolved_variant_hash.alt, vc.end, resolved_variant_hash.svlen))
        logging.debug("loci_hashes = %d, unknown_variants = %d", len(loci_hashes), len(unknown_variants))

        if unknown_variants:
            self._insert_new_variants(unknown_variants)

        self._file_batch_id += 1
