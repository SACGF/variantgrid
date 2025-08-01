import logging
import re
from functools import cached_property
from typing import Optional, Iterable, Union, Any

import django
import pydantic
from bioutils.sequences import reverse_complement
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models, IntegrityError
from django.db.models import Value, QuerySet
from django.db.models.deletion import CASCADE, DO_NOTHING
from django.db.models.fields import TextField
from django.db.models.functions.text import Concat
from django.db.models.query_utils import Q, FilteredRelation
from django.dispatch import receiver
from django.urls.base import reverse
from django.utils.text import slugify
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager
from pydantic import field_validator

from flags.models import FlagCollection, flag_collection_extra_info_signal, FlagInfos
from flags.models.models import FlagsMixin, FlagTypeContext
from library.django_utils.django_object_managers import ObjectManagerCachingRequest
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.genomics import format_chrom
from library.genomics.vcf_enums import VCFSymbolicAllele, INFO_LIFTOVER_SWAPPED_REF_ALT
from library.guardian_utils import admin_bot
from library.preview_request import PreviewModelMixin, PreviewKeyValue
from library.utils import FormerTuple, sha256sum_str
from snpdb.models import Wiki
from snpdb.models.models_clingen_allele import ClinGenAllele
from snpdb.models.models_enums import AlleleConversionTool, AlleleOrigin, ProcessingStatus
from snpdb.models.models_genome import Contig, GenomeBuild, GenomeBuildContig

LOCUS_PATTERN = re.compile(r"^([^:]+)\s*:\s*(\d+)[,\s]*([GATC]+)$", re.IGNORECASE)
LOCUS_NO_REF_PATTERN = re.compile(r"^([^:]+)\s*:\s*(\d+)$")
VARIANT_PATTERN = re.compile(r"^(MT|(?:chr)?(?:[XYM]|\d+))\s*:\s*(\d+)[,\s]*([GATC]+)>(=|[GATC]+)$", re.IGNORECASE)
# This is our internal format for symbolic (ie <DEL>/<DUP> etc)
VARIANT_SYMBOLIC_PATTERN = re.compile(r"^(MT|(?:chr)?(?:[XYM]|\d+))\s*:\s*(\d+)\s*-\s*(\d+)\s*<(DEL|DUP|INS|INV|CNV)>$", re.IGNORECASE)
# matches anything hgvs-like before any fixes
HGVS_UNCLEANED_PATTERN = re.compile(r"(^(N[MC]_|ENST)\d+.*:|[cnmg]\.|[^:]:[cnmg]).*\d+", re.IGNORECASE)


class Allele(FlagsMixin, PreviewModelMixin, models.Model):
    """ Genome build independent - ie GRCh37 and GRCh38 variants for same change point to same allele
        This is generally done via ClinGen Allele Registry, but sometimes that can fail.
        Linked against Variant with VariantAllele below """

    clingen_allele = models.OneToOneField(ClinGenAllele, null=True, on_delete=CASCADE)

    objects = ObjectManagerCachingRequest()

    class Meta:
        base_manager_name = 'objects'

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-a p-1 text-light border rounded bg-dark"

    @property
    def preview(self) -> 'PreviewData':
        return self.preview_with(
            identifier=f"Allele ({self.pk})",
            summary_extra=[PreviewKeyValue("ClinGenAllele ID", str(self.clingen_allele) if self.clingen_allele else "Unknown")]
        )

    def get_absolute_url(self):
        # will show allele if there is one, otherwise go to variant page
        return reverse('view_allele', kwargs={"allele_id": self.id})

    def flag_type_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk="allele")

    @property
    def compact_str(self):
        if clingen := self.clingen_allele:
            return str(clingen)
        else:
            return str(self)

    @property
    def metrics_logging_key(self):
        return "allele_id", self.pk

    @cached_property
    def clingen_error(self):
        error = None
        if va := self.variantallele_set.filter(clingen_error__isnull=False).first():
            error = va.clingen_error
        return error

    def variant_alleles(self) -> QuerySet['VariantAllele']:
        return self.variantallele_set.select_related('variant__locus', 'variant__locus__contig', 'variant__locus__ref', 'variant__alt').order_by("genome_build__name")

    @cached_property
    def grch37(self) -> Optional['Variant']:
        return self.variant_for_build_optional(genome_build=GenomeBuild.grch37())

    @cached_property
    def grch38(self) -> Optional['Variant']:
        return self.variant_for_build_optional(genome_build=GenomeBuild.grch38())

    @cached_property
    def variants(self):
        return Variant.objects.filter(pk__in=self.variant_alleles().values_list('variant', flat=True))

    def variant_for_any_build(self, preferred_genome_build: GenomeBuild, best_attempt=True) -> 'Variant':
        vas = self.variant_alleles()
        va = None
        if preferred_genome_build:
            va = vas.filter(genome_build=preferred_genome_build).first()
            if not va and not best_attempt:
                raise ValueError(f'Could not find a variant in allele {self.id} for build {preferred_genome_build}')
        if not va:
            va = vas.first()
        if va:
            return va.variant
        raise ValueError(f'Could not find any variants in allele {self.id}')

    def variant_for_build_optional(self, genome_build: GenomeBuild) -> Optional['Variant']:
        if va := self.variant_alleles().filter(genome_build=genome_build).first():
            return va.variant

    def variant_for_build(self, genome_build: GenomeBuild) -> 'Variant':
        v = self.variant_for_build_optional(genome_build)
        if v is None:
            raise ValueError(f'Could not find a variant in allele {self.id} for build {genome_build}')
        return v

    def merge(self, allele_linking_tool, other_allele: "Allele") -> bool:
        """ Merge other_allele into this allele """

        if self == other_allele:
            raise ValueError(f"Attempt to merge {self} to itself!")

        can_merge = True
        merge_log_message = f"{other_allele} merge into {self}"
        other_clingen_allele = other_allele.clingen_allele
        if other_clingen_allele and self.clingen_allele:
            can_merge = False
            merge_log_message = f"Error performing {merge_log_message}: both have ClinGen Alleles!"

        AlleleMergeLog.objects.create(old_allele=other_allele,
                                      new_allele=self,
                                      allele_linking_tool=allele_linking_tool,
                                      success=can_merge,
                                      message=merge_log_message)

        if can_merge:
            if other_clingen_allele:
                # Move across ClinGen Allele (may not have been possible to retrieve in all builds, but at least one
                # links there, and can't have another, so it'll work)
                other_allele.clingen_allele = None
                other_allele.save()

                self.clingen_allele = other_clingen_allele
                self.save()

            if other_fc := other_allele.flag_collection:
                other_fc.flag_set.update(collection=self.flag_collection_safe)
                existing_fc_cc_names = self.flag_collection.clinicalcontext_set.values_list("name", flat=True)
                other_fc.clinicalcontext_set.exclude(name__in=existing_fc_cc_names).update(flag_collection=self.flag_collection)
                other_fc.classification_set.update(flag_collection=self.flag_collection)
            existing_allele_cc_names = self.clinicalcontext_set.values_list("name", flat=True)
            other_allele.clinicalcontext_set.exclude(name__in=existing_allele_cc_names).update(allele=self)
            for va in other_allele.variantallele_set.all():
                try:
                    va.allele = self
                    va.clingen_error = None  # clear any errors
                    va.allele_linking_tool = allele_linking_tool
                    va.save()
                except IntegrityError:
                    logging.warning("VariantAllele exists with allele/build/variant of %s/%s/%s - deleting this one",
                                    va.allele, va.genome_build, va.variant)
                    va.delete()

        return can_merge

    @property
    def build_names(self) -> str:
        return ", ".join(sorted(self.variantallele_set.values_list("genome_build__name", flat=True)))

    @staticmethod
    def missing_variants_for_build(genome_build) -> QuerySet['Allele']:
        alleles_with_variants_qs = Allele.objects.filter(variantallele__isnull=False)
        return alleles_with_variants_qs.filter(~Q(variantallele__genome_build=genome_build))

    def __str__(self):
        name = f"Allele {self.pk}"
        if self.clingen_allele:
            name += f" ({self.clingen_allele})"
        return name

    def __format__(self, format_spec: str):
        if format_spec == 'CA' and (cligen_allele := self.clingen_allele):
            return str(cligen_allele)
        else:
            return f"Allele {self.pk}"


@receiver(flag_collection_extra_info_signal, sender=FlagCollection)
def get_extra_info(flag_infos: FlagInfos, user: User, **kwargs):  # pylint: disable=unused-argument
    alleles = Allele.objects.filter(flag_collection__in=flag_infos.ids)
    allele: Allele
    for allele in alleles:
        flag_infos.set_extra_info(allele.flag_collection_id, {
            'label': f'Allele {allele.id}'
        }, source_object=allele)


class AlleleMergeLog(TimeStampedModel):
    """ Keep track of calls to Allele.merge() """
    old_allele = models.ForeignKey(Allele, related_name="old_allele_merge", on_delete=CASCADE)
    new_allele = models.ForeignKey(Allele, related_name="new_allele_merge", on_delete=CASCADE)
    allele_linking_tool = models.CharField(max_length=2, choices=AlleleConversionTool.choices)
    success = models.BooleanField(default=True)
    message = models.TextField(null=True)


class VariantCoordinate(FormerTuple, pydantic.BaseModel):
    """ This stores coordinates, when you want to use it, be sure to call either:
            * as_external_explicit() - if you want to interface with eg VCF
            * as_internal_symbolic() - to store in internal Database
    """
    chrom: str
    position: int
    ref: str
    alt: str
    svlen: Optional[int] = None

    @field_validator('svlen')
    @classmethod
    def validate_svlen(cls, value):
        if value is not None and settings.VARIANT_SYMBOLIC_ALT_ENABLED is False:
            raise ValueError('Symbolic variants disabled via settings.')
        return value

    @property
    def as_tuple(self) -> tuple:
        return self.chrom, self.position, self.ref, self.alt, self.svlen

    @property
    def end(self) -> int:
        """
            This corresponds to VCF INFO["END"] which is defined in spec as:

            END position of the longest variant described in this record. The END of each allele is defined as:
                * Non-symbolic alleles: POS + length of REF allele − 1.
                * <INS> symbolic structural variant alleles: POS + length of REF allele − 1.
                * <DEL>, <DUP>, <INV>, and <CNV> symbolic structural variant alleles:, POS + SVLEN
        """
        if self.is_symbolic():
            # We don't support <INV> so don't need to worry about it
            return self.position + abs(self.svlen)
        return self.position + len(self.ref) - 1

    def __lt__(self, other):
        return self.as_tuple < other.as_tuple

    def __str__(self):
        return self.format()

    def __hash__(self):
        return self.as_tuple.__hash__()

    def format(self):
        if Sequence.allele_is_symbolic(self.alt):
            range_end = self.position + abs(self.svlen)
            return f"{self.chrom}:{self.position}-{range_end} {self.alt}"
        else:
            return f"{self.chrom}:{self.position} {self.ref}>{self.alt}"

    @staticmethod
    def from_variant_match(match, genome_build: Optional[GenomeBuild] = None):
        chrom = match.group(1)
        if genome_build:
            chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        position = int(match.group(2))
        ref = match.group(3)
        alt = match.group(4)
        vc = VariantCoordinate.from_explicit_no_svlen(chrom, position, ref, alt)
        return vc.as_internal_canonical_form(genome_build)

    @staticmethod
    def from_symbolic_match(match, genome_build):
        chrom, start, range_end, alt = match.groups()
        chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        alt = alt.upper()
        start = int(start)
        range_end = int(range_end)
        if range_end < start + 1:
            raise ValueError("End coordinate must be at least 1 base after start")

        if alt == "DEL":
            svlen = start - range_end
        else:
            svlen = range_end - start
        alt = "<" + alt + ">"  # captured what was inside of brackets ie <()>
        contig_sequence = genome_build.genome_fasta.fasta[chrom]
        # 0-based
        ref = contig_sequence[start-1:start].upper()
        vc = VariantCoordinate(chrom=chrom, position=start, ref=ref, alt=alt, svlen=svlen)
        return vc.as_internal_canonical_form(genome_build)

    @staticmethod
    def from_string(variant_string: str, genome_build):
        """ Pass in genome build to be able to set REF from symbolic (will be N otherwise) """
        if full_match := VARIANT_PATTERN.fullmatch(variant_string):
            return VariantCoordinate.from_variant_match(full_match, genome_build)
        elif full_match := VARIANT_SYMBOLIC_PATTERN.fullmatch(variant_string):
            return VariantCoordinate.from_symbolic_match(full_match, genome_build)
        regex_patterns = ", ".join((str(s) for s in (VARIANT_PATTERN, VARIANT_SYMBOLIC_PATTERN)))
        raise ValueError(f"{variant_string=} did not match against {regex_patterns=}")

    @staticmethod
    def from_explicit_no_svlen(chrom: str, position: int, ref: str, alt: str) -> 'VariantCoordinate':
        """ Initialise w/o providing an end """

        ref = ref.strip().upper()
        alt = alt.strip().upper()
        if Sequence.allele_is_symbolic(alt):
            raise ValueError("Must pass 'svlen' when using symbolic alt")
        return VariantCoordinate(chrom=chrom, position=position, ref=ref, alt=alt)

    def is_symbolic(self):
        return Sequence.allele_is_symbolic(self.alt)

    def calculated_reference(self, genome_build) -> str:
        contig_sequence = genome_build.genome_fasta.fasta[self.chrom]
        # reference sequence is 0-based
        ref_sequence = contig_sequence[self.position - 1:self.end]
        return ref_sequence.upper()

    def as_external_explicit(self, genome_build) -> 'VariantCoordinate':
        """ explicit ref/alt """
        if self.is_symbolic():
            if self.svlen is None:
                raise ValueError(f"{self} has 'svlen' = None")

            ref_sequence = self.calculated_reference(genome_build)
            if self.alt == VCFSymbolicAllele.DEL:
                ref = ref_sequence
                alt = ref_sequence[0]
            elif self.alt == VCFSymbolicAllele.DUP:
                ref = ref_sequence[0]
                alt = ref_sequence
            elif self.alt == VCFSymbolicAllele.INV:
                ref = ref_sequence
                alt = reverse_complement(ref_sequence)
            else:
                raise ValueError(f"Unknown symbolic alt of '{self.alt}'")

            vc = VariantCoordinate(chrom=self.chrom, position=self.position, ref=ref, alt=alt)
            # print(f"Symbolic = {self} -> {vc}")
            return vc

        if self.alt == Variant.REFERENCE_ALT:
            # Convert from our internal format (alt='=' for ref) to explicit
            alt = self.ref
        else:
            alt = self.alt
        return VariantCoordinate(chrom=self.chrom, position=self.position, ref=self.ref, alt=alt)

    @property
    def max_sequence_length(self) -> int:
        if self.is_symbolic():
            return abs(self.svlen)
        else:
            return max(len(self.ref), len(self.alt))

    def as_internal_symbolic(self, genome_build: GenomeBuild,
                             min_symbolic_alt_size=settings.VARIANT_SYMBOLIC_ALT_SIZE) -> 'VariantCoordinate':
        """ Internal format - alt can be <DEL> or <DUP>
            Uses our internal reference representation
        """
        if self.is_symbolic():
            return self

        ref = self.ref
        alt = self.alt  # default, only change if symbolic
        svlen = None
        if self.alt in (Variant.REFERENCE_ALT, self.ref):
            alt = Variant.REFERENCE_ALT
        else:
            # Could be symbolic
            ref_length = len(ref)
            alt_length = len(self.alt)
            diff = alt_length - ref_length
            if ref_length == 1:
                if alt_length >= min_symbolic_alt_size:
                    # Possible dup
                    # TODO: Can probably remove HGVS dependency from here - just look directly at sequence
                    from genes.hgvs import HGVSMatcher
                    matcher = HGVSMatcher(genome_build)
                    hgvs_variant = matcher.variant_coordinate_to_hgvs_variant(self)
                    if hgvs_variant.mutation_type == 'dup':
                        ref = self.ref[0]
                        alt = VCFSymbolicAllele.DUP
                        svlen = diff
            elif alt_length == 1 and self.alt == self.ref[0]:
                if ref_length >= min_symbolic_alt_size:
                    ref = self.ref[0]
                    alt = VCFSymbolicAllele.DEL
                    svlen = diff

            elif ref_length == alt_length:
                if ref_length > min_symbolic_alt_size:
                    if self.ref == reverse_complement(self.alt):
                        ref = self.ref[0]
                        alt = VCFSymbolicAllele.INV
                        svlen = len(self.ref) - 1  # explicit inv had same length ref/alt, now we have len(ref) == 1
        return VariantCoordinate(chrom=self.chrom, position=self.position, ref=ref, alt=alt, svlen=svlen)

    def as_internal_canonical_form(self, genome_build: GenomeBuild) -> 'VariantCoordinate':
        """ Make sure we only have 1 representation for a variant """

        # Easiest way is to just convert to symbolic then check svlen
        vc_symbolic = self.as_internal_symbolic(genome_build)
        if vc_symbolic.svlen and abs(vc_symbolic.svlen) >= settings.VARIANT_SYMBOLIC_ALT_SIZE:
            vc = vc_symbolic
        elif self.is_symbolic():
            vc = self.as_external_explicit(genome_build)
        else:
            vc = self

        if vc.alt in (Variant.REFERENCE_ALT, vc.ref):
            vc.alt = Variant.REFERENCE_ALT
        return vc

    def as_contig_accession(self, genome_build: GenomeBuild) -> 'VariantCoordinate':
        contig = genome_build.chrom_contig_mappings[self.chrom]
        return VariantCoordinate(chrom=contig.refseq_accession,
                                 position=self.position, ref=self.ref, alt=self.alt, svlen=self.svlen)


class Sequence(models.Model):
    """
        We want to guarantee seq is unique (so Locus/Variant can have unique constraints)
        Postgres by default uses indexes for constraints, and large TextFields give error of:
        "index row requires x bytes, maximum size is 8191"

        The easiest solution is to md5sum seq and make the constraint on that. Another possible solution is to use
        Gist indexes but that requires installing the btree_gist extension (requires postgres Admin rights).
        Django 3 has ExclusionConstraint, Postgres contrib has BtreeGistExtension to add via migration

        Note: Even after the introduction of symbolic alts, there are still some really long sequences in here (~10kb)
        due to large substitutions which we don't represent symbolically
    """
    seq = models.TextField()
    seq_sha256_hash = models.TextField(null=True, unique=True)

    def save(self, force_insert=False, force_update=False, using=None, update_fields=None):
        if not self.seq_sha256_hash:
            self.seq_sha256_hash = sha256sum_str(self.seq)
        super().save(force_insert=force_insert, force_update=force_update, using=using, update_fields=update_fields)

    @staticmethod
    def abbreviate(s: str, max_length: int = 20):
        if len(s) > max_length:
            s = f"{s[:3]}...{s[-3:]}"
        return s

    def __len__(self) -> int:
        return len(self.seq)

    def __str__(self):
        return self.abbreviate(self.seq)

    @staticmethod
    def get_pk_by_seq(q=None):
        qs = Sequence.objects.all()
        if q:
            qs = qs.filter(q)
        return dict(qs.values_list("seq", "pk"))

    @staticmethod
    def allele_is_symbolic(seq: Union[str, 'Sequence']) -> bool:
        return seq.startswith("<") and seq.endswith(">")

    def is_standard_sequence(self) -> bool:
        """ only contains G/A/T/C/N """
        return not re.match(r"[^GATCN]", self.seq)

    def is_symbolic(self) -> bool:
        return self.allele_is_symbolic(self.seq)

    def startswith(self, prefix: str) -> bool:
        """ To match str method, so allele_is_symbolic works above """
        return self.seq.startswith(prefix)


class Locus(models.Model):
    """ 1 per line in a VCF file (multiple Variants with different alt alleles point to the same locus)
        There is only 1 Locus for a given chrom/position/ref per database (handled via insertion queues) """

    contig = models.ForeignKey(Contig, on_delete=CASCADE)
    position = models.IntegerField(db_index=True)
    ref = models.ForeignKey(Sequence, on_delete=CASCADE)

    class Meta:
        unique_together = ("contig", "position", "ref")

    @property
    def chrom(self):
        return self.contig.name

    def __str__(self):
        return f"{self.chrom}:{self.position} {self.ref}"


class Variant(PreviewModelMixin, models.Model):
    """ Variants represent the different alleles at a locus
        Usually 2+ per line in a VCF file (ref + >= 1 alts pointing to the same locus for the row)
        There is only 1 Variant for a given locus/alt per database (handled via insertion queues) """

    REFERENCE_ALT = "="
    _BASES = "GATC"
    _REGEX_2_PLUS_BASES = "^[GATC]{2,}$"

    locus = models.ForeignKey(Locus, on_delete=CASCADE)
    alt = models.ForeignKey(Sequence, on_delete=CASCADE)
    # end is a calculated field, but we store as an optimisation for overlap queries (start = locus.position)
    end = models.IntegerField()
    svlen = models.IntegerField(null=True, blank=True)  # For symbolic variants, difference in length between REF and ALT alleles

    class Meta:
        unique_together = ("locus", "alt", "svlen")  # Possible to have eg CNV with same alt = <INS> but diff end

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-v p-1 text-light border rounded bg-dark"

    @staticmethod
    def get_chrom_q(chrom):
        return Q(locus__contig__name__iexact=chrom) | Q(locus__contig__ucsc_name__iexact=chrom)

    @staticmethod
    def get_contigs_q(genome_build: GenomeBuild) -> Q:
        """ Restrict to contigs in a genome build """
        return Q(locus__contig__genomebuildcontig__genome_build=genome_build)

    @staticmethod
    def get_reference_q() -> Q:
        return Q(alt__seq=Variant.REFERENCE_ALT)

    @staticmethod
    def get_no_reference_q() -> Q:
        return ~Variant.get_reference_q()

    @staticmethod
    def get_snp_q() -> Q:
        return Q(locus__ref__seq__in=Variant._BASES) & Q(alt__seq__in=Variant._BASES)

    @staticmethod
    def get_insertion_q() -> Q:
        return Q(locus__ref__seq__in=Variant._BASES) & Q(alt__seq__regex=Variant._REGEX_2_PLUS_BASES)

    @staticmethod
    def get_deletion_q() -> Q:
        return Q(locus__ref__seq__regex=Variant._REGEX_2_PLUS_BASES) & Q(alt__seq__in=Variant._BASES)

    @staticmethod
    def get_indel_q() -> Q:
        return Variant.get_insertion_q() | Variant.get_deletion_q()

    @staticmethod
    def get_complex_subsitution_q() -> Q:
        return Q(locus__ref__seq__regex=Variant._REGEX_2_PLUS_BASES) & Q(alt__seq__regex=Variant._REGEX_2_PLUS_BASES)

    @staticmethod
    def get_symbolic_q() -> Q:
        return Q(svlen__isnull=False)

    @staticmethod
    def annotate_variant_string(qs, name="variant_string", path_to_variant=""):
        """ Return a "1:123321 G>C" style string in a query """
        kwargs = {name: Concat(f"{path_to_variant}locus__contig__name", Value(":"),
                               f"{path_to_variant}locus__position", Value(" "),
                               f"{path_to_variant}locus__ref__seq", Value(">"),
                               f"{path_to_variant}alt__seq", output_field=TextField())}
        return qs.annotate(**kwargs)

    @staticmethod
    def validate(genome_build, chrom, position) -> list[str]:
        errors = []
        try:
            contig = genome_build.chrom_contig_mappings[chrom]
            position = int(position)
            if not (0 < position < contig.length):
                errors.append(f'position "{position}" is outside contig "{contig}" length={contig.length}')
        except KeyError:
            errors.append(f"Chromsome/contig '{chrom}' not a valid in genome build {genome_build}")
        return errors

    @staticmethod
    def format_tuple(chrom, position, ref, alt, svlen, abbreviate=False) -> str:
        if abbreviate and not Sequence.allele_is_symbolic(alt):
            ref = Sequence.abbreviate(ref)
            alt = Sequence.abbreviate(alt)
        vc = VariantCoordinate(chrom=chrom, position=position, ref=ref, alt=alt, svlen=svlen)
        return vc.format()

    @staticmethod
    def get_from_string(variant_string: str, genome_build: GenomeBuild) -> Optional['Variant']:
        variant_coordinate = VariantCoordinate.from_string(variant_string, genome_build)
        try:
            return Variant.get_from_variant_coordinate(variant_coordinate, genome_build)
        except Variant.DoesNotExist:
            return None

    @staticmethod
    def qs_from_variant_coordinate(variant_coordinate: VariantCoordinate, genome_build: GenomeBuild) -> QuerySet['Variant']:
        variant_coordinate = variant_coordinate.as_internal_symbolic(genome_build)
        params = ["locus__contig__name", "locus__position", "locus__ref__seq", "alt__seq", "svlen"]
        return Variant.objects.filter(locus__contig__genomebuildcontig__genome_build=genome_build,
                                      **dict(zip(params, variant_coordinate)))

    @staticmethod
    def get_from_variant_coordinate(variant_coordinate: VariantCoordinate, genome_build: GenomeBuild) -> 'Variant':
        return Variant.qs_from_variant_coordinate(variant_coordinate, genome_build).get()

    @cached_property
    def genome_builds(self) -> set['GenomeBuild']:
        gbc_qs = GenomeBuildContig.objects.filter(genome_build__in=GenomeBuild.builds_with_annotation(),
                                                  contig__locus__variant=self)
        return {gbc.genome_build for gbc in gbc_qs}

    @property
    def any_genome_build(self) -> GenomeBuild:
        return next(iter(self.genome_builds))

    @cached_property
    def coordinate(self) -> VariantCoordinate:
        locus = self.locus
        contig = locus.contig
        return VariantCoordinate(chrom=contig.name, position=locus.position,
                                 ref=locus.ref.seq, alt=self.alt.seq, svlen=self.svlen)

    @staticmethod
    def is_ref_alt_reference(ref, alt):
        return alt in (ref, '.')

    @property
    def is_reference(self) -> bool:
        return self.alt.seq == self.REFERENCE_ALT

    @cached_property
    def vcf_alt(self) -> str:
        """ Return the base as a string (not as our special REFERENCE_ALT) """
        if self.is_reference:
            return self.locus.ref.seq
        return self.alt.seq

    @property
    def is_standard_variant(self) -> bool:
        """ Variant alt sequence is standard [GATCN] (ie not special or reference) """
        # locus.ref should always be standard...
        return self.alt.is_standard_sequence()

    @property
    def is_indel(self) -> bool:
        return self.is_insertion or self.is_deletion

    @property
    def is_insertion(self) -> bool:
        if self.alt.is_symbolic:
            return self.alt.seq == "<INS>"
        return self.alt.seq != Variant.REFERENCE_ALT and len(self.locus.ref.seq) < len(self.alt.seq)

    @property
    def is_deletion(self) -> bool:
        if self.alt.is_symbolic:
            return self.alt.seq == VCFSymbolicAllele.DEL
        return self.alt.seq != Variant.REFERENCE_ALT and len(self.locus.ref) > len(self.alt.seq)

    @property
    def is_symbolic(self) -> bool:
        return self.locus.ref.is_symbolic() or self.alt.is_symbolic()

    @property
    def can_make_g_hgvs(self) -> bool:
        """ Can't form ones with some symbolic variants (eg <INS>) """
        if self.is_symbolic:
            return self.alt.seq in {VCFSymbolicAllele.DEL, VCFSymbolicAllele.DUP, VCFSymbolicAllele.INV}
        return True

    @property
    def _clingen_allele_size(self) -> int:
        if self.is_symbolic:
            allele_size = abs(self.svlen)
            # Bug in ClinGenAlleleRegistry - dup seems to calculate max size wrong leading to a
            # max size of 5k, @see https://github.com/BRL-BCM/Allele-Registry/issues/6
            if self.alt.seq == VCFSymbolicAllele.DUP:
                allele_size *= 2
        else:
            allele_size = len(self.locus.ref) + len(self.alt)
        return allele_size

    @property
    def can_have_clingen_allele(self) -> bool:
        return self._clingen_allele_size <= ClinGenAllele.CLINGEN_ALLELE_MAX_ALLELE_SIZE and self.can_make_g_hgvs

    @property
    def can_have_annotation(self) -> bool:
        return not self.is_reference

    @property
    def can_have_c_hgvs(self) -> bool:
        return self.can_have_annotation and self.svlen is None or abs(self.svlen) <= settings.HGVS_MAX_SEQUENCE_LENGTH

    def as_tuple(self) -> tuple[str, int, str, str, int]:
        return self.locus.contig.name, self.locus.position, self.locus.ref.seq, self.alt.seq, self.svlen

    def is_abbreviated(self):
        return str(self) != self.full_string

    @cached_property
    def full_string(self):
        """ No abbreviation """
        return self.format_tuple(*self.as_tuple())

    def __str__(self):
        return self.format_tuple(self.locus.contig.name, self.locus.position, self.locus.ref.seq, self.alt.seq, self.svlen)

    def get_absolute_url(self):
        # will show allele if there is one, otherwise go to variant page
        return reverse('view_allele_from_variant', kwargs={"variant_id": self.pk})

    @cached_property
    def allele(self) -> Optional[Allele]:
        va = VariantAllele.objects.filter(variant=self).first()
        if va:
            return va.allele
        return None

    @property
    def metrics_logging_key(self) -> tuple[str, Any]:
        if allele := self.allele:
            return "allele_id", allele.pk
        return "variant_id", self.pk

    @property
    def equivalent_variants(self) -> Iterable['Variant']:
        allele = self.allele
        if not allele:
            return [self]
        return Variant.objects.filter(variantallele__allele=allele)

    def get_canonical_transcript_annotation(self, genome_build) -> Optional['VariantTranscriptAnnotation']:
        from annotation.models import VariantAnnotationVersion
        vav = VariantAnnotationVersion.latest(genome_build)
        return self.varianttranscriptannotation_set.filter(version=vav, canonical=True).first()

    def get_best_variant_transcript_annotation(self, genome_build) -> Optional['VariantTranscriptAnnotation']:
        from annotation.models import VariantAnnotationVersion
        vav = VariantAnnotationVersion.latest(genome_build)
        if can := self.varianttranscriptannotation_set.filter(version=vav, canonical=True).first():
            return can
        if version := self.varianttranscriptannotation_set.filter(version=vav).first():
            return version
        if any_at_all := self.varianttranscriptannotation_set.first():
            return any_at_all

    def get_canonical_c_hgvs(self, genome_build):
        c_hgvs = None
        if cta := self.get_canonical_transcript_annotation(genome_build):
            c_hgvs = cta.get_hgvs_c_with_symbol()
        return c_hgvs

    @property
    def start(self):
        return self.locus.position

    @property
    def length(self) -> int:
        return self.end - self.locus.position

    @property
    def sort_string(self) -> str:
        padded_contig = self.locus.contig.name or ''
        if padded_contig.isnumeric():
            padded_contig = padded_contig.rjust(2, '0')
        else:
            padded_contig = padded_contig.rjust(2, '9')
        padded_position = f'{self.locus.position or 0:09d}'
        return f"{padded_contig}-{padded_position}"


class VariantWiki(Wiki):
    variant = models.OneToOneField(Variant, on_delete=CASCADE)


class VariantAllele(TimeStampedModel):
    """ It's possible for multiple variants from the same genome build to
        resolve to the same allele (due to our normalization not being the same as ClinGen
        or 2 loci in a genome build being represented by 1 loci in the build being used
        by ClinGen) - but it's not likely. It's a bug to have the same 3 variant/build/allele
        so we can add that unique_together constraint

        We only expect to store Alleles for a small fraction of Variants
        So don't want them on the Variant object - instead do 1-to-1 """

    # Some builds share contigs (eg GRCh37/38 share MT and some unplaced scaffolds) - in those cases
    # we'll have the same variant linked through different VariantAlleles (so it can't be 1-to-1)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    allele = models.ForeignKey(Allele, on_delete=CASCADE)
    origin = models.CharField(max_length=1, choices=AlleleOrigin.choices)
    allele_linking_tool = models.CharField(max_length=2, choices=AlleleConversionTool.choices)
    # We call CAR to populate ClinGenAllele on Allele - if there's any error have to store it here
    # as we need a PK back to create ClinGenAllele() model, and can fail on 1 build succeed on another
    clingen_error = models.JSONField(null=True)  # null on success

    class Meta:
        unique_together = ("variant", "genome_build", "allele")

    @property
    def canonical_c_hgvs(self):
        return self.variant.get_canonical_c_hgvs(self.genome_build)

    def needs_clingen_call(self):
        if settings.CLINGEN_ALLELE_REGISTRY_LOGIN and self.allele.clingen_allele is None:
            if self.clingen_error:
                # Retry if server was down
                return self.clingen_error.get("errorType") == ClinGenAllele.CLINGEN_ALLELE_SERVER_ERROR_TYPE
            return True
        return False

    def __str__(self):
        s = f"{self.allele} - {self.variant_id}({self.genome_build}"
        if linking_tool := self.get_allele_linking_tool_display():
            s += f"/{linking_tool=}"
        return s


class VariantCollection(RelatedModelsPartitionModel):
    """ A set of variants - usually used as a cached result """

    RECORDS_BASE_TABLE_NAMES = ["snpdb_variantcollectionrecord"]
    RECORDS_FK_FIELD_TO_THIS_MODEL = "variant_collection_id"
    PARTITION_LABEL_TEXT = "variant_collection"
    name = models.TextField(null=True)
    count = models.IntegerField(null=True)
    status = models.CharField(max_length=1, choices=ProcessingStatus.choices, default=ProcessingStatus.CREATED)

    @property
    def variant_collection_alias(self):
        return f"variantcollection_{self.pk}"

    def get_annotation_kwargs(self, **kwargs):
        vcr_condition = Q(variantcollectionrecord__variant_collection=self)
        return {self.variant_collection_alias: FilteredRelation('variantcollectionrecord', condition=vcr_condition)}

    def get_arg_q_dict(self) -> dict[Optional[str], set[Q]]:
        if self.status != ProcessingStatus.SUCCESS:
            raise ValueError(f"{self}: status {self.get_status_display()} != SUCCESS")

        q = Q(**{f"{self.variant_collection_alias}__isnull": False})
        return {self.variant_collection_alias: {str(q): q}}

    def __lt__(self, other):
        return self.pk < other.pk

    def __str__(self):
        return f"VariantCollection: {self.pk} ({self.name})"


class VariantCollectionRecord(models.Model):
    variant_collection = models.ForeignKey(VariantCollection, on_delete=DO_NOTHING)  # handled via drop partition
    variant = models.ForeignKey(Variant, on_delete=CASCADE)


class AlleleSource(models.Model):
    """ Provides a source of alleles for liftover pipelines. """
    objects = InheritanceManager()

    def get_genome_build(self):
        return None

    def get_variants_qs(self):
        return Variant.objects.none()

    def get_allele_qs(self):
        return Allele.objects.filter(variantallele__variant__in=self.get_variants_qs())

    def liftover_complete(self, genome_build: GenomeBuild):
        """ This is called at the end of a liftover pipeline (once per build) """


class VariantAlleleSource(AlleleSource):
    variant_allele = models.ForeignKey(VariantAllele, on_delete=CASCADE)

    def get_genome_build(self):
        return self.variant_allele.genome_build

    def get_variants_qs(self):
        return Variant.objects.filter(variantallele=self.variant_allele)


class VariantAlleleCollectionSource(AlleleSource):
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)

    def get_genome_build(self):
        return self.genome_build

    def get_variants_qs(self):
        return Variant.objects.filter(variantallele__in=self.get_variant_allele_ids())

    def get_variant_allele_ids(self):
        return self.variantallelecollectionrecord_set.values_list("variant_allele", flat=True)


class VariantAlleleCollectionRecord(models.Model):
    collection = models.ForeignKey(VariantAlleleCollectionSource, on_delete=CASCADE)
    variant_allele = models.ForeignKey(VariantAllele, on_delete=CASCADE)


liftover_run_complete_signal = django.dispatch.Signal()


class LiftoverRun(TimeStampedModel):
    """ Liftover pipeline involves reading through a VCF where ID is set to Allele.pk and then creating
        VariantAllele entries for the variant/allele

        Some AlleleConversionTools (eg ClinGen AlleleRegistry) we can write the VCF in the desired genome build
        For others (NCBI Remap) we need to write the source genome build VCF first

        Alleles must have already been created - allele_source used to retrieve them

        The VCF (in genome_build build) is set in UploadedFile for the UploadPipeline """
    user = models.ForeignKey(User, on_delete=CASCADE)
    allele_source = models.ForeignKey(AlleleSource, null=True, on_delete=CASCADE)
    conversion_tool = models.CharField(max_length=2, choices=AlleleConversionTool.choices)
    source_vcf = models.TextField(null=True)
    source_genome_build = models.ForeignKey(GenomeBuild, null=True, on_delete=CASCADE,
                                            related_name="liftover_source_genome_build")
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)  # destination

    @staticmethod
    def get_clingen_auto_fail_liftover_run(genome_build: GenomeBuild) -> 'LiftoverRun':
        """ We can know straight away we're going to fail a ClinGen run (missing ClinGenAllele, no data for the build)
            but would still like to record the failure - so put them all in 1 LiftoverRun per build """
        kwargs = {
            "user": admin_bot(),
            "genome_build": genome_build,
            "conversion_tool": AlleleConversionTool.CLINGEN_ALLELE_REGISTRY,
        }
        # There should have been 1 created in 'one_off_legacy_populate_allele_liftover'
        lr = LiftoverRun.objects.filter(**kwargs).order_by("pk").first()
        if not lr:
            # Otherwise create it
            lr = LiftoverRun.objects.create(**kwargs)
        return lr

    def get_absolute_url(self):
        return reverse("view_liftover_run", kwargs={"liftover_run_id": self.pk})

    def get_allele_source(self) -> AlleleSource:
        """ Returns subclass instance """

        if self.allele_source is None:
            raise ValueError("Shouldn't call get_allele_source() on new liftovers")

        return AlleleSource.objects.get_subclass(pk=self.allele_source_id)

    def get_allele_qs(self) -> QuerySet:
        if self.allele_source:
            allele_qs = self.get_allele_source().get_allele_qs()
        else:
            allele_qs = Allele.objects.filter(alleleliftover__liftover=self)
        return allele_qs

    def complete(self):
        liftover_run_complete_signal.send_robust(sender=self.__class__, instance=self)

    def error(self):
        self.alleleliftover_set.update(status=ProcessingStatus.ERROR,
                                       error={"message": f"{self} failed"})

    def __str__(self):
        source = ""
        if self.source_genome_build:
            source = f"from {self.source_genome_build.name} "
        return f"Liftover({self.pk}) {source}to {self.genome_build} via {self.get_conversion_tool_display()}"


class AlleleLiftover(models.Model):
    ERROR_JSON_MESSAGE_KEY = "message"

    liftover = models.ForeignKey(LiftoverRun, on_delete=CASCADE)
    allele = models.ForeignKey(Allele, on_delete=CASCADE)
    # There will only ever be 1 successful AlleleLiftover for a VariantAllele - the one that populated it
    variant_allele = models.OneToOneField(VariantAllele, null=True, blank=True, on_delete=CASCADE)
    status = models.CharField(max_length=1, choices=ProcessingStatus.choices, default=ProcessingStatus.CREATED)
    data = models.JSONField(null=True)  # Per-tool extra info about liftover
    error = models.JSONField(null=True)  # Only set on error - uses ERROR_JSON_MESSAGE_KEY key in dict

    def set_info(self, info: dict):
        print(f"Setting AlleleLiftover({self.pk}): {info=}")
        data = self.data or {}
        data[self._get_data_key()] = info
        self.data = data

    def _get_data_key(self) -> str:
        return slugify(self.liftover.get_conversion_tool_display()) + "-vcf-info"

    def data_tidy(self) -> str | dict:
        if self.data:
            key = self._get_data_key()
            if value := self.data.get(key):
                if INFO_LIFTOVER_SWAPPED_REF_ALT in value:
                    return "Swapped Ref/Alt due to SWAP=1"
                return value
            else:
                return self.data
        return ""

    def error_tidy(self) -> str | dict:
        # If the JSON is just message=, grab the message
        if error_json := self.error:
            if message := error_json.get(self.ERROR_JSON_MESSAGE_KEY):
                if len(error_json.keys()) == 1:
                    return message
            return error_json
        return ""

    class Meta:
        unique_together = ('liftover', 'allele')

    def __str__(self):
        s = f"{self.allele}/{self.liftover}: {self.get_status_display()}"
        if self.error:
            if msg := self.error.get(self.ERROR_JSON_MESSAGE_KEY):
                s += f" error: {msg}"
        return s

    @staticmethod
    def has_existing_failure(allele, dest_genome_build, conversion_tool) -> bool:
        return allele.alleleliftover_set.filter(liftover__genome_build=dest_genome_build,
                                                liftover__conversion_tool=conversion_tool,
                                                status=ProcessingStatus.ERROR).exists()

    @staticmethod
    def get_last_failed_liftover_run(allele, genome_build) -> Optional['LiftoverRun']:
        last_failed_liftover_run = None
        allele_liftover_qs = AlleleLiftover.objects.filter(allele=allele,
                                                           status=ProcessingStatus.ERROR,
                                                           liftover__genome_build=genome_build)
        if al := allele_liftover_qs.order_by("liftover__modified").last():
            last_failed_liftover_run = al.liftover
        return last_failed_liftover_run

    @staticmethod
    def get_unfinished_liftover_run(allele, genome_build) -> Optional['LiftoverRun']:
        unfinished_liftover_run = None
        if al := AlleleLiftover.objects.filter(allele=allele,
                                               status__in=ProcessingStatus.RUNNING_STATES,
                                               liftover__genome_build=genome_build).first():
            unfinished_liftover_run = al.liftover
        return unfinished_liftover_run
