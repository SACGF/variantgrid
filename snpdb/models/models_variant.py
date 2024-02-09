import logging
import re
from functools import cached_property
from typing import Optional, Iterable, Union, Any

import pydantic
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models, IntegrityError
from django.db.models import Value, QuerySet, F
from django.db.models.deletion import CASCADE, DO_NOTHING
from django.db.models.fields import TextField
from django.db.models.functions import Greatest
from django.db.models.functions.text import Concat
from django.db.models.query_utils import Q, FilteredRelation
from django.dispatch import receiver
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager

from flags.models import FlagCollection, flag_collection_extra_info_signal, FlagInfos
from flags.models.models import FlagsMixin, FlagTypeContext
from library.django_utils.django_object_managers import ObjectManagerCachingRequest
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.genomics import format_chrom
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.preview_request import PreviewModelMixin, PreviewKeyValue
from library.utils import md5sum_str, FormerTuple
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
        if va := self.variantallele_set.filter(error__isnull=False).first():
            error = va.error
        return error

    def variant_alleles(self) -> QuerySet['VariantAllele']:
        return self.variantallele_set.select_related('variant__locus', 'variant__locus__contig', 'variant__locus__ref', 'variant__alt').order_by("genome_build__name")

    @cached_property
    def grch37(self) -> Optional['Variant']:
        try:
            return self.variant_for_build(genome_build=GenomeBuild.grch37(), best_attempt=False)
        except ValueError:
            return None

    @cached_property
    def grch38(self) -> Optional['Variant']:
        try:
            return self.variant_for_build(genome_build=GenomeBuild.grch38(), best_attempt=False)
        except ValueError:
            return None

    @cached_property
    def variants(self):
        return Variant.objects.filter(pk__in=self.variant_alleles().values_list('variant', flat=True))

    def variant_for_build(self, genome_build: GenomeBuild, best_attempt=True) -> 'Variant':
        vas = self.variant_alleles()
        va = None
        if genome_build:
            va = vas.filter(genome_build=genome_build).first()
            if not va and not best_attempt:
                raise ValueError(f'Could not find a variant in allele {self.id} for build {genome_build}')
        if not va:
            va = vas.first()
        if va:
            return va.variant
        raise ValueError(f'Could not find any variants in allele {self.id}')

    def variant_for_build_optional(self, genome_build: GenomeBuild) -> Optional['Variant']:
        if va := self.variant_alleles().filter(genome_build=genome_build).first():
            return va.variant

    def get_liftover_tuple(self, genome_build: GenomeBuild,
                           hgvs_matcher=None) -> tuple[AlleleConversionTool, Union[int, 'VariantCoordinate']]:
        """ Used by to write VCF coordinates during liftover. Can be slow (API call)
            If you know a VariantAllele exists for your build, use variant_for_build(genome_build).as_tuple()

            Optionally pass in hgvs_matcher to save re-instantiating it all the time """

        from annotation.models import VariantAnnotationVersion
        from snpdb.models.models_dbsnp import DbSNP
        from genes.hgvs import get_hgvs_variant_coordinate

        # Check if the other build shares existing contig
        genome_build_contigs = set(c.pk for c in genome_build.chrom_contig_mappings.values())
        for variant_allele in self.variantallele_set.all():
            if variant_allele.variant.locus.contig_id in genome_build_contigs:
                conversion_tool = AlleleConversionTool.SAME_CONTIG
                # Return variant_id so we can create it directly
                return conversion_tool, variant_allele.variant_id

        conversion_tool = None
        g_hgvs = None
        if self.clingen_allele:
            try:
                g_hgvs = self.clingen_allele.get_g_hgvs(genome_build)
                conversion_tool = AlleleConversionTool.CLINGEN_ALLELE_REGISTRY
            except ValueError:  # Various contig errors all subclass from this
                pass
        if g_hgvs is None:
            if settings.LIFTOVER_DBSNP_ENABLED:
                va = self.variantallele_set.all().first()
                if va is None:
                    raise ValueError("Allele contains no VariantAlleles at all! Cannot liftover")
                if dbsnp := DbSNP.get_for_variant(va.variant, VariantAnnotationVersion.latest(va.genome_build)):
                    g_hgvs = dbsnp.get_g_hgvs(genome_build, alt=va.variant.alt)
                    conversion_tool = AlleleConversionTool.DBSNP

        variant_coordinate = None
        if g_hgvs:
            if hgvs_matcher:
                variant_coordinate = hgvs_matcher.get_variant_coordinate(g_hgvs)
            else:
                variant_coordinate = get_hgvs_variant_coordinate(g_hgvs, genome_build)

        return conversion_tool, variant_coordinate

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
                    va.error = None  # clear any errors
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
    start: int
    end: int
    ref: str
    alt: str

    @property
    def as_tuple(self) -> tuple:
        return self.chrom, self.start, self.end, self.ref, self.alt

    def __str__(self):
        return self.format()

    def format(self):
        if Sequence.allele_is_symbolic(self.alt):
            return f"{self.chrom}:{self.start}-{self.end} {self.alt}"
        else:
            return f"{self.chrom}:{self.start} {self.ref}>{self.alt}"

    @staticmethod
    def from_variant_match(match, genome_build=None):
        chrom = match.group(1)
        if genome_build:
            chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        start = int(match.group(2))
        ref = match.group(3).strip().upper()
        alt = match.group(4).strip().upper()
        return VariantCoordinate.from_start_only(chrom, start, ref, alt)

    @staticmethod
    def from_symbolic_match(match, genome_build=None):
        chrom, start, end, alt = match.groups()
        if genome_build:
            chrom = format_chrom(chrom, genome_build.reference_fasta_has_chr)
        start = int(start)
        end = int(end)
        alt = "<" + alt + ">"  # captured what was inside of brackets ie <()>
        if genome_build:
            contig_sequence = genome_build.genome_fasta.fasta[chrom]
            ref = contig_sequence[start:start + 1].upper()
        else:
            ref = "N"
        return VariantCoordinate(chrom=chrom, start=start, end=end, ref=ref, alt=alt)

    @staticmethod
    def from_string(variant_string: str, genome_build=None):
        """ Pass in genome build to be able to set REF from symbolic (will be N otherwise) """
        if full_match := VARIANT_PATTERN.fullmatch(variant_string):
            return VariantCoordinate.from_variant_match(full_match, genome_build)
        elif full_match := VARIANT_SYMBOLIC_PATTERN.fullmatch(variant_string):
            return VariantCoordinate.from_symbolic_match(full_match, genome_build)
        regex_patterns = ", ".join((str(s) for s in (VARIANT_PATTERN, VARIANT_SYMBOLIC_PATTERN)))
        raise ValueError(f"{variant_string=} did not match against {regex_patterns=}")

    @staticmethod
    def from_start_only(chrom: str, start: int, ref: str, alt: str):
        """ Initialise w/o providing an end """

        if Sequence.allele_is_symbolic(alt):
            raise ValueError("Must pass 'end' when using symbolic alt")
        end = start + abs(len(ref) - len(alt))
        return VariantCoordinate(chrom=chrom, start=start, end=end, ref=ref, alt=alt)

    def is_symbolic(self):
        return Sequence.allele_is_symbolic(self.alt)

    def as_external_explicit(self, genome_build) -> 'VariantCoordinate':
        """ explicit ref/alt """
        if self.is_symbolic():
            contig_sequence = genome_build.genome_fasta.fasta[self.chrom]
            ref_sequence = contig_sequence[self.start-1:self.end].upper()
            if self.alt == VCFSymbolicAllele.DEL:
                ref = ref_sequence
                alt = ref_sequence[0]
            elif self.alt == VCFSymbolicAllele.DUP:
                ref = ref_sequence[0]
                alt = ref_sequence
            else:
                raise ValueError(f"Unknown symbolic alt of '{self.alt}'")

            return VariantCoordinate(chrom=self.chrom, start=self.start, end=self.end, ref=ref, alt=alt)

        if self.alt == Variant.REFERENCE_ALT:
            # Convert from our internal format (alt='=' for ref) to explicit
            alt = self.ref
        else:
            alt = self.alt
        return VariantCoordinate(chrom=self.chrom, start=self.start, end=self.end, ref=self.ref, alt=alt)

    def as_internal_symbolic(self):
        """ Internal format - alt can be <DEL> or <DUP>
            Uses our internal reference representation
        """
        if self.is_symbolic():
            return self

        ref = self.ref
        if self.alt == self.ref:
            alt = Variant.REFERENCE_ALT
        else:
            ref_length = len(ref)
            alt_length = len(self.alt)
            diff = alt_length - ref_length
            svlen = abs(diff)
            if svlen >= settings.VARIANT_SYMBOLIC_ALT_SIZE:
                # TODO Check against existing end?
                end = self.start + svlen
                if self.end != end:
                    raise ValueError(f"{self}: end={self.end}, calculated end={end}")

                if diff > 0:
                    alt = VCFSymbolicAllele.DUP
                else:
                    ref = self.ref[0]
                    alt = VCFSymbolicAllele.DEL
            else:
                alt = self.alt
        return VariantCoordinate(chrom=self.chrom, start=self.start, end=self.end, ref=ref, alt=alt)


class Sequence(models.Model):
    """
        We want to guarantee seq is unique (so Locus/Variant can have unique constraints)
        Postgres by default uses indexes for constraints, and large TextFields give error of:
        "index row requires x bytes, maximum size is 8191"

        The easiest solution is to md5sum seq and make the constraint on that. Another possible solution is to use
        Gist indexes but that requires installing the btree_gist extension (requires postgres Admin rights).
        Django 3 has ExclusionConstraint, Postgres contrib has BtreeGistExtension to add via migration
    """
    seq = models.TextField()
    seq_md5_hash = models.CharField(max_length=32, unique=True)
    length = models.IntegerField()  # TODO: I think we should remove this as we now have Variant.end

    def save(self, force_insert=False, force_update=False, using=None, update_fields=None):
        if not self.seq_md5_hash:
            self.seq_md5_hash = md5sum_str(self.seq)
        super().save(force_insert=force_insert, force_update=force_update, using=using, update_fields=update_fields)

    @staticmethod
    def abbreviate(s: str, max_length: int = 20):
        if len(s) > max_length:
            s = f"{s[:3]}...{s[-3:]}"
        return s

    def __len__(self) -> int:
        return self.length

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
    locus = models.ForeignKey(Locus, on_delete=CASCADE)
    alt = models.ForeignKey(Sequence, on_delete=CASCADE)
    # end depends on length of ref/alt so can't be on a locus
    end = models.IntegerField(null=True)

    class Meta:
        unique_together = ("locus", "alt", "end")  # Possible to have eg CNV with same alt = <INS> but diff end

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
    def get_no_reference_q():
        return ~Q(alt__seq=Variant.REFERENCE_ALT)

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
    def format_tuple(chrom, start, end, ref, alt, abbreviate=False) -> str:
        if abbreviate and not Sequence.allele_is_symbolic(alt):
            ref = Sequence.abbreviate(ref)
            alt = Sequence.abbreviate(alt)
        vc = VariantCoordinate(chrom=chrom, start=start, end=end, ref=ref, alt=alt)
        return vc.format()

    @staticmethod
    def get_from_string(variant_string: str, genome_build: GenomeBuild) -> Optional['Variant']:
        variant_coordinate = VariantCoordinate.from_string(variant_string)
        try:
            return Variant.get_from_variant_coordinate(variant_coordinate, genome_build)
        except Variant.DoesNotExist:
            return None

    @staticmethod
    def get_from_variant_coordinate(variant_coordinate: VariantCoordinate, genome_build: GenomeBuild) -> 'Variant':
        variant_coordinate = variant_coordinate.as_internal_symbolic()
        params = ["locus__contig__name", "locus__position", "end", "locus__ref__seq", "alt__seq"]
        return Variant.objects.get(locus__contig__genomebuildcontig__genome_build=genome_build,
                                   **dict(zip(params, variant_coordinate)))

    @cached_property
    def genome_builds(self) -> set['GenomeBuild']:
        gbc_qs = GenomeBuildContig.objects.filter(genome_build__in=GenomeBuild.builds_with_annotation(),
                                                  contig__locus__variant=self)
        return {gbc.genome_build for gbc in gbc_qs}

    @cached_property
    def coordinate(self) -> VariantCoordinate:
        locus = self.locus
        contig = locus.contig
        return VariantCoordinate(chrom=contig.name, start=locus.position, end=self.end,
                                 ref=locus.ref.seq, alt=self.alt.seq)

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
        return self.alt.seq != Variant.REFERENCE_ALT and self.locus.ref.length < self.alt.length

    @property
    def is_deletion(self) -> bool:
        if self.alt.is_symbolic:
            return self.alt.seq == VCFSymbolicAllele.DEL
        return self.alt.seq != Variant.REFERENCE_ALT and self.locus.ref.length > self.alt.length

    @property
    def is_symbolic(self) -> bool:
        return self.locus.ref.is_symbolic() or self.alt.is_symbolic()

    @property
    def can_make_g_hgvs(self) -> bool:
        """ Can't form ones with some symbolic variants (eg <INS>) """
        if self.is_symbolic:
            return self.alt.seq in {VCFSymbolicAllele.DEL, VCFSymbolicAllele.DUP}
        return True

    @property
    def can_have_clingen_allele(self) -> bool:
        return self.length <= settings.CLINGEN_ALLELE_REGISTRY_MAX_ALLELE_SIZE and self.can_make_g_hgvs
    @property
    def can_have_annotation(self) -> bool:
        return not self.is_reference

    def as_tuple(self) -> tuple[str, int, int, str, str]:
        return self.locus.contig.name, self.locus.position, self.end, self.locus.ref.seq, self.alt.seq

    def is_abbreviated(self):
        return str(self) != self.full_string

    @cached_property
    def full_string(self):
        """ No abbreviation """
        return self.format_tuple(*self.as_tuple())

    def __str__(self):
        return self.format_tuple(self.locus.contig.name, self.locus.position, self.end, self.locus.ref.seq, self.alt.seq)

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
            c_hgvs = cta.hgvs_c
        return c_hgvs

    @property
    def start(self):
        return self.locus.position

    @property
    def length(self) -> int:
        return self.end - self.locus.position

    @staticmethod
    def calculate_end_from_lengths(position, ref_length, alt_length) -> int:
        return position + abs(ref_length - alt_length)

    @staticmethod
    def calculate_end(position, ref, alt) -> int:
        return Variant.calculate_end_from_lengths(position, len(ref), len(alt))

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
    error = models.JSONField(null=True)  # Only set on error

    class Meta:
        unique_together = ("variant", "genome_build", "allele")

    @property
    def canonical_c_hgvs(self):
        return self.variant.get_canonical_c_hgvs(self.genome_build)

    def needs_clingen_call(self):
        if settings.CLINGEN_ALLELE_REGISTRY_LOGIN and self.allele.clingen_allele is None:
            if self.error:
                # Retry if server was down
                return self.error.get("errorType") == ClinGenAllele.CLINGEN_ALLELE_SERVER_ERROR_TYPE
            return True
        return False

    def __str__(self):
        return f"{self.allele} - {self.variant_id}({self.genome_build}/{self.allele_linking_tool})"


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

    @staticmethod
    def get_liftover_for_allele(allele, genome_build) -> Optional['Liftover']:
        """ Only works if liftover was done via VariantAlleleSource """
        allele_sources_qs = VariantAlleleSource.objects.filter(variant_allele__allele=allele)
        return Liftover.objects.filter(allele_source__in=allele_sources_qs, genome_build=genome_build).first()


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


class Liftover(TimeStampedModel):
    """ Liftover pipeline involves reading through a VCF where ID is set to Allele.pk and then creating
        VariantAllele entries for the variant/allele

        Some AlleleConversionTools (eg ClinGen AlleleRegistry) we can write the VCF in the desired genome build
        For others (NCBI Remap) we need to write the source genome build VCF first

        Alleles must have already been created - allele_source used to retrieve them

        The VCF (in genome_build build) is set in UploadedFile for the UploadPipeline """
    user = models.ForeignKey(User, on_delete=CASCADE)
    allele_source = models.ForeignKey(AlleleSource, on_delete=CASCADE)
    conversion_tool = models.CharField(max_length=2, choices=AlleleConversionTool.choices)
    source_vcf = models.TextField(null=True)
    source_genome_build = models.ForeignKey(GenomeBuild, null=True, on_delete=CASCADE,
                                            related_name="liftover_source_genome_build")
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)  # destination

    def get_allele_source(self) -> AlleleSource:
        """ Returns subclass instance """
        return AlleleSource.objects.get_subclass(pk=self.allele_source_id)

    def get_allele_qs(self) -> QuerySet:
        return self.get_allele_source().get_allele_qs()

    def complete(self):
        self.get_allele_source().liftover_complete(genome_build=self.genome_build)

    def __str__(self):
        source = ""
        if self.source_genome_build:
            source = f"from {self.source_genome_build.name} "
        return f"Liftover {source}to {self.genome_build} via {self.get_conversion_tool_display()}"


class LiftoverError(models.Model):
    liftover = models.ForeignKey(Liftover, on_delete=CASCADE)
    allele = models.ForeignKey(Allele, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, null=True, on_delete=CASCADE)  # Optional, if got a variant but invalid
    error_message = models.TextField()

    class Meta:
        unique_together = ('liftover', 'allele')

    def __str__(self):
        return f"{self.allele} failed {self.liftover}: {self.error_message}"
