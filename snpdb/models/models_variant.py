from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.db.models import Value as V, QuerySet
from django.db.models.deletion import CASCADE, DO_NOTHING
from django.db.models.fields import TextField
from django.db.models.functions.text import Concat
from django.db.models.query_utils import Q, FilteredRelation
from django.dispatch import receiver
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from lazy import lazy
from model_utils.managers import InheritanceManager
from typing import Optional, Pattern, Tuple
import collections
import django.dispatch
import re

from flags.models import FlagCollection, flag_collection_extra_info_signal, FlagInfos
from flags.models.models import FlagsMixin, FlagTypeContext
from library.django_utils.django_partition import RelatedModelsPartitionModel
from library.utils import md5sum_str
from snpdb.models import Wiki
from snpdb.models.flag_types import allele_flag_types
from snpdb.models.models_clingen_allele import ClinGenAllele
from snpdb.models.models_genome import Contig, GenomeBuild
from snpdb.models.models_enums import AlleleConversionTool, AlleleOrigin, ProcessingStatus

LOCUS_PATTERN = re.compile(r"^([^:]+):(\d+)[,\s]*([GATC]+)$")
LOCUS_NO_REF_PATTERN = r"^([^:]+):(\d+)$"
VARIANT_PATTERN = re.compile(r"^([^:]+):(\d+)[,\s]*([GATC]+)>([GATC]+)$")
REF_VARIANT_PATTERN = re.compile(r"^([^:]+):(\d+)[,\s]*([GATC]+)>= \(ref\)$")

allele_validate_signal = django.dispatch.Signal(providing_args=["allele"])


class Allele(FlagsMixin, models.Model):
    """ Genome build independent - ie GRCh37 and GRCh38 variants for same change point to same allele
        This is generally done via ClinGen Allele Registry, but sometimes that can fail.
        Linked against Variant with VariantAllele below """

    clingen_allele = models.OneToOneField(ClinGenAllele, null=True, on_delete=CASCADE)

    def get_absolute_url(self):
        # will show allele if there is one, otherwise go to variant page
        return reverse('view_allele', kwargs={"pk": self.id})

    def flag_type_context(self) -> FlagTypeContext:
        return FlagTypeContext.objects.get(pk="allele")

    @lazy
    def clingen_error(self):
        va = self.variantallele_set.filter(error__isnull=False).first()
        error = None
        if va:
            error = va.error
        return error

    def variant_alleles(self):
        return self.variantallele_set.order_by("genome_build__name")

    @lazy
    def grch37(self) -> Optional['Variant']:
        try:
            return self.variant_for_build(genome_build=GenomeBuild.grch37(), best_attempt=False)
        except ValueError:
            return None

    @lazy
    def grch38(self) -> Optional['Variant']:
        try:
            return self.variant_for_build(genome_build=GenomeBuild.grch38(), best_attempt=False)
        except ValueError:
            return None

    @lazy
    def variants(self):
        return Variant.objects.filter(pk__in=self.variant_alleles().values_list('variant', flat=True))

    def variant_for_build(self, genome_build: GenomeBuild, best_attempt=True) -> 'Variant':
        vas = self.variant_alleles()
        va = None
        if genome_build:
            va = vas.filter(genome_build=genome_build).first()
            if not va and not best_attempt:
                raise ValueError(f'Could not find a variant in allele {self.id} for build {genome_build.name}')
        if not va:
            va = vas.first()
        if va:
            return va.variant
        raise ValueError(f'Could not find any variants in allele {self.id}')

    def get_liftover_variant_tuple(self, genome_build: GenomeBuild) -> Tuple[str, 'VariantCoordinate']:
        """ Used by to write VCF coordinates during liftover. Can be slow (API call)
            If you know a VariantAllele exists for your build, use variant_for_build(genome_build).as_tuple() """

        from snpdb.models.models_dbsnp import DbSNP
        from genes.hgvs import get_hgvs_variant_tuple

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
                dbsnp = DbSNP.get_for_variant(va.variant, va.genome_build.latest_variant_annotation_version)
                if dbsnp:
                    g_hgvs = dbsnp.get_g_hgvs(genome_build, alt=va.variant.alt)
                    conversion_tool = AlleleConversionTool.DBSNP

        variant_tuple = None
        if g_hgvs:
            variant_tuple = get_hgvs_variant_tuple(g_hgvs, genome_build)

        return conversion_tool, variant_tuple

    def merge(self, conversion_tool, other_allele: "Allele"):
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
                                      conversion_tool=conversion_tool,
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

            other_fc = other_allele.flag_collection
            other_fc.flag_set.update(collection=self.flag_collection)
            other_fc.flagwatch_set.update(flag_collection=self.flag_collection)
            existing_fc_cc_names = self.flag_collection.clinicalcontext_set.values_list("name", flat=True)
            other_fc.clinicalcontext_set.exclude(name__in=existing_fc_cc_names).update(flag_collection=self.flag_collection)
            other_fc.classification_set.update(flag_collection=self.flag_collection)
            existing_allele_cc_names = self.clinicalcontext_set.values_list("name", flat=True)
            other_allele.clinicalcontext_set.exclude(name__in=existing_allele_cc_names).update(allele=self)
            other_allele.variantallele_set.update(allele=self, conversion_tool=conversion_tool)

    def __str__(self):
        name = f"Allele {self.pk}"
        if self.clingen_allele:
            name += f" ({self.clingen_allele})"

        builds = ", ".join(sorted(self.variantallele_set.values_list("genome_build__name", flat=True)))
        name += f" ({builds})"
        return name

    def validate(self, liftover_complete=True):
        """
        :param liftover_complete: If False does not check for missing representations
        """
        if liftover_complete:
            v37 = self.variant_alleles().filter(genome_build=GenomeBuild.grch37()).first()
            v38 = self.variant_alleles().filter(genome_build=GenomeBuild.grch38()).first()

            if v37:
                self.close_open_flags_of_type(allele_flag_types.missing_37)
            else:
                self.flag_collection_safe.get_or_create_open_flag_of_type(flag_type=allele_flag_types.missing_37, only_if_new=True)

            if v38:
                self.close_open_flags_of_type(allele_flag_types.missing_38)
            else:
                self.flag_collection_safe.get_or_create_open_flag_of_type(flag_type=allele_flag_types.missing_38, only_if_new=True)

        allele_validate_signal.send(sender=Allele, allele=self)


@receiver(flag_collection_extra_info_signal, sender=FlagCollection)
def get_extra_info(flag_infos: FlagInfos, user: User, **kwargs):
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
    conversion_tool = models.CharField(max_length=2, choices=AlleleConversionTool.CHOICES)
    success = models.BooleanField(default=True)
    message = models.TextField(null=True)


VariantCoordinate = collections.namedtuple('VariantCoordinate', 'chrom pos ref alt')


class Sequence(models.Model):
    """
        Issue #2140: We want to guarantee seq is unique (so Locus/Variant can have unique constraints)
        Postgres by default uses indexes for constraints, and large TextFields give error of:
        "index row requires x bytes, maximum size is 8191"

        The easiest solution is to md5sum seq and make the constraint on that. Another possible solution is to use
        Gist indexes but that requires installing the btree_gist extension (requires postgres Admin rights).
        Django 3 has ExclusionConstraint, Postgres contrib has BtreeGistExtension to add via migration
    """
    seq = models.TextField()
    seq_md5_hash = models.CharField(max_length=32, unique=True)
    length = models.IntegerField()

    def save(self, *args, **kwargs):
        if not self.seq_md5_hash:
            self.seq_md5_hash = md5sum_str(self.seq)
        super().save(*args, **kwargs)

    @staticmethod
    def abbreviate(s: str, max_length: int = 20):
        if len(s) > max_length:
            s = f"{s[:3]}...{s[-3:]}"
        return s

    def __str__(self):
        return self.abbreviate(self.seq)

    @staticmethod
    def get_pk_by_seq(q=None):
        qs = Sequence.objects.all()
        if q:
            qs = qs.filter(q)
        return dict(qs.values_list("seq", "pk"))

    def is_standard_sequence(self):
        """ only contains G/A/T/C/N """
        return not re.match(r"[^GATCN]", self.seq)


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


class Variant(models.Model):
    """ Variants represent the different alleles at a locus
        Usually 2+ per line in a VCF file (ref + >= 1 alts pointing to the same locus for the row)
        There is only 1 Variant for a given locus/alt per database (handled via insertion queues) """

    REFERENCE_ALT = "="
    locus = models.ForeignKey(Locus, on_delete=CASCADE)
    alt = models.ForeignKey(Sequence, on_delete=CASCADE)

    class Meta:
        unique_together = ("locus", "alt")

    @staticmethod
    def get_chrom_q(chrom):
        return Q(locus__contig__name__iexact=chrom) | Q(locus__contig__ucsc_name__iexact=chrom)

    @staticmethod
    def get_contigs_q(genome_build: GenomeBuild):
        return Q(locus__contig__in=genome_build.contigs)

    @staticmethod
    def get_no_reference_q():
        return ~Q(alt__seq=Variant.REFERENCE_ALT)

    @staticmethod
    def annotate_variant_string(qs, name="variant_string", path_to_variant=""):
        """ Return a "1:123321 G>C" style string in a query """
        kwargs = {name: Concat(f"{path_to_variant}locus__contig__name", V(":"),
                               f"{path_to_variant}locus__position", V(" "),
                               f"{path_to_variant}locus__ref__seq", V(">"),
                               f"{path_to_variant}alt__seq", output_field=TextField())}
        return qs.annotate(**kwargs)

    @staticmethod
    def format_tuple(chrom, position, ref, alt, abbreviate=False) -> str:
        if abbreviate:
            ref = Sequence.abbreviate(ref)
            alt = Sequence.abbreviate(alt)
        return f"{chrom}:{position} {ref}>{alt}"

    @staticmethod
    def get_tuple_from_string(variant_string: str, genome_build: GenomeBuild,
                              regex_pattern: Pattern[str] = VARIANT_PATTERN) -> VariantCoordinate:
        """ regex_pattern - has to have 4 groups, returns (chrom, position, ref, alt) """
        variant_tuple = None
        if m := regex_pattern.match(variant_string):
            chrom, position, ref, alt = m.groups()
            contig = genome_build.chrom_contig_mappings[chrom]
            variant_tuple = VariantCoordinate(contig.name, int(position), ref, alt)
        return variant_tuple

    @staticmethod
    def get_from_string(variant_string: str, genome_build: GenomeBuild,
                        regex_pattern=VARIANT_PATTERN) -> Optional['Variant']:
        variant_tuple = Variant.get_tuple_from_string(variant_string, genome_build, regex_pattern=regex_pattern)
        try:
            return Variant.get_from_tuple(variant_tuple, genome_build)
        except Variant.DoesNotExist:
            return None

    @staticmethod
    def get_from_tuple(variant_tuple: VariantCoordinate, genome_build: GenomeBuild) -> 'Variant':
        params = ["locus__contig__name", "locus__position", "locus__ref__seq", "alt__seq"]
        return Variant.objects.get(locus__contig__genomebuildcontig__genome_build=genome_build,
                                   **dict(zip(params, variant_tuple)))

    @lazy
    def genome_build(self) -> Optional['GenomeBuild']:
        try:  # Short cut
            return self.variantallele.genome_build
        except VariantAllele.DoesNotExist:
            pass

        q = Q(genome_build__name__in=settings.ANNOTATION)
        build_contig = self.locus.contig.genomebuildcontig_set.filter(q).first()
        if build_contig:
            return build_contig.genome_build
        return None

    @lazy
    def coordinate(self) -> VariantCoordinate:
        locus = self.locus
        contig = locus.contig
        return VariantCoordinate(chrom=contig.name, pos=locus.position, ref=locus.ref.seq, alt=self.alt.seq)

    @property
    def is_reference(self) -> bool:
        return self.alt.seq == self.REFERENCE_ALT

    @property
    def is_standard_variant(self) -> bool:
        """ Variant alt sequence is standard [GATCN] (ie not special or reference) """
        # locus.ref should always be standard...
        return self.alt.is_standard_sequence()

    @property
    def is_deletion(self) -> bool:
        return self.alt.seq != Variant.REFERENCE_ALT and self.locus.ref.length > self.alt.length

    @property
    def can_have_annotation(self) -> bool:
        return self.is_standard_variant

    def as_tuple(self) -> VariantCoordinate:
        return self.locus.contig.name, self.locus.position, self.locus.ref.seq, self.alt.seq

    def is_abbreviated(self):
        return str(self) != self.full_string

    @lazy
    def full_string(self):
        """ No abbreviation """
        return self.format_tuple(*self.as_tuple())

    def __str__(self):
        return self.format_tuple(self.locus.contig.name, self.locus.position, self.locus.ref, self.alt)

    def get_absolute_url(self):
        # will show allele if there is one, otherwise go to variant page
        return reverse('view_allele_from_variant', kwargs={"variant_id": self.pk})

    @lazy
    def allele(self) -> Optional[Allele]:
        va = VariantAllele.objects.filter(variant=self).first()
        if va:
            return va.allele
        return None

    @property
    def equivalent_variants(self):
        allele = self.allele
        if not allele:
            return self
        return VariantAllele.objects.filter(allele=allele).values_list('variant', flat=True)

    @lazy
    def variant_annotation(self) -> Optional['VariantAnnotation']:
        vav = self.genome_build.latest_variant_annotation_version
        return self.variantannotation_set.filter(version=vav).first()

    @lazy
    def canonical_transcript_annotation(self) -> Optional['VariantTranscriptAnnotation']:
        vav = self.genome_build.latest_variant_annotation_version
        return self.varianttranscriptannotation_set.filter(version=vav, canonical=True).first()

    @lazy
    def canonical_c_hgvs(self):
        c_hgvs = None
        if self.canonical_transcript_annotation:
            c_hgvs = self.canonical_transcript_annotation.hgvs_c
        return c_hgvs


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

    variant = models.OneToOneField(Variant, on_delete=CASCADE)
    genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE)
    allele = models.ForeignKey(Allele, on_delete=CASCADE)
    origin = models.CharField(max_length=1, choices=AlleleOrigin.CHOICES)
    conversion_tool = models.CharField(max_length=2, choices=AlleleConversionTool.CHOICES)
    error = models.JSONField(null=True)  # Only set on error

    def __str__(self):
        return f"{self.allele} - {self.variant_id}({self.genome_build}/{self.conversion_tool})"


class VariantCollection(RelatedModelsPartitionModel):
    """ A set of variants - usually used as a cached result """

    RECORDS_BASE_TABLE_NAMES = ["snpdb_variantcollectionrecord"]
    RECORDS_FK_FIELD_TO_THIS_MODEL = "variant_collection_id"
    PARTITION_LABEL_TEXT = "variant_collection"
    name = models.TextField(null=True)
    count = models.IntegerField(null=True)
    status = models.CharField(max_length=1, choices=ProcessingStatus.CHOICES, default=ProcessingStatus.CREATED)

    @property
    def variant_collection_alias(self):
        return f"variantcollection_{self.pk}"

    def get_annotation_kwargs(self):
        vcr_condition = Q(variantcollectionrecord__variant_collection=self)
        return {self.variant_collection_alias: FilteredRelation('variantcollectionrecord', condition=vcr_condition)}

    def get_q(self):
        if self.status != ProcessingStatus.SUCCESS:
            raise ValueError(f"{self}: status {self.get_status_display()} != SUCCESS")

        return Q(**{f"{self.variant_collection_alias}__isnull": False})

    def __str__(self):
        return f"VariantCollection: {self.pk} ({self.name})"


class VariantCollectionRecord(models.Model):
    variant_collection = models.ForeignKey(VariantCollection, on_delete=DO_NOTHING)  # handled via drop partition
    variant = models.ForeignKey(Variant, on_delete=CASCADE)


class AlleleSource(models.Model):
    """ Provides a source of alleles for liftover pipelines. """
    objects = InheritanceManager()

    def get_allele_qs(self):
        return Allele.objects.none()

    def liftover_complete(self, genome_build: GenomeBuild):
        """ This is called at the end of a liftover pipeline (once per build) """
        pass


class Liftover(TimeStampedModel):
    """ Liftover pipeline involves reading through a VCF where ID is set to Allele.pk and then creating
        VariantAllele entries for the variant/allele

        Some AlleleConversionTools (eg ClinGen AlleleRegistry) we can write the VCF in the desired genome build
        For others (NCBI Remap) we need to write the source genome build VCF first

        Alleles must have already been created - allele_source used to retrieve them

        The VCF (in genome_build build) is set in UploadedFile for the UploadPipeline """
    user = models.ForeignKey(User, on_delete=CASCADE)
    allele_source = models.ForeignKey(AlleleSource, on_delete=CASCADE)
    conversion_tool = models.CharField(max_length=2, choices=AlleleConversionTool.CHOICES)
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
        return f"Liftover {source}to {self.genome_build.name} via {self.get_conversion_tool_display()}"


class LiftoverError(models.Model):
    liftover = models.ForeignKey(Liftover, on_delete=CASCADE)
    allele = models.ForeignKey(Allele, on_delete=CASCADE)
    variant = models.ForeignKey(Variant, null=True, on_delete=CASCADE)  # Optional, if got a variant but invalid
    error_message = models.TextField()

    class Meta:
        unique_together = ('liftover', 'allele')

    def __str__(self):
        return f"{self.allele} failed {self.liftover}: {self.error_message}"
