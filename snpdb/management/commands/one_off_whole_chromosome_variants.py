"""
Moves stored whole-chromosome <DUP>/<DEL> Variants to their whole-chromosome form (<CHROM_GAIN>/<CHROM_LOSS> at
position 1, svlen the contig's length - 1), imported before the canonicaliser knew one (#2908,
@see snpdb.whole_chromosome_variants).

Each is changed in place, keeping its pk so genotypes, tags, alleles and classifications follow it. Its
annotation stays as it is: it was already a STRUCTURAL_VARIANT too long for VEP, which is where the
whole-chromosome form goes. Where the whole-chromosome Variant already exists the pair is reported and both are
left alone. --dry-run lists what would change. Registered as a deploy-time ManualOperation
(snpdb/migrations/0019_whole_chromosome_variants.py).
"""
from django.core.management.base import BaseCommand
from django.db import transaction
from django.db.models import F, Q, QuerySet

from library.genomics.vcf_enums import VCFSymbolicAllele
from library.utils import sha256sum_str
from snpdb.models import GenomeBuild, Locus, Sequence, Variant, VariantCoordinate


def _whole_chromosome_candidates(genome_build: GenomeBuild) -> QuerySet[Variant]:
    """ A <DUP>/<DEL> over half its contig - a whole chromosome is always longer """
    half_contig = F("locus__contig__length") / 2
    sv_alts = Sequence.objects.filter(seq__in=[VCFSymbolicAllele.DUP, VCFSymbolicAllele.DEL])
    qs = Variant.objects.filter(Variant.get_contigs_q(genome_build), alt__in=sv_alts)
    qs = qs.filter(Q(svlen__gt=half_contig) | Q(svlen__lt=-half_contig))
    return qs.select_related("locus__contig", "locus__ref", "alt").order_by("pk")


def _get_sequence(seq: str) -> Sequence:
    return Sequence.objects.get_or_create(seq_sha256_hash=sha256sum_str(seq), defaults={"seq": seq})[0]


def _move_variant(variant: Variant, whole_chromosome: VariantCoordinate):
    locus, _ = Locus.objects.get_or_create(contig=variant.locus.contig, position=whole_chromosome.position,
                                           ref=_get_sequence(whole_chromosome.ref))
    Variant.objects.filter(pk=variant.pk).update(locus=locus, alt=_get_sequence(whole_chromosome.alt),
                                                 svlen=whole_chromosome.svlen, end=whole_chromosome.end)


class Command(BaseCommand):
    help = "Move stored whole-chromosome <DUP>/<DEL> variants to <CHROM_GAIN>/<CHROM_LOSS>, keeping their pks"
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report what would change, change nothing")

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        moved = 0
        for genome_build in GenomeBuild.builds_with_annotation():
            for variant in _whole_chromosome_candidates(genome_build):
                whole_chromosome = variant.coordinate.as_internal_canonical_form(genome_build)
                if not whole_chromosome.is_whole_chromosome:
                    continue
                description = f"{genome_build} variant {variant.pk} {variant} -> {whole_chromosome}"
                if existing := Variant.qs_from_variant_coordinate(whole_chromosome, genome_build).first():
                    self.stdout.write(f"{description}: already exists as variant {existing.pk} - left both")
                    continue
                self.stdout.write(f"{description}{' (dry run)' if dry_run else ''}")
                if not dry_run:
                    with transaction.atomic():
                        _move_variant(variant, whole_chromosome)
                moved += 1
        self.stdout.write(f"{moved} variants {'to move' if dry_run else 'moved'}")
