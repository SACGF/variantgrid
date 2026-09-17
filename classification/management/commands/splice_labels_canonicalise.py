"""Re-point existing splice Variants at the canonical label for their junction.

A splice event's label is canonical - lower-case tokens joined by underscores (@see
genes.gene_splice.canonical_splice_label) - so every way of writing one junction's name is one
Variant. Variants loaded before that hold the label as it was written (V7, VIII, EX14SKIP) and a
junction named by its breakpoints holds them without the build that makes them mean something.

Only the alt moves: the Variant keeps its pk, so classifications, VariantAllele, samples and
annotation follow it. The ImportedAlleleInfo coordinate strings that name the old alt are rewritten
with it, so a re-match finds the same Variant. A canonical alt some other Variant already holds is
reported rather than merged - merging two Variants is not this command's job.

Writes snpdb_variant, so run it with a person watching (--dry-run first).
"""
from django.core.management.base import BaseCommand
from django.db import transaction
from django.db.models import Value
from django.db.models.functions import Replace

from classification.models import ImportedAlleleInfo
from genes.gene_splice import canonical_splice_label
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from library.utils import sha256sum_str
from snpdb.models import CohortGenotypeCollection, GenomeBuild, Sequence, Variant

SPLICE_ALT_PREFIX = f"<{GeneLevelSymbolicAlt.SPLICE}:"


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help="List what would be re-pointed, and change nothing")

    def handle(self, *args, **options):
        dry_run = options['dry_run']
        changed = 0
        for variant in splice_variants():
            kind, namespace, gene_id, label = GeneLevelSymbolicAlt.parse(variant.alt.seq)
            canonical = canonical_label_for(variant, label)
            if canonical is None:
                print(f"{variant.pk}\t{variant.alt.seq}\tno canonical label - the build its breakpoints "
                      f"were called in could not be determined from the VCFs it was loaded from")
                continue

            new_alt = GeneLevelSymbolicAlt.format(kind, namespace, gene_id, canonical)
            old_alt = variant.alt.seq
            if new_alt == old_alt:
                continue

            clash = Variant.objects.filter(locus=variant.locus, alt__seq=new_alt, svlen=variant.svlen) \
                                   .exclude(pk=variant.pk).first()
            if clash:
                print(f"{variant.pk}\t{old_alt}\t-> {new_alt} already held by variant {clash.pk} - "
                      f"left for you to merge")
                continue

            allele_infos = ImportedAlleleInfo.objects.filter(variant_coordinate__contains=old_alt)
            print(f"{variant.pk}\t{old_alt}\t-> {new_alt}\t"
                  f"{allele_infos.count()} imported allele info coordinate(s)")
            if not dry_run:
                repoint(variant, new_alt, old_alt)
            changed += 1

        print(f"{'Would re-point' if dry_run else 'Re-pointed'} {changed} splice variant(s)")


def splice_variants():
    """ Every splice Variant, oldest first - the alt is how a splice event is stored """
    variant_qs = Variant.objects.filter(Variant.get_gene_level_q(), alt__seq__startswith=SPLICE_ALT_PREFIX)
    for variant in variant_qs.select_related("locus", "alt").order_by("pk"):
        if GeneLevelSymbolicAlt.parse(variant.alt.seq):
            yield variant


def canonical_label_for(variant: Variant, label: str):
    """ The canonical form of a label read off an alt. Only a label made of breakpoints needs a
        build, which is the one the VCF it was loaded from was called in """
    if canonical := canonical_splice_label(label):
        return canonical
    genome_builds = genome_builds_for_variant(variant)
    if len(genome_builds) != 1:
        return None
    return canonical_splice_label(label, genome_builds.pop())


def genome_builds_for_variant(variant: Variant) -> set[GenomeBuild]:
    """ The builds of the VCFs this Variant was loaded from - a gene-level Variant sits on the contig
        every build shares, so the VCF is the only record of which build its coordinates came from """
    collections = CohortGenotypeCollection.objects.filter(cohortgenotype__variant=variant,
                                                          cohort__vcf__isnull=False) \
                                                  .select_related("cohort__vcf__genome_build")
    return {c.cohort.vcf.genome_build for c in collections if c.cohort.vcf.genome_build}


@transaction.atomic
def repoint(variant: Variant, new_alt: str, old_alt: str):
    """ The Variant keeps its pk, so everything linked to it follows """
    sequence, _ = Sequence.objects.get_or_create(seq=new_alt,
                                                 defaults={"seq_sha256_hash": sha256sum_str(new_alt)})
    variant.alt = sequence
    variant.save(update_fields=["alt"])
    ImportedAlleleInfo.objects.filter(variant_coordinate__contains=old_alt) \
                              .update(variant_coordinate=Replace("variant_coordinate",
                                                                 Value(old_alt), Value(new_alt)))
