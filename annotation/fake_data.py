"""
Fake annotation: the annotation versions tests build on (get_fake_annotation_version and friends), and the
'create_fake_data' steps "annotation" (a fake annotation version for a build that has none) and "variants"
(annotated variants in the context's genes).
"""
import copy
import os
import random
from dataclasses import dataclass
from typing import Optional
from uuid import uuid4

from django.conf import settings
from django.db.models.fields import IntegerField, TextField
from django.utils import timezone

from annotation.models import ClinVarReviewStatus, GeneAnnotationRelease
from annotation.models.damage_enums import PathogenicityImpact
from annotation.models.models import (
    AnnotationRangeLock,
    AnnotationRun,
    AnnotationVersion,
    ClinVar,
    ClinVarCitation,
    ClinVarCitationsCollection,
    ClinVarVersion,
    GeneAnnotationVersion,
    HumanProteinAtlasAnnotationVersion,
    SubVersionPartition,
    VariantAnnotation,
    VariantAnnotationVersion,
)
from annotation.models.models_citations import CitationIdNormalized, CitationSource
from annotation.models.models_enums import AnnotationStatus
from annotation.vep_config import VEPConfig
from genes.hgvs import HGVSMatcher
from genes.models import GeneAnnotationImport, TranscriptVersion
from genes.models_enums import AnnotationConsortium
from library.django_utils.django_partition import temporary_db_table
from library.fake_data import FakeData, FakeDataContext, register
from library.genomics.vcf_enums import VariantClass
from ontology.tests.test_data_ontology import (
    create_ontology_test_data,
    create_test_ontology_version,
)
from snpdb.models import Variant
from snpdb.models.models_genome import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import (
    slowly_create_loci_and_variants_for_vcf,
    slowly_create_test_variant,
)

# VEP version each columns_version's fixtures were generated against (see the ##VEP= header in
# annotation/tests/test_data/test_columns_version*.vep_annotated.vcf). VEPColumnDefs are gated on
# vep_version as well as columns_version, so this has to be pinned alongside the data files below -
# a developer whose env settings point at an older local VEP install would otherwise silently drop
# columns from the importer bindings and leave the scores unset.
FIXTURE_VEP_VERSIONS = {
    1: "110",
    2: "110",
    3: "112",
    4: "115",
    5: "116",
}


def get_fake_annotation_settings_dict(columns_version: int) -> dict:
    TEST_IMPORT_PROCESSING_DIR = os.path.join(settings.IMPORT_PROCESSING_DIR, "test", str(uuid4()))

    TEST_ANNOTATION = copy.deepcopy(settings.ANNOTATION)
    # phastCons/phyloP custom tracks: v1-v3 fixtures were generated without the bigwig data, so disable
    # to keep importer bindings clean. v4 fixtures include them, so pin the paths (has_data_files gates
    # the columns on a non-None vep_config path) - a developer's local override nulling the bigwig tracks
    # would otherwise drop the columns and leave the scores unset. The files aren't opened: the importer
    # reads the scores straight from the pre-annotated CSQ.
    if columns_version < 4:
        TEST_ANNOTATION[settings.BUILD_GRCH37]["vep_config"].update({
            "phastcons100way": None,
            "phastcons46way": None,
            "phylop100way": None,
            "phylop46way": None,
        })
        TEST_ANNOTATION[settings.BUILD_GRCH38]["vep_config"].update({
            "phastcons100way": None,
            "phastcons30way": None,
            "phylop100way": None,
            "phylop30way": None,
        })
    else:
        TEST_ANNOTATION[settings.BUILD_GRCH37]["vep_config"].update({
            "phastcons100way": "annotation_data/GRCh37/hg19.100way.phastCons.bw",
            "phastcons46way": "annotation_data/GRCh37/hg19.phastCons46way.placental.bw",
            "phylop100way": "annotation_data/GRCh37/hg19.100way.phyloP100way.bw",
            "phylop46way": "annotation_data/GRCh37/hg19.phyloP46way.placental.bw",
        })
        TEST_ANNOTATION[settings.BUILD_GRCH38]["vep_config"].update({
            "phastcons100way": "annotation_data/GRCh38/hg38.phastCons100way.bw",
            "phastcons30way": "annotation_data/GRCh38/hg38.phastCons30way.bw",
            "phylop100way": "annotation_data/GRCh38/hg38.phyloP100way.bw",
            "phylop30way": "annotation_data/GRCh38/hg38.phyloP30way.bw",
        })

    # columns_version 4 fixtures were generated against gnomAD 4.1 (the FILTER column shifts
    # VEPConfig.gnomad4_minor_version, which gates which gnomAD CSQ fields apply); earlier
    # fixtures used gnomAD 4.0. Plugin/version-only data files (denovo_db, mave, etc.) are
    # gated through vep_columns column-defs (min_columns_version), so we don't need to null
    # them out per-version here.
    if columns_version >= 4:
        gnomad4_path = "annotation_data/GRCh38/gnomad4.1_GRCh38_contigs.vcf.gz"
    else:
        gnomad4_path = "annotation_data/GRCh38/gnomad4.0_GRCh38_combined_af.vcf.bgz"

    # Pin gnomAD so a developer's local override doesn't shift VEP CSQ fields and break fixture parsing.
    TEST_ANNOTATION[settings.BUILD_GRCH38]["vep_config"]["gnomad4"] = gnomad4_path

    # COSMIC names the sample count INFO field differently per release (#1673) and the fixtures carry
    # whichever CSQ column the release they were generated against produced - CNT for v1/v2,
    # SAMPLE_COUNT from v3. Pin the release so the columns the importer expects match the fixture.
    if columns_version < 3:
        cosmic_paths = {
            settings.BUILD_GRCH37: "annotation_data/GRCh37/CosmicCodingMuts_v95_20211101_grch37.normal.vcf.gz",
            settings.BUILD_GRCH38: "annotation_data/GRCh38/CosmicCodingMuts_v95_20211101_grch38.normal.vcf.gz",
        }
    else:
        cosmic_paths = {
            settings.BUILD_GRCH37: "annotation_data/GRCh37/Cosmic_GenomeScreensMutant_v99_GRCh37.vcf.gz",
            settings.BUILD_GRCH38: "annotation_data/GRCh38/Cosmic_GenomeScreensMutant_v99_GRCh38.vcf.gz",
        }
    for build_name, cosmic_path in cosmic_paths.items():
        TEST_ANNOTATION[build_name]["vep_config"]["cosmic"] = cosmic_path

    # Same for the columns_version 5 plugin data (#1638) - a deployment pinned below cv5 calls
    # _disable_columns_version_5_plugins(), which nulls these and would drop the columns entirely.
    if columns_version >= 5:
        TEST_ANNOTATION[settings.BUILD_GRCH37]["vep_config"]["protvar"] = \
            "annotation_data/all_builds/ProtVar_data.db"
        TEST_ANNOTATION[settings.BUILD_GRCH38]["vep_config"].update({
            "protvar": "annotation_data/all_builds/ProtVar_data.db",
            "open_targets": "annotation_data/GRCh38/open_targets_26.03_vep.tsv.bgz",
            "eve": "annotation_data/GRCh38/eve_merged.vcf.gz",
            "popeve": "annotation_data/GRCh38/grch38_popEVE_ukbb_20250715.vcf.gz",
            "promoter_ai": "annotation_data/GRCh38/promoterAI_tss500.tsv.bgz",
        })

    ANNOTATION_COLUMNS = copy.deepcopy(TEST_ANNOTATION)
    ANNOTATION_COLUMNS[settings.BUILD_GRCH37]["columns_version"] = columns_version
    ANNOTATION_COLUMNS[settings.BUILD_GRCH38]["columns_version"] = columns_version

    return {
        "IMPORT_PROCESSING_DIR": TEST_IMPORT_PROCESSING_DIR,
        "VARIANT_ZYGOSITY_GLOBAL_COLLECTION": "global",
        "ANNOTATION_VEP_FAKE_VERSION": True,
        "ANNOTATION_VEP_VERSION": FIXTURE_VEP_VERSIONS[columns_version],
        # AnnotSV is off in the shipped defaults - pin it so a developer who enables it locally doesn't
        # trip the SV guards. Tests that want it on override at the method level.
        "ANNOTATION_ANNOTSV_ENABLED": False,
        # Gene level is on in the shipped defaults - pin it so a developer on settings that turn it
        # off still exercises the pipeline.
        "VARIANT_GENE_LEVEL_ENABLED": True,
        "ANNOTATION": ANNOTATION_COLUMNS,
    }


def get_fake_vep_version(genome_build: GenomeBuild, annotation_consortium, columns_version: int):
    # We need to use a later assembly of GRCh37 as the 1st one didn't have MT in it
    if genome_build.name == "GRCh37":
        assembly = "GRCh37.p13"
    else:
        assembly = genome_build.name

    fake_version = {"id": None,
                    "genome_build": genome_build,
                    "assembly": assembly,
                    "annotation_consortium": annotation_consortium,
                    "columns_version": columns_version}
    for f in VariantAnnotationVersion._meta.fields:  # @UndefinedVariable
        if f.name in fake_version:
            continue  # already set
        if isinstance(f, IntegerField):
            value = -1
        elif isinstance(f, TextField):
            # Need proper gnomAD for get_classified_high_frequency_variants_qs
            if f.name == 'gnomad':
                if genome_build.name == 'GRCh37':
                    value = "2.1.1"
                else:
                    value = "3.1"
            elif f.name == 'dbnsfp':
                value = '4.0a'
            else:
                value = "fake"
        else:
            continue
        fake_version[f.name] = value
    # Real COSMIC release from the (test-pinned) vep_config, as the sample count INFO field is
    # gated on it - see the cosmic_count VEPColumnDefs
    fake_version["cosmic"] = VEPConfig(genome_build).cosmic_version
    return fake_version


def get_fake_annotation_version(genome_build: GenomeBuild) -> AnnotationVersion:
    if not settings.UNIT_TEST:
        raise ValueError("Called get_fake_annotation_version while not in a test!")
    return create_fake_annotation_version(genome_build)


def create_fake_annotation_version(genome_build: GenomeBuild) -> AnnotationVersion:
    """ Outside tests only 'create_fake_data annotation' calls this, and only for a build with no annotation """
    gene_annotation_import = GeneAnnotationImport.objects.get_or_create(genome_build=genome_build,
                                                                        annotation_consortium=AnnotationConsortium.ENSEMBL,
                                                                        url="fake")[0]
    # Dotted version, like the real NCBI/Ensembl releases - it ends up in labels and column names
    gene_annotation_release = GeneAnnotationRelease.objects.get_or_create(version="42.20240101",  # TextField
                                                                          genome_build=genome_build,
                                                                          annotation_consortium=AnnotationConsortium.ENSEMBL,
                                                                          defaults={
                                                                              "gene_annotation_import": gene_annotation_import,
                                                                          })[0]

    create_ontology_test_data()
    ontology_version = create_test_ontology_version()

    # Each sub-version save() would otherwise bump AnnotationVersion, so we'd build and discard 4 of
    # them on the way to the one we create below.
    with SubVersionPartition.defer_new_sub_version():
        # gnomad_import_date is a default, not a lookup key - as a key every call created a new row
        gene_annotation_version = GeneAnnotationVersion.objects.get_or_create(gene_annotation_release=gene_annotation_release,
                                                                              ontology_version=ontology_version,
                                                                              defaults={"gnomad_import_date": timezone.now()})[0]

        vav_kwargs = get_fake_vep_version(genome_build, AnnotationConsortium.ENSEMBL, 2)
        vav_kwargs["gene_annotation_release"] = gene_annotation_release
        vav_defaults = {k: vav_kwargs.pop(k) for k in list(vav_kwargs) if k != "genome_build"}
        vav_defaults["status"] = VariantAnnotationVersion.Status.ACTIVE
        variant_annotation_version, _ = VariantAnnotationVersion.objects.get_or_create(
            genome_build=genome_build,
            status=VariantAnnotationVersion.Status.ACTIVE,
            defaults=vav_defaults,
        )
        clinvar_version = ClinVarVersion.objects.get_or_create(filename="fake_clinvar.vcf",
                                                               sha256_hash="not_a_real_hash",
                                                               genome_build=genome_build)[0]
        human_protein_atlas_version = HumanProteinAtlasAnnotationVersion.objects.get_or_create(filename="fake_hpa",
                                                                                               sha256_hash="not_a_real_hash",
                                                                                               hpa_version=0.42)[0]

    av, _ = AnnotationVersion.objects.get_or_create(genome_build=genome_build,
                                                    variant_annotation_version=variant_annotation_version,
                                                    gene_annotation_version=gene_annotation_version,
                                                    clinvar_version=clinvar_version,
                                                    human_protein_atlas_version=human_protein_atlas_version,
                                                    ontology_version=ontology_version)
    return av


def retire_seeded_annotation_version(genome_build: GenomeBuild):
    """ variantgrid/test_runner.py:VariantGridTestRunner seeds an ACTIVE VariantAnnotationVersion per build,
        and one_active_vav_per_build allows only one - so a test that needs its own (a particular consortium,
        columns_version or unpinned fields) retires the seeded one first. """
    VariantAnnotationVersion.objects.filter(genome_build=genome_build,
                                            status=VariantAnnotationVersion.Status.ACTIVE) \
                                    .update(status=VariantAnnotationVersion.Status.HISTORICAL)


def create_fake_variants(genome_build: GenomeBuild):
    build_lc = genome_build.name.lower()
    vcf_filename = os.path.join(settings.BASE_DIR, f"annotation/tests/test_data/test_columns_version1_{build_lc}.vep_annotated.vcf")
    slowly_create_loci_and_variants_for_vcf(genome_build, vcf_filename, get_variant_id_from_info=True)


def create_fake_clinvar_data(clinvar_version: ClinVarVersion):
    create_fake_variants(clinvar_version.genome_build)
    variant = Variant.objects.filter(Variant.get_no_reference_q()).first()

    clinvar_variation_id = 42
    clinvar_allele_id = 42

    defaults = {
        "clinvar_variation_id": clinvar_variation_id,
        "clinvar_allele_id": clinvar_allele_id,
        "preferred_disease_name": "smelly feet",
        "disease_database_name": "blah",
        "review_status": ClinVarReviewStatus.CRITERIA_PROVIDED_MULTIPLE_SUBMITTERS_NO_CONFLICTS,
        "clinical_significance": "Pathogenic",
        "highest_pathogenicity": 5
    }
    ClinVar.objects.get_or_create(version=clinvar_version, variant=variant, defaults=defaults)
    citation = CitationIdNormalized.from_parts(source=CitationSource.PUBMED, index=20613862).for_bulk_create()
    citation.save()
    cvcc, _ = ClinVarCitationsCollection.objects.get_or_create(pk=1)

    ClinVarCitation.objects.get_or_create(clinvar_citations_collection=cvcc,
                                          clinvar_variation_id=clinvar_variation_id,
                                          clinvar_allele_id=clinvar_allele_id,
                                          citation=citation)


def create_fake_variant_annotation(variant, variant_annotation_version: VariantAnnotationVersion) -> VariantAnnotation:
    matcher = HGVSMatcher(variant_annotation_version.genome_build)
    defaults = {
        "hgvs_g": matcher.variant_to_g_hgvs(variant)
        # ??
    }
    annotation_range_lock, _ = AnnotationRangeLock.objects.get_or_create(version=variant_annotation_version,
                                                                         min_variant=variant,
                                                                         max_variant=variant,
                                                                         count=1)
    annotation_run, _ = AnnotationRun.objects.get_or_create(annotation_range_lock=annotation_range_lock)
    va, _ = VariantAnnotation.objects.get_or_create(variant=variant, version=variant_annotation_version,
                                                    annotation_run=annotation_run, defaults=defaults)
    return va


def get_variant_ids_by_gene(genome_build: GenomeBuild, genes: list[str],
                            without_alleles: bool = False) -> dict[str, list[int]]:
    """ without_alleles restricts to variants no allele has claimed yet, so fake data can make its own
        alleles and take them away again without touching anything real """
    gene_symbol_field = "transcript_version__gene_version__gene_symbol_id"
    transcript_versions_qs = TranscriptVersion.objects.filter(genome_build=genome_build,
                                                              gene_version__gene_symbol__in=genes)
    variant_annotation_qs = VariantAnnotation.objects.filter(version=VariantAnnotationVersion.latest(genome_build),
                                                             transcript_version__in=transcript_versions_qs)
    if without_alleles:
        variant_annotation_qs = variant_annotation_qs.filter(variant__variantallele__isnull=True)
    variant_ids_by_gene = {}
    for gene_symbol, variant_id in variant_annotation_qs.values_list(gene_symbol_field, "variant_id"):
        variant_ids_by_gene.setdefault(gene_symbol, []).append(variant_id)
    return variant_ids_by_gene


FAKE_VARIANTS_PIPELINE_COMMAND = "create_fake_data variants"
""" Marks the AnnotationRun the fake variants' annotation hangs off, so delete finds exactly those rows """

COMPLEMENT = str.maketrans("ACGT", "TGCA")
AMINO_ACIDS = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
               "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]


@dataclass(frozen=True)
class FakeVariantKind:
    consequence: str
    impact: str
    variant_class: str
    weight: float


FAKE_VARIANT_KINDS = [
    FakeVariantKind("missense_variant", PathogenicityImpact.MODERATE, VariantClass.SNV, 55),
    FakeVariantKind("synonymous_variant", PathogenicityImpact.LOW, VariantClass.SNV, 22),
    FakeVariantKind("stop_gained", PathogenicityImpact.HIGH, VariantClass.SNV, 7),
    FakeVariantKind("frameshift_variant", PathogenicityImpact.HIGH, VariantClass.DELETION, 9),
    FakeVariantKind("frameshift_variant", PathogenicityImpact.HIGH, VariantClass.INSERTION, 7),
]


def is_fake_annotation_version(variant_annotation_version: VariantAnnotationVersion) -> bool:
    """ What create_fake_annotation_version makes - anything else is real annotation, which fake data never
        writes annotation rows into """
    release = variant_annotation_version.gene_annotation_release
    return release is not None and release.gene_annotation_import.url == "fake"


@register
class FakeAnnotation(FakeData):
    name = "annotation"
    help = "A fake annotation version, when the build has no variant annotation version"
    requires = ("genes",)

    def create(self, context: FakeDataContext, **options):
        genome_build = context.genome_build
        if vav := VariantAnnotationVersion.objects.filter(genome_build=genome_build).order_by("pk").last():
            context.stdout.write(f"{genome_build} already has {vav}")
            return
        annotation_version = create_fake_annotation_version(genome_build)
        context.stdout.write(f"Created fake {annotation_version}")

    def delete(self, context: FakeDataContext, **options):
        context.stdout.write("Leaving the fake annotation version - everything annotated hangs off it")


@register
class FakeVariants(FakeData):
    name = "variants"
    help = ("SNVs and 1bp indels in the genes' coding exons, annotated on a fake annotation version - "
            "on real annotation nothing is made, and later steps use real variants no allele has claimed")
    requires = ("annotation",)

    @classmethod
    def add_arguments(cls, parser):
        parser.add_argument("--variants", type=int, default=300, help="How many, on a fake annotation version")

    def create(self, context: FakeDataContext, **options):
        genome_build = context.genome_build
        vav = VariantAnnotationVersion.latest(genome_build)
        if vav is None:
            raise ValueError(f"{genome_build} has no active variant annotation version")
        if not is_fake_annotation_version(vav):
            context.variant_ids_by_gene = get_variant_ids_by_gene(genome_build, context.genes, without_alleles=True)
            num_variants = sum(len(variant_ids) for variant_ids in context.variant_ids_by_gene.values())
            context.stdout.write(f"{vav} is real annotation - using {num_variants} of its variants")
            return

        fake_annotation_qs = VariantAnnotation.objects.filter(
            version=vav, annotation_run__pipeline_command=FAKE_VARIANTS_PIPELINE_COMMAND)
        if fake_annotation_qs.exists():
            context.stdout.write("Fake variants already exist")
        else:
            self._create_annotated_variants(context, vav, options["variants"])

        for gene_symbol, variant_id in fake_annotation_qs.values_list("symbol", "variant_id"):
            context.variant_ids_by_gene.setdefault(gene_symbol, []).append(variant_id)

    def _create_annotated_variants(self, context: FakeDataContext, vav: VariantAnnotationVersion,
                                   num_variants: int):
        rng = random.Random(context.seed)
        genome_build = context.genome_build
        transcript_versions = list(TranscriptVersion.objects.filter(
            genome_build=genome_build, gene_version__gene_symbol__in=context.genes).select_related(
            "gene_version__gene", "transcript", "contig"))
        if not transcript_versions:
            raise ValueError(f"No transcripts for {context.genes} in {genome_build}")

        exon_sequences = {tv: _exon_sequences(genome_build, tv) for tv in transcript_versions}
        fake_variants = []
        for _ in range(num_variants * 3):  # a spot outside the CDS, or on an N, costs a try
            if len(fake_variants) == num_variants:
                break
            tv = rng.choice(transcript_versions)
            if fake_variant := _random_coding_variant(rng, genome_build, tv, exon_sequences[tv]):
                fake_variants.append(fake_variant)
        variants = [slowly_create_test_variant(fv.chrom, fv.position, fv.ref, fv.alt, genome_build)
                    for fv in fake_variants]

        variant_pks = [variant.pk for variant in variants]
        range_lock = AnnotationRangeLock.objects.create(version=vav, min_variant_id=min(variant_pks),
                                                        max_variant_id=max(variant_pks), count=len(variants))
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=range_lock,
                                                      status=AnnotationStatus.FINISHED,
                                                      pipeline_command=FAKE_VARIANTS_PIPELINE_COMMAND,
                                                      count=len(variants))
        annotations = {}  # One per variant - two fake variants can land on the same spot
        for variant, fake_variant in zip(variants, fake_variants):
            tv = fake_variant.transcript_version
            annotations[variant.pk] = VariantAnnotation(
                version=vav, variant=variant, annotation_run=annotation_run,
                gene=tv.gene_version.gene, transcript=tv.transcript, transcript_version=tv,
                symbol=tv.gene_version.gene_symbol_id, canonical=True,
                consequence=fake_variant.kind.consequence, impact=fake_variant.kind.impact,
                variant_class=fake_variant.kind.variant_class,
                hgvs_g=fake_variant.hgvs_g, hgvs_c=fake_variant.hgvs_c, hgvs_p=fake_variant.hgvs_p,
                gnomad_af=_random_population_af(rng))
        # Into the version's partition, which is the table analysis nodes join to
        representative_table = vav.get_partition_table(VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantAnnotation, representative_table):
            VariantAnnotation.objects.bulk_create(annotations.values())
        context.stdout.write(f"Created {len(annotations)} annotated variants")

    def delete(self, context: FakeDataContext, **options):
        """ The variants themselves stay - they are only coordinates, and anything real may share them """
        range_locks_qs = AnnotationRangeLock.objects.filter(
            annotationrun__pipeline_command=FAKE_VARIANTS_PIPELINE_COMMAND)
        deleted, _ = VariantAnnotation.objects.filter(annotation_run__annotation_range_lock__in=range_locks_qs) \
            .delete()
        range_locks_qs.delete()  # cascades the run
        context.stdout.write(f"Deleted {deleted} fake variant annotations")


@dataclass
class FakeCodingVariant:
    transcript_version: TranscriptVersion
    kind: FakeVariantKind
    chrom: str
    position: int
    ref: str
    alt: str
    hgvs_g: str
    hgvs_c: str
    hgvs_p: str


def _exon_sequences(genome_build: GenomeBuild, transcript_version: TranscriptVersion) -> list[tuple[list, str]]:
    """ (cdot exon, its reference sequence from the base before it) - read whole, so a test run records every
        region the step can read for the sparse CI fastas (@see claude/guides/testing.md) """
    contig_sequence = genome_build.genome_fasta.fasta[transcript_version.contig.name]
    # cdot exons are [0-based genomic start, genomic end, exon number, 1-based transcript start, end, gap]
    return [(exon, contig_sequence[exon[0] - 1:exon[1]].upper())
            for exon in transcript_version.data["genome_builds"][genome_build.name]["exons"]]


def _random_coding_variant(rng: random.Random, genome_build: GenomeBuild, transcript_version: TranscriptVersion,
                           exon_sequences: list[tuple[list, str]]) -> Optional[FakeCodingVariant]:
    """ A variant at a random coding base of the transcript, with its HGVS worked out from the exons.
        None if the spot falls outside the CDS or on an N """
    data = transcript_version.data
    build_data = data["genome_builds"][genome_build.name]
    (alt_start_i, alt_end_i, _exon_number, tx_start, _tx_end, _gap), sequence = rng.choice(exon_sequences)
    position = rng.randint(alt_start_i + 1, alt_end_i)
    plus_strand = build_data["strand"] == "+"
    tx_position = tx_start + (position - 1 - alt_start_i if plus_strand else alt_end_i - position)
    c_position = tx_position - data["start_codon"]
    if not 1 <= c_position <= data["stop_codon"] - data["start_codon"]:
        return None

    ref = sequence[position - alt_start_i]
    anchor = sequence[position - alt_start_i - 1]  # the base before, for a VCF style deletion
    if ref not in "ACGT" or anchor not in "ACGT":
        return None

    def transcript_bases(bases: str) -> str:
        return bases if plus_strand else bases.translate(COMPLEMENT)[::-1]

    kind = rng.choices(FAKE_VARIANT_KINDS, weights=[k.weight for k in FAKE_VARIANT_KINDS])[0]
    codon = (c_position + 2) // 3
    ref_aa = rng.choice(AMINO_ACIDS)
    accession = transcript_version.accession
    contig = transcript_version.contig
    g_prefix = f"{contig.refseq_accession}:g."

    if kind.variant_class == VariantClass.SNV:
        alt = rng.choice([base for base in "ACGT" if base != ref])
        hgvs_g = f"{g_prefix}{position}{ref}>{alt}"
        hgvs_c = f"{accession}:c.{c_position}{transcript_bases(ref)}>{transcript_bases(alt)}"
        protein_change = {"missense_variant": rng.choice([aa for aa in AMINO_ACIDS if aa != ref_aa]),
                          "synonymous_variant": "=", "stop_gained": "Ter"}[kind.consequence]
    elif kind.variant_class == VariantClass.DELETION:
        hgvs_g = f"{g_prefix}{position}del"
        hgvs_c = f"{accession}:c.{c_position}del{transcript_bases(ref)}"
        # VCF style - the deleted base, anchored on the one before it
        ref, alt, position = anchor + ref, anchor, position - 1
        protein_change = "fs"
    else:
        inserted = rng.choice("ACGT")
        alt = ref + inserted
        hgvs_g = f"{g_prefix}{position}_{position + 1}ins{inserted}"
        # the insertion follows the base in genomic order, which is before it on the minus strand
        c_range = f"{c_position}_{c_position + 1}" if plus_strand else f"{c_position - 1}_{c_position}"
        hgvs_c = f"{accession}:c.{c_range}ins{transcript_bases(inserted)}"
        protein_change = "fs"

    return FakeCodingVariant(transcript_version=transcript_version, kind=kind, chrom=contig.name, position=position,
                             ref=ref, alt=alt, hgvs_g=hgvs_g, hgvs_c=hgvs_c,
                             hgvs_p=f"p.{ref_aa}{codon}{protein_change}")


def _random_population_af(rng: random.Random) -> Optional[float]:
    """ Most are rare or absent from gnomAD, a few common enough for a population filter to drop """
    roll = rng.random()
    if roll < 0.3:
        return None
    if roll < 0.9:
        return 10 ** rng.uniform(-5.5, -3)
    return 10 ** rng.uniform(-2, -0.5)
