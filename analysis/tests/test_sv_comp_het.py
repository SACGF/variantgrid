"""Compound het where a structural variant crosses multiple genes (issue #940).

A long SV is skipped by VEP (TOO_LONG) so it has no VariantTranscriptAnnotation - its gene
overlaps only exist as VariantGeneOverlap rows, so that's what comp het has to group genes on.
"""
from django.contrib.auth.models import User
from django.db import connection
from django.test import TestCase, override_settings

from analysis.models import Analysis, QuadNode, SampleNode, TrioNode
from analysis.models.enums import NodeStatus, QuadInheritance, TrioInheritance
from analysis.models.nodes.sources.quad_node import QuadCompHet
from analysis.models.nodes.sources.trio_node import CompHet
from annotation.fake_annotation import get_fake_annotation_version
from annotation.models import AnnotationRun, VariantGeneOverlap
from annotation.models.models import VariantAnnotationVersion, VariantTranscriptAnnotation
from annotation.tests.test_data_fake_genes import create_fake_transcript_version, _create_fake_gene_version, \
    _insert_transcript_data
from genes.models_enums import AnnotationConsortium
from library.utils import sha256sum_str
from snpdb.models import GenomeBuild, Locus, Sequence, Variant
from snpdb.models.models_cohort import CohortGenotype, CohortGenotypeCollection
from snpdb.tests.utils.fake_cohort_data import create_fake_quad, create_fake_trio
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


def _create_second_gene_on_chr21(genome_build, release):
    """ Gene B - downstream of RUNX1 (chr21) so an SV can span both """
    gene_version = _create_fake_gene_version(genome_build, "ENSG00000000001", "GENEB",
                                             AnnotationConsortium.ENSEMBL)
    data = {
        "id": "ENST00000000001.1",
        "hgnc": "99999",
        "biotype": ["protein_coding"],
        "gene_name": "GENEB",
        "genome_builds": {
            genome_build.name: {
                "url": "fake",
                "exons": [[35100000, 35110000, 0, 1, 10000, None]],
                "contig": "21",
                "strand": "+",
                "cds_start": 35100000,
                "cds_end": 35110000,
            }
        },
    }
    return _insert_transcript_data(genome_build, data, gene_version, release)


def _create_test_sv(genome_build, chrom: str, position: int, end: int) -> Variant:
    contig = genome_build.chrom_contig_mappings[chrom]
    ref_seq, _ = Sequence.objects.get_or_create(seq="N", seq_sha256_hash=sha256sum_str("N"))
    alt_seq, _ = Sequence.objects.get_or_create(seq="<DEL>", seq_sha256_hash=sha256sum_str("<DEL>"))
    locus, _ = Locus.objects.get_or_create(contig=contig, position=position, ref=ref_seq)
    variant, _ = Variant.objects.get_or_create(locus=locus, alt=alt_seq, svlen=-(end - position),
                                               defaults={"end": end})
    return variant


def _move_annotation_to_partitions(variant_annotation_version):
    """ Annotation tables are partitioned by inheritance - ORM inserts land in the base table
        while node querysets are rewritten to read from the partition """
    with connection.cursor() as cursor:
        for base_table in (VariantAnnotationVersion.TRANSCRIPT_ANNOTATION,
                           VariantAnnotationVersion.VARIANT_GENE_OVERLAP):
            partition = variant_annotation_version.get_partition_table(base_table_name=base_table)
            cursor.execute(f'INSERT INTO "{partition}" SELECT * FROM ONLY "{base_table}" WHERE version_id = %s',
                           [variant_annotation_version.pk])
            cursor.execute(f'DELETE FROM ONLY "{base_table}" WHERE version_id = %s',
                           [variant_annotation_version.pk])


class AbstractSVCompHetTest:
    """ Gene A (RUNX1) - mother SNV + father SNV, both annotated per-transcript by VEP
        Gene B (GENEB) - mother SNV + father SV, where the SV only has VariantGeneOverlap rows

        MOTHER_SIDE/FATHER_SIDE are samples_zygosity strings for a comp het hit from that parent """
    MOTHER_SIDE = None
    FATHER_SIDE = None

    @classmethod
    def _create_family(cls, user):
        """ Returns (family, cohort_genotype_collection) """
        raise NotImplementedError()

    @classmethod
    def _make_comp_het_node(cls):
        raise NotImplementedError()

    @classmethod
    def _inheritance_handler(cls, node):
        raise NotImplementedError()

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='testuser_sv_comp_het')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.grch37)
        cls.vav = cls.annotation_version.variant_annotation_version
        release = cls.annotation_version.gene_annotation_version.gene_annotation_release

        cls.tv_a = create_fake_transcript_version(cls.grch37, release=release)  # RUNX1 chr21:34.7M-35.0M
        cls.tv_b = _create_second_gene_on_chr21(cls.grch37, release)  # GENEB chr21:35.1M-35.11M
        cls.gene_a = cls.tv_a.gene_version.gene
        cls.gene_b = cls.tv_b.gene_version.gene

        cls.family, cls.cgc = cls._create_family(user)
        cls.analysis = Analysis(genome_build=cls.grch37, annotation_version=cls.annotation_version)
        cls.analysis.set_defaults_and_save(user)
        cls.annotation_run = AnnotationRun.objects.create()

        # Gene A - regular SNV comp het pair
        cls.snv_a_mum = cls._make_snv(34_800_000, cls.MOTHER_SIDE, cls.tv_a)
        cls.snv_a_dad = cls._make_snv(34_810_000, cls.FATHER_SIDE, cls.tv_a)

        # Gene B - mother's SNV
        cls.snv_b_mum = cls._make_snv(35_105_000, cls.MOTHER_SIDE, cls.tv_b)

        # Father's DEL spanning gene A + gene B. VEP skipped it (TOO_LONG) so it has no
        # VariantTranscriptAnnotation - just VariantGeneOverlap rows
        cls.sv = _create_test_sv(cls.grch37, "21", 34_700_000, 35_200_000)
        cls._make_cg(cls.sv, cls.FATHER_SIDE)
        for gene in (cls.gene_a, cls.gene_b):
            VariantGeneOverlap.objects.create(version=cls.vav, annotation_run=cls.annotation_run,
                                              variant=cls.sv, gene=gene)

        _move_annotation_to_partitions(cls.vav)

    @classmethod
    def _make_snv(cls, position: int, samples_zygosity: str, transcript_version) -> Variant:
        variant = slowly_create_test_variant("21", position, "A", "T", cls.grch37)
        cls._make_cg(variant, samples_zygosity)
        gene = transcript_version.gene_version.gene
        VariantTranscriptAnnotation.objects.create(version=cls.vav, variant=variant,
                                                   annotation_run=cls.annotation_run,
                                                   gene=gene, transcript=transcript_version.transcript,
                                                   transcript_version=transcript_version)
        VariantGeneOverlap.objects.create(version=cls.vav, annotation_run=cls.annotation_run,
                                          variant=variant, gene=gene)
        return variant

    @classmethod
    def _make_cg(cls, variant, samples_zygosity):
        n = len(samples_zygosity)
        CohortGenotype.objects.create(
            collection=cls.cgc, variant=variant,
            ref_count=samples_zygosity.count('R'),
            het_count=samples_zygosity.count('E'),
            hom_count=samples_zygosity.count('O'),
            samples_zygosity=samples_zygosity,
            samples_allele_depth=[20] * n, samples_allele_frequency=[100] * n,
            samples_read_depth=[30] * n, samples_genotype_quality=[30] * n,
            samples_phred_likelihood=[0] * n,
        )

    def _comp_het_node(self):
        parent = SampleNode.objects.create(analysis=self.analysis, sample=self.family.proband.sample)
        parent.count = parent.get_queryset().count()
        parent.status = NodeStatus.READY
        parent.save()

        node = self._make_comp_het_node()
        node.add_parent(parent)
        return node.__class__.objects.get(pk=node.pk)

    def _two_hit_genes(self):
        node = self._comp_het_node()
        _, two_hit_genes = self._inheritance_handler(node)._get_comp_het_q_and_two_hit_genes()
        return two_hit_genes

    def _comp_het_variants(self):
        arg_q_dict = self._comp_het_node()._get_node_arg_q_dict()
        qs = Variant.objects.annotate(**self.cgc.get_annotation_kwargs())
        for alias_qs in arg_q_dict.values():
            for q in alias_qs.values():
                qs = qs.filter(q)
        return set(qs.values_list('pk', flat=True))

    def test_two_hit_genes_from_snv_pair(self):
        self.assertIn(self.gene_a.pk, self._two_hit_genes(), "Gene A - mother SNV + father SNV")

    def test_two_hit_genes_from_snv_and_sv(self):
        """ Issue #940 - the SV is gene B's only father-side hit """
        self.assertIn(self.gene_b.pk, self._two_hit_genes(), "Gene B - mother SNV + father SV")

    def test_snv_pair_returned(self):
        variant_ids = self._comp_het_variants()
        self.assertIn(self.snv_a_mum.pk, variant_ids, "Mother's SNV in gene A")
        self.assertIn(self.snv_a_dad.pk, variant_ids, "Father's SNV in gene A")

    def test_sv_returned_for_gene_with_snv_pair(self):
        """ The SV also overlaps gene A, which has 2 hits from the SNV pair """
        self.assertIn(self.sv.pk, self._comp_het_variants(), "Multi-gene SV overlapping gene A")

    def test_snv_and_sv_pair_returned(self):
        variant_ids = self._comp_het_variants()
        self.assertIn(self.snv_b_mum.pk, variant_ids, "Mother's SNV in gene B")
        self.assertIn(self.sv.pk, variant_ids, "Father's SV in gene B")


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestTrioSVCompHet(AbstractSVCompHetTest, TestCase):
    MOTHER_SIDE = "EER"  # proband HET, mother HET, father HOM_REF
    FATHER_SIDE = "ERE"  # proband HET, mother HOM_REF, father HET

    @classmethod
    def _create_family(cls, user):
        trio = create_fake_trio(user, cls.grch37)
        return trio, CohortGenotypeCollection.objects.get(cohort=trio.cohort)

    @classmethod
    def _make_comp_het_node(cls):
        return TrioNode.objects.create(analysis=cls.analysis, trio=cls.family,
                                       inheritance=TrioInheritance.COMPOUND_HET)

    @classmethod
    def _inheritance_handler(cls, node):
        return CompHet(node)


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestQuadSVCompHet(AbstractSVCompHetTest, TestCase):
    # Quad comp het requires the sibling to share each hit
    MOTHER_SIDE = "EERE"  # proband HET, mother HET, father HOM_REF, sibling HET
    FATHER_SIDE = "EREE"  # proband HET, mother HOM_REF, father HET, sibling HET

    @classmethod
    def _create_family(cls, user):
        quad = create_fake_quad(user, cls.grch37, sibling_affected=True)
        return quad, CohortGenotypeCollection.objects.get(cohort=quad.cohort)

    @classmethod
    def _make_comp_het_node(cls):
        return QuadNode.objects.create(analysis=cls.analysis, quad=cls.family,
                                       inheritance=QuadInheritance.COMPOUND_HET)

    @classmethod
    def _inheritance_handler(cls, node):
        return QuadCompHet(node)
