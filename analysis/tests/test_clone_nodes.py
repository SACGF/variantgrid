from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import (
    AlleleFrequencyNode,
    Analysis,
    CohortNode,
    CohortNodeZygosityFiltersCollection,
    FilterNode,
    FilterNodeItem,
    GeneListNode,
    IntersectionNode,
    MOINode,
    NodeAlleleFrequencyFilter,
    NodeVCFFilter,
    PedigreeNode,
    PhenotypeNode,
    PopulationNode,
    SampleNode,
    TagNode,
    TrioNode,
    VariantTag,
)
from annotation.fake_annotation import create_fake_variants, get_fake_annotation_version
from annotation.models import AnnotationRun, VariantGeneOverlap
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import GeneList, GeneListGeneSymbol
from ontology.models import OntologyTerm
from patients.models_enums import GnomADPopulation
from pedigree.models import PedigreeInheritance
from snpdb.models import GenomeBuild, ImportStatus, Tag, Variant
from snpdb.tests.utils.fake_cohort_data import create_fake_pedigree, create_fake_trio


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestCloneAnalysisNodes(TestCase):
    """ It's only worth testing AnalysisNodes with related objects, as normal fields are just copied OK """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version_grch37 = get_fake_annotation_version(cls.grch37)
        gene_annotation_release = cls.annotation_version_grch37.gene_annotation_version.gene_annotation_release
        cls.transcript_version = create_fake_transcript_version(cls.grch37,
                                                                release=gene_annotation_release)
        cls.gene_symbol = cls.transcript_version.gene_version.gene_symbol
        cls.trio = create_fake_trio(user, cls.grch37)
        cls.pedigree = create_fake_pedigree(user, cls.grch37)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

        cls.sample = cls.trio.get_samples()[0]

        # Gene List
        cls.gene_list = GeneList.objects.get_or_create(name="fake list",
                                                       user=cls.analysis.user,
                                                       import_status=ImportStatus.SUCCESS)[0]
        GeneListGeneSymbol.objects.get_or_create(gene_list=cls.gene_list, gene_symbol=cls.gene_symbol)

        # Need some overlapping variants so gene list will work
        create_fake_variants(cls.grch37)
        # Note: Variant probably doesn't overlap with gene, just want a random one
        variant = Variant.objects.filter(Variant.get_no_reference_q()).first()
        annotation_run = AnnotationRun.objects.create()
        VariantGeneOverlap.objects.create(version=cls.annotation_version_grch37.variant_annotation_version,
                                          annotation_run=annotation_run,
                                          gene=cls.transcript_version.gene,
                                          variant=variant)

        # Tag that variant
        cls.tag = Tag.objects.get_or_create(pk="foo")[0]
        VariantTag.objects.create(genome_build=cls.grch37, analysis=cls.analysis,
                                  variant=variant, tag=cls.tag, user=user)

    def _test_clone_node(self, node):
        clone = node.save_clone()
        node_arg_q_dict = str(node._get_node_arg_q_dict())
        clone_arg_q_dict = str(clone._get_node_arg_q_dict())
        class_name = node.get_class_name()
        msg = f"Clone of {class_name} has same get_q()"
        self.assertEqual(node_arg_q_dict, clone_arg_q_dict, msg)

    @staticmethod
    def _set_af_for_node(node):
        naff = NodeAlleleFrequencyFilter.objects.get_or_create(node=node)[0]
        naff.nodeallelefrequencyrange_set.get_or_create(min=15, max=100)

    @staticmethod
    def _set_vcf_filter_for_node(node, vcf):
        vcf_filter = vcf.vcffilter_set.first()
        NodeVCFFilter.objects.create(node=node, vcf_filter=vcf_filter)

    def test_clone_allele_frequency_node(self):
        af_node = AlleleFrequencyNode.objects.create(analysis=self.analysis, sample=self.sample)
        self._set_af_for_node(af_node)
        self._test_clone_node(af_node)

    def _get_cohort_node(self):
        cohort = self.trio.cohort
        cohort_node = CohortNode.objects.create(analysis=self.analysis, cohort=cohort,
                                                accordion_panel=CohortNode.PER_SAMPLE_ZYGOSITY)

        self._set_af_for_node(cohort_node)
        self._set_vcf_filter_for_node(cohort_node, cohort.vcf)
        return cohort_node

    def _get_sample_node(self):
        return SampleNode.objects.create(analysis=self.analysis, sample=self.sample)

    def test_clone_cohort_node(self):
        cohort_node = self._get_cohort_node()
        fc = CohortNodeZygosityFiltersCollection.get_for_node_and_cohort(cohort_node, cohort_node.cohort)
        zf1 = fc.cohortnodezygosityfilter_set.order_by("pk").first()
        zf1.zygosity_hom = False
        zf1.save()

        zf2 = fc.cohortnodezygosityfilter_set.order_by("pk").last()
        zf2.zygosity_het = False
        zf2.save()

        self._test_clone_node(cohort_node)

    def test_clone_filter_node(self):
        filter_node = FilterNode.objects.create(analysis=self.analysis)
        FilterNodeItem.objects.create(filter_node=filter_node, sort_order=1,
                                      operation="eq", field="id", data="42")

        self._test_clone_node(filter_node)

    def test_clone_gene_list_node(self):
        gene_list_node = GeneListNode.objects.create(analysis=self.analysis)
        gene_list_node.genelistnodegenelist_set.create(gene_list_node=gene_list_node, gene_list=self.gene_list)
        clone = gene_list_node.save_clone()
        self.assertEqual([gln_gl.gene_list for gln_gl in clone.genelistnodegenelist_set.all()],
                         [self.gene_list], "Clone gets its own link to the same GeneList")

    def test_clone_intersection_node(self):
        # Resolved entries are concrete fields, so the clone carries them without save_clone() plumbing
        intersection_node = IntersectionNode.objects.create(analysis=self.analysis,
                                                            accordion_panel=IntersectionNode.VARIANTS,
                                                            variant_text="1:10000-20000",
                                                            variant_regions=["1:10000-20000"])
        clone = intersection_node.save_clone()
        self.assertEqual(clone.variant_text, "1:10000-20000")
        self.assertEqual(clone.variant_regions, ["1:10000-20000"])

    def test_clone_intersection_node_contigs(self):
        intersection_node = IntersectionNode.objects.create(analysis=self.analysis,
                                                            accordion_panel=IntersectionNode.CONTIG)
        contig = self.analysis.genome_build.contigs.first()
        intersection_node.intersectionnodecontig_set.create(contig=contig)
        intersection_node = IntersectionNode.objects.get(pk=intersection_node.pk)  # contig_ids is cached
        clone = intersection_node.save_clone()
        self.assertEqual([inc.contig for inc in clone.intersectionnodecontig_set.all()], [contig])

    def test_clone_moi_node(self):
        # TODO: Need to also test patient setup w/ontology etc
        moi_node = MOINode.objects.create(analysis=self.analysis, sample=self.sample)

        ontology_term = OntologyTerm.objects.first()
        moi_node.moinodeontologyterm_set.create(ontology_term=ontology_term)
        self._test_clone_node(moi_node)

    def test_clone_pedigree_node(self):
        pedigree_node = PedigreeNode.objects.create(analysis=self.analysis, pedigree=self.pedigree,
                                                    inheritance_model=PedigreeInheritance.AUTOSOMAL_DOMINANT)
        self._set_af_for_node(pedigree_node)
        self._set_vcf_filter_for_node(pedigree_node, self.trio.cohort.vcf)
        self._test_clone_node(pedigree_node)

    def test_clone_phenotype_node(self):
        # TODO: Need to insert genes against the ontology term so that it actually does filtering
        phenotype_node = PhenotypeNode.objects.create(analysis=self.analysis)
        ontology_term = OntologyTerm.objects.first()
        phenotype_node.phenotypenodeontologyterm_set.create(ontology_term=ontology_term)
        self._test_clone_node(phenotype_node)

    def test_clone_population_node(self):
        population_node = PopulationNode.objects.create(analysis=self.analysis)
        population_node.populationnodegnomadpopulation_set.create(population=GnomADPopulation.NON_FINNISH_EUROPEAN)
        self._test_clone_node(population_node)

    def test_clone_sample_node(self):
        # Test allele frequency, VCF filters
        sample_node = self._get_sample_node()
        self._set_af_for_node(sample_node)
        self._set_vcf_filter_for_node(sample_node, self.sample.vcf)
        self._test_clone_node(sample_node)

    def test_clone_tag_node(self):
        tag_node = TagNode.objects.create(analysis=self.analysis)
        tag_node.tagnodetag_set.create(tag=self.tag)
        self._test_clone_node(tag_node)

    def test_clone_trio_node(self):
        trio_node = TrioNode.objects.create(analysis=self.analysis, trio=self.trio)
        self._set_af_for_node(trio_node)
        self._set_vcf_filter_for_node(trio_node, self.trio.cohort.vcf)
        self._test_clone_node(trio_node)
