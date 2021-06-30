from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, CohortNode, CohortNodeZygosityFiltersCollection, FilterNode, FilterNodeItem, \
    SampleNode, GeneListNode, IntersectionNode, PhenotypeNode, PopulationNode, TagNode, VariantTag, \
    AlleleFrequencyNode, NodeAlleleFrequencyFilter, NodeVCFFilter, TrioNode
from annotation.fake_annotation import get_fake_annotation_version, create_fake_variants
from annotation.models import VariantGeneOverlap, AnnotationRun
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import GeneList, GeneListGeneSymbol
from ontology.models import OntologyTerm
from patients.models_enums import GnomADPopulation
from snpdb.models import GenomeBuild, ImportStatus, Variant, GenomicInterval, Tag
from snpdb.tests.test_data import create_fake_trio


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestCloneAnalysisNodes(TestCase):
    """ It's only worth testing AnalysisNodes with related objects, as normal fields are just copied OK """

    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version_grch37 = get_fake_annotation_version(cls.grch37)
        gene_annotation_release = cls.annotation_version_grch37.gene_annotation_version.gene_annotation_release
        cls.transcript_version = create_fake_transcript_version(cls.grch37,
                                                                release=gene_annotation_release)
        cls.gene_symbol = cls.transcript_version.gene_version.gene_symbol
        cls.trio = create_fake_trio(user, cls.grch37)

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
        node_q = str(node._get_node_q())
        clone_q = str(clone._get_node_q())
        class_name = node.get_class_name()
        msg = f"Clone of {class_name} has same get_q()"
        print(f"Type '{class_name}' Node: {node_q} <=> {clone_q}")
        self.assertEqual(node_q, clone_q, msg)

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

    def test_clone_cohort_node(self):
        cohort = self.trio.cohort
        cohort_node = CohortNode.objects.create(analysis=self.analysis, cohort=cohort,
                                                accordion_panel=CohortNode.PER_SAMPLE_ZYGOSITY)

        self._set_af_for_node(cohort_node)
        self._set_vcf_filter_for_node(cohort_node, cohort.vcf)

        fc = CohortNodeZygosityFiltersCollection.get_for_node_and_cohort(cohort_node, cohort)
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
        # TODO: This isn't actually working properly
        self._test_clone_node(gene_list_node)

    def test_clone_intersection_node(self):
        # Doesn't have any related objects, but does need to make its own copy of GenomicInterval so test that
        genomic_interval = GenomicInterval.objects.create(chrom="1", start=10000, end=20000)
        intersection_node = IntersectionNode.objects.create(analysis=self.analysis, genomic_interval=genomic_interval)
        # Save the PK as the actual object gets changed (known issue - ok if we don't save the node)
        original_genomic_interval_id = genomic_interval.pk
        clone = intersection_node.save_clone()
        self.assertNotEqual(original_genomic_interval_id, clone.genomic_interval.pk,
                            "Intersection node made copy of genomic interval")

    # TODO: test_clone_pedigree_node - same as sample/trio etc

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
        sample_node = SampleNode.objects.create(analysis=self.analysis, sample=self.sample)
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

