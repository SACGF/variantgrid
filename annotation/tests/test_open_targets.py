from django.test import TestCase
from django.utils import timezone

from annotation.models.models import VariantAnnotation, VariantAnnotationVersion
from genes.models import Gene, GeneAnnotationImport, GeneSymbol, GeneVersion
from genes.models_enums import AnnotationConsortium
from ontology.models import OntologyImport, OntologyService, OntologyTerm
from snpdb.models.models_genome import GenomeBuild


def _variant_annotation(**kwargs) -> VariantAnnotation:
    """ The open_targets_* properties only read fields (and the build, to look up gene symbols) """
    version = VariantAnnotationVersion(genome_build=GenomeBuild.get_name_or_alias("GRCh38"))
    return VariantAnnotation(version=version, **kwargs)


class OpenTargetsRecordsTests(TestCase):
    """ #1822 - the OpenTargets plugin emits every column for every association, '&'-joined, so the
        open_targets_* fields are parallel arrays that have to be zipped back into records. """

    def test_issue_example_qtl_only(self):
        """ The strings from #1822: three tuQTL records for one gene in three tissues """
        va = _variant_annotation(
            open_targets_study_type="tuqtl&tuqtl&tuqtl",
            open_targets_study_id="study1&study2&study3",
            open_targets_variant_id="1_51427_T_G&1_51427_T_G&1_51427_T_G",
            open_targets_gwas_gene_id="NA&NA&NA",
            open_targets_gwas_diseases="NA&NA&NA",
            open_targets_qtl_gene_id="ENSG00000175756&ENSG00000175756&ENSG00000175756",
            open_targets_qtl_biosample="breast_epithelium&gastroesophageal_sphincter&fibroblast",
        )
        self.assertEqual(va.open_targets_gwas_genes, [])
        qtl_genes = va.open_targets_qtl_genes
        self.assertEqual(len(qtl_genes), 1)
        qtl = qtl_genes[0]
        self.assertEqual(qtl["gene_id"], "ENSG00000175756")
        self.assertEqual(qtl["study_type_label"], "tuQTL")
        self.assertEqual(qtl["biosamples"],
                         ["breast epithelium", "fibroblast", "gastroesophageal sphincter"])
        self.assertEqual(qtl["study_count"], 3)
        self.assertNotIn("NA", str(qtl))

    def test_no_open_targets_data(self):
        self.assertEqual(_variant_annotation().open_targets_records, [])
        self.assertEqual(_variant_annotation().open_targets_gwas_genes, [])
        self.assertEqual(_variant_annotation().open_targets_qtl_genes, [])
        self.assertIsNone(_variant_annotation().open_targets_variant)

    def test_mixed_gwas_and_qtl(self):
        """ Two GWAS genes scored differently across studies, plus eQTL/pQTL records for one gene """
        va = _variant_annotation(
            open_targets_study_type="gwas&gwas&gwas&eqtl&pqtl",
            open_targets_study_id="GCST1&GCST1&GCST2&eqtl_study&pqtl_study",
            open_targets_variant_id="&".join(["1_26695422_G_C"] * 5),
            open_targets_gwas_gene_id="ENSG00000117713&ENSG00000000001&ENSG00000117713&NA&NA",
            open_targets_gwas_l2g_scores="0.24&0.9&0.76&NA&NA",
            open_targets_gwas_diseases="EFO_0004527&EFO_0004527&MONDO_0005148&NA&NA",
            open_targets_qtl_gene_id="NA&NA&NA&ENSG00000117713&ENSG00000117713",
            open_targets_qtl_biosample="NA&NA&NA&liver&blood_plasma",
        )
        gwas_genes = va.open_targets_gwas_genes
        self.assertEqual([g["gene_id"] for g in gwas_genes],
                         ["ENSG00000000001", "ENSG00000117713"])  # sorted by L2G score desc
        self.assertEqual(gwas_genes[0]["l2g_score"], 0.9)

        top_gene = gwas_genes[1]
        self.assertEqual(top_gene["l2g_score"], 0.76)  # max across its 2 studies
        self.assertEqual([d["id"] for d in top_gene["diseases"]], ["EFO_0004527", "MONDO_0005148"])
        self.assertEqual(top_gene["study_count"], 2)

        qtl_genes = va.open_targets_qtl_genes
        self.assertEqual([(g["gene_id"], g["study_type"]) for g in qtl_genes],
                         [("ENSG00000117713", "eqtl"), ("ENSG00000117713", "pqtl")])
        self.assertEqual(qtl_genes[0]["biosamples"], ["liver"])
        self.assertEqual(qtl_genes[1]["biosamples"], ["blood plasma"])

    def test_l2g_scores_null_before_backfill(self):
        """ open_targets_gwas_l2g_scores is added by #1822 - old rows are null until backfilled """
        va = _variant_annotation(
            open_targets_study_type="gwas&gwas",
            open_targets_study_id="GCST1&GCST2",
            open_targets_gwas_gene_id="ENSG00000117713&ENSG00000117713",
            open_targets_gwas_diseases="EFO_0004527&EFO_0004527",
            open_targets_gwas_l2g_score=0.76,
        )
        gwas_genes = va.open_targets_gwas_genes
        self.assertEqual(len(gwas_genes), 1)
        self.assertIsNone(gwas_genes[0]["l2g_score"])
        self.assertEqual(gwas_genes[0]["study_count"], 2)

    def test_variant_link_uses_a_single_id(self):
        va = _variant_annotation(
            open_targets_study_type="gwas&gwas",
            open_targets_variant_id="1_51427_T_G&1_51427_T_G",
        )
        self.assertEqual(va.open_targets_variant,
                         {"id": "1_51427_T_G",
                          "url": "https://platform.opentargets.org/variant/1_51427_T_G"})

    def test_biosample_comma_does_not_break_alignment(self):
        """ VEP escapes ',' to '&' (and the space after it to '_'), so a biosample name containing a
            comma splits into an extra element that has to be rejoined """
        va = _variant_annotation(
            open_targets_study_type="eqtl&eqtl",
            open_targets_study_id="study1&study2",
            open_targets_qtl_gene_id="ENSG00000175756&ENSG00000175756",
            open_targets_qtl_biosample="CD4-positive&_alpha-beta_T_cell&liver",
        )
        qtl_genes = va.open_targets_qtl_genes
        self.assertEqual(len(qtl_genes), 1)
        self.assertEqual(qtl_genes[0]["biosamples"], ["CD4-positive, alpha-beta T cell", "liver"])
        self.assertEqual(qtl_genes[0]["study_count"], 2)

    def test_multi_disease_record_assigned_within_its_run(self):
        """ Diseases are '|'-separated within a record, which VEP also escapes to '&', so a run of
            consecutive GWAS records can carry more diseases than it has records """
        va = _variant_annotation(
            open_targets_study_type="gwas&eqtl&gwas",
            open_targets_study_id="GCST1&eqtl_study&GCST2",
            open_targets_gwas_gene_id="ENSG00000117713&NA&ENSG00000000001",
            open_targets_gwas_diseases="OBA_2050103&EFO_0007702&NA&EFO_0004527",
            open_targets_qtl_gene_id="NA&ENSG00000117713&NA",
        )
        diseases_by_gene = {g["gene_id"]: [d["id"] for d in g["diseases"]]
                            for g in va.open_targets_gwas_genes}
        self.assertEqual(diseases_by_gene["ENSG00000117713"], ["EFO_0007702", "OBA_2050103"])
        self.assertEqual(diseases_by_gene["ENSG00000000001"], ["EFO_0004527"])


class OpenTargetsDiseaseNameTests(TestCase):
    """ Disease ids are ontology CURIEs with '_' - MONDO/HP we store locally so can name """

    def setUp(self):
        o_import = OntologyImport.objects.create(import_source=OntologyService.MONDO, filename="test",
                                                 context="test", hash="N/A", processor_version=1,
                                                 completed=True, processed_date=timezone.now())
        OntologyTerm.objects.create(id="MONDO:0005148", ontology_service=OntologyService.MONDO,
                                    index=5148, name="type 2 diabetes mellitus", from_import=o_import)

    def test_mondo_named_others_show_id(self):
        va = _variant_annotation(
            open_targets_study_type="gwas&gwas&gwas",
            open_targets_study_id="GCST1&GCST2&GCST3",
            open_targets_gwas_gene_id="&".join(["ENSG00000117713"] * 3),
            open_targets_gwas_diseases="MONDO_0005148&EFO_0004527&MONDO_9999999",
        )
        diseases = va.open_targets_gwas_genes[0]["diseases"]
        self.assertEqual([d["name"] for d in diseases],
                         ["EFO_0004527", "type 2 diabetes mellitus", "MONDO_9999999"])
        self.assertEqual(diseases[1]["url"], "https://platform.opentargets.org/disease/MONDO_0005148")


class OpenTargetsGeneLinkTests(TestCase):
    """ Open Targets uses Ensembl gene ids, which is what Gene.pk holds - so genes we've imported get
        their symbol and a link to our gene page, and the rest link out to Ensembl """

    def setUp(self):
        genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        gene = Gene.objects.create(identifier="ENSG00000117713",
                                   annotation_consortium=AnnotationConsortium.ENSEMBL)
        import_source = GeneAnnotationImport.objects.create(
            url="fake", genome_build=genome_build, annotation_consortium=AnnotationConsortium.ENSEMBL)
        for version, symbol in [(1, "OLD_SYMBOL"), (2, "ARID1A")]:
            GeneVersion.objects.create(gene=gene, version=version, genome_build=genome_build,
                                       gene_symbol=GeneSymbol.objects.get_or_create(symbol=symbol)[0],
                                       import_source=import_source)

    def test_known_gene_uses_latest_symbol_and_gene_page(self):
        va = _variant_annotation(
            open_targets_study_type="gwas&gwas",
            open_targets_study_id="GCST1&GCST2",
            open_targets_gwas_gene_id="ENSG00000117713&ENSG00000999999",
            open_targets_gwas_diseases="EFO_0004527&EFO_0004527",
        )
        by_id = {g["gene_id"]: g for g in va.open_targets_gwas_genes}
        self.assertEqual(by_id["ENSG00000117713"]["gene_symbol"], "ARID1A")
        self.assertEqual(by_id["ENSG00000117713"]["gene_url"], "/genes/view_gene/ENSG00000117713")
        self.assertEqual(by_id["ENSG00000999999"]["gene_symbol"], "ENSG00000999999")
        self.assertEqual(by_id["ENSG00000999999"]["gene_url"],
                         "https://ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000999999")
