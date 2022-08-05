import os

from django.conf import settings
from django.core.management import BaseCommand

from annotation.models import DBNSFPGeneAnnotationVersion, DBNSFPGeneAnnotation


class Command(BaseCommand):
    COLUMNS_TO_FIELDS = {
        "Pathway(BioCarta)_full": "pathway_biocarta_full",
        "Pathway(ConsensusPathDB)": "pathway_consensus_pathdb",
        "Pathway(KEGG)_id": "pathway_kegg_id",
        "Pathway(KEGG)_full": "pathway_kegg_full",
        "Trait_association(GWAS)": "gwas_trait_association",
        "GO_biological_process": "go_biological_process",
        "GO_cellular_component": "go_cellular_component",
        "GO_molecular_function": "go_molecular_function",
        "Interactions(BioGRID)": "interactions_biogrid",
        "Interactions(ConsensusPathDB)": "interactions_consensus_pathdb",
        "gnomAD_pLI": "gnomad_pli",
        "gnomAD_pRec": "gnomad_prec",
        "gnomAD_pNull": "gnomad_pnull",
        "LoFtool_score": "loftool",
        "Essential_gene_CRISPR": "essential_gene_crispr",
        "Essential_gene_CRISPR2": "essential_gene_crispr2",
        "Essential_gene_gene-trap": "essential_gene_gene_trap",
    }

    # This is for the eventual variant grid columns data migration
    __DESCRIPTIONS = {
        "Pathway(BioCarta)_short": "Short name of the Pathway(s) the gene belongs to (from BioCarta)",
        "Pathway(BioCarta)_full": "Full name(s) of the Pathway(s) the gene belongs to (from BioCarta)",
        "Pathway(ConsensusPathDB)": "Pathway(s) the gene belongs to (from ConsensusPathDB)",
        "Pathway(KEGG)_id": "ID(s) of the Pathway(s) the gene belongs to (from KEGG)",
        "Pathway(KEGG)_full": "Full name(s) of the Pathway(s) the gene belongs to (from KEGG)",
        "Trait_association(GWAS)": "Trait(s) the gene associated with (from GWAS catalog)",
        "GO_biological_process": "GO terms for biological process",
        "GO_cellular_component": "GO terms for cellular component",
        "GO_molecular_function": "GO terms for molecular function",
        "Interactions(BioGRID)": "The number of other genes this gene interacting with (from BioGRID)",
        "Interactions(ConsensusPathDB)": "The number of other genes this gene interacting with (from ConsensusPathDB).",
        "gnomAD_pLI": '"The probability of being loss-of-function intolerant (intolerant of both heterozygous and homozygous lof variants)" based on gnomAD 2.1 data',
        "gnomAD_pRec": '"the probability of being intolerant of homozygous, but not heterozygous lof variants" based on gnomAD 2.1 data',
        "gnomAD_pNull": '"the probability of being tolerant of both heterozygous and homozygous lof variants" based on gnomAD 2.1 data',
        "LoFtool_score": "a percentile score for gene intolerance to functional change. The lower the score the higher gene intolerance to functional change. For details see doi: 10.1093/bioinformatics/btv602.",
        "Essential_gene_CRISPR": 'Essential ("E") or Non-essential phenotype-changing ("N") based on large scale CRISPR experiments. from doi: 10.1126/science.aac7041',
        "Essential_gene_CRISPR2": 'Essential ("E"), context-Specific essential ("S"), or Non-essential phenotype-changing ("N") based on large scale CRISPR experiments. from http://dx.doi.org/10.1016/j.cell.2015.11.015',
        "Essential_gene_gene-trap": 'Essential ("E"), HAP1-Specific essential ("H"), KBM7-Specific essential ("K"), or Non-essential phenotype-changing ("N"), based on large scale mutagenesis experiments. from doi: 10.1126/science.aac7557',
    }

    def add_arguments(self, parser):
        parser.add_argument('--force', action="store_true", help="Force create new GeneAnnotation for same release")

    def handle(self, *args, **options):
        pass