"""
Imports OMIM gene level annotation data.

This creates MIMGene links for Ensembl genes. RefSeq genes are done via
OMIM_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt

Source: http://www.ensembl.org/biomart/
Export columns: "MIM morbid description" "MIM morbid accession" "Gene stable ID"
"""

from django.core.management.base import BaseCommand

from annotation.models import MIMMorbid, MIMMorbidAlias, MIMGene, HumanPhenotypeOntology
from annotation.models.models_phenotype_match import TextPhenotypeMatch
from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from library.log_utils import console_logger
import pandas as pd


DO_INSERT = True
BULK_INSERT_SIZE = 10000


def load_diseases_to_genes_to_phenotypes(filename):
    df = pd.read_csv(filename, index_col=None, comment='#', sep='\t',
                     names=['disease_id', 'gene_symbol', 'gene_id', 'hpo_id', 'hpo_name'],
                     dtype={"gene_id": int})
    df["hpo_id"] = df["hpo_id"].str.replace(HumanPhenotypeOntology.PREFIX, "").astype(int)  # Strip HP: from column
    df["disease_id"] = df["disease_id"].str.replace(MIMMorbid.PREFIX, "").astype(int)  # Strip OMIM: from column
    return df


class Command(BaseCommand):
    MIM_DESCRIPTION = "MIM morbid description"
    MIM_ACCESSION = "MIM morbid accession"
    ENSEMBL_GENE_ID = "Gene stable ID"

    def add_arguments(self, parser):
        parser.add_argument('--disease_gene_txt', required=True)
        parser.add_argument('--biomart', required=True)

    @staticmethod
    def insert_description(logger, description_series):
        mim_morbid_list = []
        mim_alias_list = []

        for mim_accession_id, description in description_series.items():
            descriptions_list = [x.strip() for x in str(description).split(";;")]

            description = descriptions_list[0]
            mm = MIMMorbid(accession=mim_accession_id, description=description)
            mim_morbid_list.append(mm)

            # Because we auto-complete MIMMorbidAlias (to get all the terms)
            # we need to create a duplicate of MIMMorbid as an alias
            mma = MIMMorbidAlias(mim_morbid_id=mim_accession_id, description=description)
            mim_alias_list.append(mma)

            for alias_description in descriptions_list[1:]:
                mma = MIMMorbidAlias(mim_morbid_id=mim_accession_id, description=alias_description)
                mim_alias_list.append(mma)

        logger.info(f"Inserting MIMMorbid: {len(mim_morbid_list)}")
        logger.info(f"Inserting Aliases: {len(mim_alias_list)}")

        MIMMorbid.objects.bulk_create(mim_morbid_list, ignore_conflicts=True)
        MIMMorbidAlias.objects.bulk_create(mim_alias_list, ignore_conflicts=True)

    @staticmethod
    def insert_mim_ensembl_genes(logging, mim_ids, known_ensembl_gene_ids, gene_series):
        mim_genes = []
        unknown_genes = []
        for mim_morbid_id, ensembl_gene_id in gene_series.iteritems():
            if ensembl_gene_id in known_ensembl_gene_ids:
                mg = MIMGene(mim_morbid_id=mim_morbid_id, gene_id=ensembl_gene_id)
                mim_genes.append(mg)
            else:
                unknown_genes.append(ensembl_gene_id)
        logging.info(f"Inserting {len(mim_genes)} Ensembl MIMGenes")
        MIMGene.objects.bulk_create(mim_genes, ignore_conflicts=True)

        if unknown_genes:
            logging.warning(f"Could't insert {len(unknown_genes)} unknown Ensembl genes")

    @staticmethod
    def insert_mim_refseq_genes(logging, mim_ids, known_refseq_gene_ids, df):
        disease_gb = df.groupby(['disease_id'])
        mim_gene_records = []
        unknown_genes = set()

        for disease_id, data in disease_gb:
            if int(disease_id) not in mim_ids:
                continue

            for gene_id in data['gene_id']:
                if str(gene_id) in known_refseq_gene_ids:  # geneIds are hashed as strings
                    mg = MIMGene(mim_morbid_id=disease_id, gene_id=gene_id)
                    mim_gene_records.append(mg)
                else:
                    unknown_genes.add(gene_id)

        if DO_INSERT and mim_gene_records:  # Leftovers
            MIMGene.objects.bulk_create(mim_gene_records, ignore_conflicts=True)
            logging.info(f"Inserted {len(mim_gene_records)}")

        if unknown_genes:
            logging.warning(f"Could't insert {len(unknown_genes)} unknown RefSeq genes")

    def handle(self, *args, **options):
        mim_biomart_filename = options['biomart']
        omim_disease_filename = options['disease_gene_txt']

        logger = console_logger()
        logger.info("Deleting existing MIM entries...")
        if TextPhenotypeMatch.objects.exists():
            logger.warning("You have Phenotype matches which will be deleted w/OMIM. Please re-match with 'python3 manage.py match_patient_phenotypes --clear'")

        MIMMorbid.objects.all().delete()  # Will get aliases too

        # Create MIMMorbid from BioMart file
        mim_biomart_df = pd.read_csv(mim_biomart_filename, sep='\t').dropna().astype({"MIM morbid accession": int})
        for expected_col in [self.MIM_DESCRIPTION, self.MIM_ACCESSION, self.ENSEMBL_GENE_ID]:
            if expected_col not in mim_biomart_df.columns:
                msg = f"file {mim_biomart_filename} missing column: '{expected_col}': columns: '{mim_biomart_df.columns}'"
                raise ValueError(msg)

        mim_biomart_df = mim_biomart_df.set_index(self.MIM_ACCESSION)
        description_series = mim_biomart_df[self.MIM_DESCRIPTION]
        Command.insert_description(logger, description_series)

        mim_ids = set(MIMMorbid.objects.all().values_list("pk", flat=True))

        logger.info("Ensembl MIMGene")
        known_ensembl_gene_ids = Gene.known_gene_ids(AnnotationConsortium.ENSEMBL)
        Command.insert_mim_ensembl_genes(logger, mim_ids, known_ensembl_gene_ids, mim_biomart_df[self.ENSEMBL_GENE_ID])

        logger.info("RefSeq MIMGene")
        known_entrez_gene_ids = Gene.known_gene_ids(AnnotationConsortium.REFSEQ)
        mim_gene_df = load_diseases_to_genes_to_phenotypes(omim_disease_filename)
        Command.insert_mim_refseq_genes(logger, mim_ids, known_entrez_gene_ids, mim_gene_df)
