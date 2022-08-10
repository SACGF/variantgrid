import logging

import pandas as pd
from django.core.management import BaseCommand, CommandError

from annotation.models import DBNSFPGeneAnnotationVersion, DBNSFPGeneAnnotation
from genes.models import GeneSymbol
from library.file_utils import file_md5sum
from library.pandas_utils import df_nan_to_none


class Command(BaseCommand):
    COLUMNS_TO_FIELDS = {
        "Essential_gene_CRISPR": "essential_gene_crispr",
        "Essential_gene_CRISPR2": "essential_gene_crispr2",
        "Essential_gene_gene-trap": "essential_gene_gene_trap",
        "Expression(GNF/Atlas)": "expression_gnf_atlas",
        "Expression(egenetics)": "expression_egenetics",
        "GDI": "gene_damage_index_score",
        "GDI-Phred": "gene_damage_index_phred",
        "GHIS": "ghis",
        "GO_biological_process": "go_biological_process",
        "GO_cellular_component": "go_cellular_component",
        "GO_molecular_function": "go_molecular_function",
        "Gene_indispensability_pred": "gene_indispensability_pred",
        "Gene_indispensability_score": "gene_indispensability_score",
        "HIPred": "hipred_prediction",
        "HIPred_score": "hipred_score",
        "Interactions(BioGRID)": "interactions_biogrid",
        "Interactions(ConsensusPathDB)": "interactions_consensus_pathdb",
        "LoFtool_score": "loftool",
        "P(HI)": "phi",
        "P(rec)": "prec",
        "Pathway(BioCarta)_full": "pathway_biocarta_full",
        "Pathway(ConsensusPathDB)": "pathway_consensus_pathdb",
        "Pathway(KEGG)_full": "pathway_kegg_full",
        "Pathway(KEGG)_id": "pathway_kegg_id",
        "Trait_association(GWAS)": "gwas_trait_association",
        "gnomAD_pLI": "gnomad_pli",
        "gnomAD_pNull": "gnomad_pnull",
        "gnomAD_pRec": "gnomad_prec",
    }

    def add_arguments(self, parser):
        parser.add_argument('--replace', action='store_true', help="Replace existing version if exists")
        parser.add_argument('--dbnsfp-version', required=True, help="dbNSFP version eg '4.3'")
        parser.add_argument('dbnsfp_gene_filename')

    def handle(self, *args, **options):
        replace = options["replace"]
        version = options["dbnsfp_version"]
        filename = options["dbnsfp_gene_filename"]
        if not filename.endswith(".gz"):
            raise CommandError("Unknown file type, expecting .gz")

        md5_hash = file_md5sum(filename)
        dbnsfp_version, created = DBNSFPGeneAnnotationVersion.objects.get_or_create(version=version,
                                                                                    md5_hash=md5_hash)
        if not created:
            if replace:
                logging.warning("--replace option used, deleting existing objects")
                if dbnsfp_version.geneannotationversion_set.exists():
                    logging.warning("Gene Annotation relies on this - you will need to re-run:")
                    logging.warning("python3 manage.py gene_annotation --add-dbnsfp-gene")
                dbnsfp_version.dbnsfpgeneannotation_set.all().delete()
            else:
                msg = f"Existing DBNSFPGeneAnnotationVersion = {version}, use --replace to delete existing objects"
                raise ValueError(msg)

        self._import_dbnsfp_gene(filename, dbnsfp_version)

    @staticmethod
    def format_field(field, value):
        LOOKUPS = {
            "gene_indispensability_pred": {"E": True, "N": False},
            "hipred_prediction": {"Y": True, "N": False},
        }
        if value is not None:
            if lookup := LOOKUPS.get(field):
                value = lookup[value]
        return value

    def _import_dbnsfp_gene(self, filename: str, dbnsfp_version):
        records = []
        df = pd.read_csv(filename, sep='\t', index_col=None)
        df = df_nan_to_none(df.replace(".", None))
        gene_symbols = GeneSymbol.get_upper_case_lookup()
        new_gene_symbols = set()

        for _, row in df.iterrows():
            gene_symbol_str = row["Gene_name"]
            gene_symbol_id = gene_symbols.get(gene_symbol_str.upper())
            if gene_symbol_id is None:
                new_gene_symbols.add(GeneSymbol(symbol=gene_symbol_str))
                gene_symbol_id = gene_symbol_str

            kwargs = {
                "version": dbnsfp_version,
                "gene_symbol_id": gene_symbol_id,
            }
            for column, field in self.COLUMNS_TO_FIELDS.items():
                value = row[column]
                kwargs[field] = self.format_field(field, value)

            records.append(DBNSFPGeneAnnotation(**kwargs))

        if new_gene_symbols:
            GeneSymbol.objects.bulk_create(new_gene_symbols, batch_size=2000)

        if records:
            logging.info("Creating %d DBNSFPGeneAnnotation records", len(records))
            DBNSFPGeneAnnotation.objects.bulk_create(records, batch_size=2000)
