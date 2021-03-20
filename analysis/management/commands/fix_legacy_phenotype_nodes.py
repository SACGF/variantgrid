""" Old PhenotypeNode ontology records needed to be saved, and restored after new ontology is loaded.
    Files may have been written in 0014_one_off_move_ontology.py
    This can be deleted once all environments have been upgraded """
from pathlib import Path

from django.conf import settings
from django.core.management import BaseCommand
from django.db import IntegrityError
import pandas as pd
from analysis.models import PhenotypeNodeOntologyTerm
from ontology.models import OntologyTerm, OntologyService


class Command(BaseCommand):
    def handle(self, *args, **options):
        omim = OntologyTerm.objects.filter(ontology_service=OntologyService.OMIM)
        if not omim.exists():
            raise ValueError("No OMIM ontology terms exist. Run 'ontology_import' first")

        hpo = OntologyTerm.objects.filter(ontology_service=OntologyService.HPO)
        if not hpo.exists():
            raise ValueError("No HPO ontology terms exist. Run 'ontology_import' first")

        # Some OMIM terms became obsolete / moved etc.
        OMIM_CHANGES = {
             614087: [227650, 609644]
        }

        migrations_dir = Path(settings.BASE_DIR) / "data/migrations"
        omim_csv = migrations_dir / "omim_legacy.csv"
        if omim_csv.exists():
            df = pd.read_csv(omim_csv)
            for _, row in df.iterrows():
                node_id = row["phenotype_node"]
                omim_id = row["mim_morbid_alias__mim_morbid"]
                omim_ids = []
                if modified_omim := OMIM_CHANGES.get(omim_id):
                    print(f"{omim_id} => {modified_omim}")
                    omim_ids = modified_omim
                else:
                    omim_ids = [omim_id]

                for omim_id in omim_ids:
                    term = f"OMIM:{omim_id}"
                    PhenotypeNodeOntologyTerm.objects.get_or_create(phenotype_node_id=node_id, ontology_term_id=term)
        else:
            print(f"{omim_csv} not found")

        hpo_csv = migrations_dir / "hpo_legacy.csv"
        if hpo_csv.exists():
            hpo_df = pd.read_csv(hpo_csv)
            for _, row in hpo_df.iterrows():
                node_id = row["phenotype_node"]
                hpo_id = row["hpo_synonym__hpo"]

                term = "HP:%07d" % hpo_id
                try:
                    PhenotypeNodeOntologyTerm.objects.get_or_create(phenotype_node_id=node_id, ontology_term_id=term)
                except IntegrityError as e:
                    print(e)
        else:
            print(f"{hpo_csv} not found")
