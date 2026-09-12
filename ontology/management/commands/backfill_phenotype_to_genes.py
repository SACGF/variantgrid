"""
Backfill HPO -> OMIM -> gene relations into the latest OntologyVersion's phenotype_to_genes import (#1862).

HPO changed phenotype_to_genes.txt to a five column format in 2023 and the loader silently imported nothing from
it, so phenotype nodes found no genes for any HPO term. Rather than importing afresh (a new OntologyImport means a
new OntologyVersion, annotation sub-version and gene annotation build), this writes the relations into the import
the latest version already points at. Downloads the current file from HPO when --phenotype_to_genes is not given.
"""
import tempfile
from pathlib import Path

import requests
from django.core.management import BaseCommand

from library.constants import MINUTE_SECS
from ontology.management.commands.ontology_import import load_phenotype_to_genes
from ontology.models import OntologyTermRelation, OntologyVersion

PHENOTYPE_TO_GENES_URL = "http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt"


def phenotype_to_genes_import_is_empty(ontology_version: OntologyVersion) -> bool:
    return not OntologyTermRelation.objects.filter(from_import_id=ontology_version.hp_phenotype_to_genes_import_id).exists()


def download_phenotype_to_genes(directory: str) -> str:
    filename = Path(directory) / Path(PHENOTYPE_TO_GENES_URL).name
    print(f"Downloading {PHENOTYPE_TO_GENES_URL}")
    with requests.get(PHENOTYPE_TO_GENES_URL, stream=True, timeout=5 * MINUTE_SECS) as r:
        r.raise_for_status()
        with open(filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                f.write(chunk)
    return str(filename)


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument('--phenotype_to_genes', required=False,
                            help="phenotype_to_genes.txt to load (downloaded from HPO if not given)")
        parser.add_argument('--force', action="store_true",
                            help="Backfill even if the latest version's import already has relations")

    def handle(self, *args, **options):
        ontology_version = OntologyVersion.objects.order_by("pk").last()
        if ontology_version is None:
            print("No OntologyVersion - import the ontology first (see annotation page)")
            return

        ontology_import = ontology_version.hp_phenotype_to_genes_import
        if not (options["force"] or phenotype_to_genes_import_is_empty(ontology_version)):
            print(f"{ontology_version} {ontology_import} already has relations - nothing to do (use --force to reload)")
            return

        with tempfile.TemporaryDirectory() as tmp_dir:
            filename = options["phenotype_to_genes"] or download_phenotype_to_genes(tmp_dir)
            print(f"Backfilling {ontology_import} from {filename}")
            load_phenotype_to_genes(filename, force=True, existing_import=ontology_import)

        num_relations = OntologyTermRelation.objects.filter(from_import=ontology_import).count()
        print(f"{ontology_import} now has {num_relations} relations.")
        print("Phenotype nodes cache 'no genes' lookups per term for a day (cached_gene_symbols_for_terms_tuple), "
              "so a term viewed before the backfill may keep its warning until then.")
        print("Gene annotation HPO/OMIM columns for this ontology version were built without them; rebuild with:")
        print("python3 manage.py gene_annotation --latest-releases --force")
