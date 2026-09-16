import warnings
from pathlib import Path

from django.conf import settings
from django.utils import timezone

from genes.models import GeneSymbol
from ontology.management.commands import ontology_import
from ontology.models import OntologyImport, OntologyService, OntologyTerm, OntologyVersion
from ontology.ontology_builder import OntologyBuilderDataUpToDateException

# HGNC terms so phenotype_to_genes.txt can link OMIM -> gene (the other symbols in the file stay unresolved)
TEST_HGNC_GENE_SYMBOLS = {986: "BCKDHA", 2548: "CUBN", 7010: "MEN1"}


def _create_hgnc_test_data():
    ontology_import, _ = OntologyImport.objects.get_or_create(import_source=OntologyService.HGNC, filename="test_hgnc",
                                                              context="test", defaults={"processed_date": timezone.now(),
                                                                                        "completed": True})
    for hgnc_id, symbol in TEST_HGNC_GENE_SYMBOLS.items():
        GeneSymbol.objects.get_or_create(symbol=symbol)
        OntologyTerm.objects.get_or_create(id=f"HGNC:{hgnc_id}", defaults={"ontology_service": OntologyService.HGNC,
                                                                          "index": hgnc_id, "name": symbol,
                                                                          "from_import": ontology_import})


def create_ontology_test_data():
    # Suppress Pronto warnings, @see https://docs.python.org/3/library/warnings.html#temporarily-suppressing-warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        test_data_dir = Path(settings.BASE_DIR) / "ontology" / "tests" / "test_data"
        biomart_filename = test_data_dir / "biomart_omim.tsv"
        hpo_filename = test_data_dir / "small.owl"
        phenotype_to_genes_filename = test_data_dir / "phenotype_to_genes.txt"

        _create_hgnc_test_data()
        # force=False so a database the test runner already seeded (see VariantGridTestRunner) costs
        # an md5 of the file rather than another pronto parse of the OWL
        for loader, filename in [(ontology_import.load_biomart, biomart_filename),
                                 (ontology_import.load_hpo, hpo_filename),
                                 (ontology_import.load_phenotype_to_genes, phenotype_to_genes_filename)]:
            try:
                loader(str(filename), False)
            except OntologyBuilderDataUpToDateException:
                pass


def create_test_ontology_version() -> OntologyVersion:
    kwargs = {}
    now = timezone.now()
    for field, (import_source, filenames) in OntologyVersion.ONTOLOGY_IMPORTS.items():
        filename = filenames[0]
        oi = OntologyImport.objects.filter(import_source=import_source, filename=filename).first()
        if not oi:
            oi, _ = OntologyImport.objects.get_or_create(import_source=import_source, filename=filename,
                                                         defaults={"processed_date": now})
        kwargs[field] = oi
    ontology_version = OntologyVersion.objects.filter(**kwargs).first()
    if ontology_version is None:
        ontology_version, _ = OntologyVersion.objects.get_or_create(**kwargs)
    return ontology_version
