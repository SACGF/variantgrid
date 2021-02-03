import os
from collections import defaultdict, Counter

from django.core.management import BaseCommand
from django.utils import timezone

from annotation.models import GeneAnnotationVersion, OntologyImport, OntologyTerm
from genes.gene_matching import GeneMatcher
from genes.models import GeneAnnotationRelease, RVIS, GnomADGeneConstraint
from library.django_utils.django_file_utils import get_import_processing_filename
from ontology.models import OntologyService, OntologySnake
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv


class Command(BaseCommand):
    GENE_ANNOTATION_HEADER = ["version_id", "gene_id", "hpo_terms", "omim_terms",
                              "rvis_oe_ratio_percentile", "gnomad_oe_lof"]
    TERM_JOIN_STRING = " | "

    def add_arguments(self, parser):
        parser.add_argument('--gene-annotation-release', required=True, type=int)
        parser.add_argument('--force', action="store_true", help="Force create new GeneAnnotation for same release")

    def handle(self, *args, **options):
        print(f"Started: {timezone.now()}")

        gar_id = options["gene_annotation_release"]
        force = options["force"]

        gene_annotation_release = GeneAnnotationRelease.objects.get(pk=gar_id)
        if not force and GeneAnnotationVersion.objects.filter(gene_annotation_release=gene_annotation_release).exists():
            raise ValueError("Existing GeneAnnotationVersion for gene_annotation_release={} exists! Use --force?")

        for ontology_service in [OntologyService.OMIM, OntologyService.HPO, OntologyService.HGNC]:
            if not OntologyTerm.objects.filter(ontology_service=ontology_service).exists():
                raise ValueError(f"No {ontology_service.label} records - please import first")

        rvis = RVIS.objects.first()
        if rvis is None:
            raise ValueError("You need to import RVIS (see annotation page)")

        gnomad_gene_constraint = GnomADGeneConstraint.objects.first()
        if gnomad_gene_constraint is None:
            raise ValueError("You need to import gnomAD Gene Constraints (see annotation page)")

        # Only 1 of each of RVIS/Gnomad (CachedWebResource - deleted upon reload)
        gav = GeneAnnotationVersion.objects.create(gene_annotation_release=gene_annotation_release,
                                                   last_ontology_import=OntologyImport.objects.order_by("pk").last(),
                                                   rvis_import_date=rvis.cached_web_resource.created,
                                                   gnomad_import_date=gnomad_gene_constraint.cached_web_resource.created)

        # 1st we need to make sure all symbols are matched in a GeneAnnotationRelease (HGNC already done)
        gene_matcher = GeneMatcher(gene_annotation_release)
        gnomad_gene_symbols = GnomADGeneConstraint.objects.all().values_list("gene_symbol_id", flat=True)
        rvis_gene_symbols = RVIS.objects.all().values_list("gene_symbol_id", flat=True)
        gene_symbols = set(gnomad_gene_symbols) | set(rvis_gene_symbols)
        gene_matcher.match_unmatched_symbols(gene_symbols)

        missing_genes = Counter()
        # The various annotations are for different genes, so group kwargs by gene
        annotation_by_gene = defaultdict(dict)
        for hgnc_ot in OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC):
            omim_snake = OntologySnake.terms_for_gene_symbol(hgnc_ot.name, OntologyService.OMIM)
            omim_terms = self.TERM_JOIN_STRING.join((str(lt) for lt in omim_snake.leafs()))
            hpo_snake = OntologySnake.terms_for_gene_symbol(hgnc_ot.name, OntologyService.HPO)
            hpo_terms = self.TERM_JOIN_STRING.join((str(lt) for lt in hpo_snake.leafs()))

            uc_symbol = hgnc_ot.name.upper()
            genes_qs = gene_annotation_release.genes_for_symbol(uc_symbol)
            if genes_qs.exists():
                for gene in genes_qs:
                    annotation_by_gene[gene]["omim_terms"] = omim_terms
                    annotation_by_gene[gene]["hpo_terms"] = hpo_terms
            else:
                print(f"Warning: {gene_annotation_release} has no match for '{uc_symbol}' (hpo={hpo_terms}, omim={omim_terms})")
                missing_genes["ontology"] += 1

        for gene_symbol_id, oe_lof in GnomADGeneConstraint.objects.all().values_list("gene_symbol_id", "oe_lof"):
            genes_qs = gene_annotation_release.genes_for_symbol(gene_symbol_id)
            if genes_qs.exists():
                for gene in genes_qs:
                    annotation_by_gene[gene]["gnomad_oe_lof"] = oe_lof
            else:
                print(f"Warning: {gene_annotation_release} has no match for '{gene_symbol_id}' - GnomADGeneConstraint")
                missing_genes["gnomad_gene_constraints"] += 1

        for gene_symbol_id, oe_ratio_percentile in RVIS.objects.all().values_list("gene_symbol_id", "oe_ratio_percentile"):
            genes_qs = gene_annotation_release.genes_for_symbol(gene_symbol_id)
            if genes_qs.exists():
                for gene in genes_qs:
                    annotation_by_gene[gene]["rvis_oe_ratio_percentile"] = oe_ratio_percentile
            else:
                print(f"Warning: {gene_annotation_release} has no match for '{gene_symbol_id}' - RVIS")
                missing_genes["rvis"] += 1

        gene_annotation_records = []
        for gene, ga_data in annotation_by_gene.items():
            ga_data["gene_id"] = gene.pk
            ga_data["version_id"] = gav.pk
            gene_annotation_records.append(tuple((ga_data.get(k) for k in self.GENE_ANNOTATION_HEADER)))

        if gene_annotation_records:
            self._write_records(gav, gene_annotation_records)

        if missing_genes:
            print("WARNING: Could not map genes to the following records: ")
            print(missing_genes)

        print(f"Finished: {timezone.now()}")

    def _write_records(self, gene_annotation_version: GeneAnnotationVersion, gene_annotation_records):
        print(f"Creating {len(gene_annotation_records)} gene annotation records")

        separator = '\t'  # Was having trouble with quoting CSVs
        self.stdout.write("Creating CSV to insert into database...\n")
        csv_filename = get_import_processing_filename(gene_annotation_version.pk, "gene_annotation.csv",
                                                      prefix='gene_annotation')
        if os.path.exists(csv_filename):
            os.remove(csv_filename)
        write_sql_copy_csv(gene_annotation_records, csv_filename, separator=separator)
        partition_table = gene_annotation_version.get_partition_table()
        self.stdout.write(f"Inserting file '{csv_filename}' into partition {partition_table}\n")
        sql_copy_csv(csv_filename, partition_table, self.GENE_ANNOTATION_HEADER, separator=separator)
        self.stdout.write("Done!\n")
