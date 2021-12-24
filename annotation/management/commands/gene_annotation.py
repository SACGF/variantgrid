import os
from collections import defaultdict, Counter

from django.core.management import BaseCommand
from django.utils import timezone

from annotation.models import GeneAnnotationVersion, OntologyImport, OntologyTerm, GenomeBuild, AnnotationVersion, \
    InvalidAnnotationVersionError, GeneAnnotation
from genes.gene_matching import GeneMatcher
from genes.models import GeneAnnotationRelease, GnomADGeneConstraint
from library.django_utils.django_file_utils import get_import_processing_filename
from ontology.models import OntologyService, OntologySnake, GeneDiseaseClassification
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv


class Command(BaseCommand):
    GENE_ANNOTATION_HEADER = ["version_id", "gene_id", "hpo_terms", "omim_terms", "mondo_terms",
                              "gene_disease_moderate_or_above", "gene_disease_supportive_or_below", "gnomad_oe_lof"]
    TERM_JOIN_STRING = " | "

    def add_arguments(self, parser):
        parser.add_argument('--force', action="store_true", help="Force create new GeneAnnotation for same release")
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--gene-annotation-release', type=int)
        group.add_argument('--missing', action="store_true",
                           help="Automatically create for latest AnnotationVersions for each build if missing")
        group.add_argument('--add-new-to-existing', action="store_true",
                           help="Add new columns gene/disease and MONDO terms to existing gene annotation")

    def handle(self, *args, **options):
        print(f"Started: {timezone.now()}")

        force = options["force"]
        gar_id = options["gene_annotation_release"]
        missing = options["missing"]
        add_new_to_existing = options["add_new_to_existing"]
        self._validate_has_required_data()

        if add_new_to_existing:
            self._add_new_columns_to_existing()
            return

        gene_symbols = set(GnomADGeneConstraint.objects.all().values_list("gene_symbol_id", flat=True))

        if gar_id:
            gene_annotation_release = GeneAnnotationRelease.objects.get(pk=gar_id)
            if not force and GeneAnnotationVersion.objects.filter(gene_annotation_release=gene_annotation_release).exists():
                raise ValueError("Existing GeneAnnotationVersion for gene_annotation_release={} exists! Use --force?")
            self._create_gene_annotation_version(gene_annotation_release, gene_symbols)
        else:
            for genome_build in GenomeBuild.builds_with_annotation():
                av = AnnotationVersion.latest(genome_build, validate=False)
                if av.gene_annotation_version:
                    print(f"Skipping {av} - already has GeneAnnotation")
                    continue

                if not av.variant_annotation_version:
                    raise InvalidAnnotationVersionError(f"AnnotationVersion {av} has no VariantAnnotationVersion set")

                gar = av.variant_annotation_version.gene_annotation_release
                if not gar:
                    raise InvalidAnnotationVersionError(f"VariantAnnotationVersion {av.variant_annotation_version} needs to be assinged an GeneAnnotationRelease")

                self._create_gene_annotation_version(gar, gene_symbols)

    @staticmethod
    def _validate_has_required_data():
        for ontology_service in [OntologyService.OMIM, OntologyService.HPO, OntologyService.HGNC]:
            if not OntologyTerm.objects.filter(ontology_service=ontology_service).exists():
                raise ValueError(f"No {ontology_service.label} records - please import first")

        if not GnomADGeneConstraint.objects.exists():
            raise ValueError("You need to import gnomAD Gene Constraints (see annotation page)")

    def _add_new_columns_to_existing(self):
        """ As we only added not changed columns, can just populate existing annotation """
        UPDATE_COLUMNS = ["mondo_terms", "gene_disease_moderate_or_above", "gene_disease_supportive_or_below"]

        print("Updating existing gene annotation...")
        print("Loading HGNC data...")

        # Loop through HGNC - save into a dict
        hgnc_data = defaultdict(dict)
        for hgnc_ot in OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC):
            gene_symbol = hgnc_ot.name
            snake = OntologySnake.terms_for_gene_symbol(hgnc_ot.name, OntologyService.MONDO, max_depth=0)
            uc_symbol = gene_symbol.upper()
            hgnc_data[uc_symbol]["mondo_terms"] = self.TERM_JOIN_STRING.join((str(lt) for lt in snake.leafs()))
            hgnc_data[uc_symbol]["gene_disease"] = self._get_gene_disease(gene_symbol, Command.TERM_JOIN_STRING)

        # Then go through
        for ga_version in GeneAnnotationVersion.objects.all():
            print(f"Updating GeneAnnotationVersion: {ga_version}")
            ga_by_gene_id = {}
            for ga in ga_version.geneannotation_set.all():
                ga_by_gene_id[ga.gene_id] = ga
            print("Loaded existing gene annotation, matching to HGNC")

            gene_annotation = []
            for uc_symbol, data in hgnc_data.items():
                for gene_id in ga_version.gene_annotation_release.genes_for_symbol(uc_symbol).values_list("pk", flat=True):
                    if ga := ga_by_gene_id.get(gene_id):
                        ga.mondo_terms = data["mondo_terms"]
                        ga.gene_disease_supportive_or_below = data["gene_disease"][0]
                        ga.gene_disease_moderate_or_above = data["gene_disease"][1]
                        gene_annotation.append(ga)

            if gene_annotation:
                print(f"Updating {len(gene_annotation)} records....")
                GeneAnnotation.objects.bulk_update(gene_annotation, UPDATE_COLUMNS, batch_size=2000)

            print("*" * 20)

        print(f"Completed: {timezone.now()}")

    @staticmethod
    def _get_gene_disease(gene_symbol, delimiter: str):
        moderate_or_above = GeneDiseaseClassification.get_above_min(GeneDiseaseClassification.MODERATE)
        supportive_or_below = [gdc.label for gdc in reversed(GeneDiseaseClassification)
                               if gdc.label not in moderate_or_above]

        diseases_supportive_or_below = []
        diseases_moderate_or_above = []
        for otr in OntologySnake.gene_disease_relations(gene_symbol,
                                                        min_classification=GeneDiseaseClassification.LIMITED):
            disease = otr.source_term.name
            moi_classifications = defaultdict(lambda: defaultdict(set))
            for source in otr.extra["sources"]:
                moi = source["mode_of_inheritance"]
                classification = source["gencc_classification"]
                submitter = source["submitter"]
                moi_classifications[moi][classification].add(submitter)

            moi_supportive_or_below = []
            moi_moderate_or_above = []
            for moi, classifications in moi_classifications.items():
                classification_submitters_supportive_or_below = []
                for classification in supportive_or_below:
                    if submitters := classifications.get(classification):
                        classification_submitters = f"{classification}: {'/'.join(sorted(submitters))}"
                        classification_submitters_supportive_or_below.append(classification_submitters)
                if classification_submitters_supportive_or_below:
                    moi_supportive_or_below.append(f"{moi} ({' '.join(classification_submitters_supportive_or_below)})")

                classification_submitters_moderate_or_above = []
                for classification in moderate_or_above:
                    if submitters := classifications.get(classification):
                        classification_submitters = f"{classification}: {'/'.join(sorted(submitters))}"
                        classification_submitters_moderate_or_above.append(classification_submitters)
                if classification_submitters_moderate_or_above:
                    moi_moderate_or_above.append(f"{moi} ({' '.join(classification_submitters_moderate_or_above)})")

            if moi_supportive_or_below:
                diseases_supportive_or_below.append(f"{disease}={', '.join(sorted(moi_supportive_or_below))}")
            if moi_moderate_or_above:
                diseases_moderate_or_above.append(f"{disease}={', '.join(sorted(moi_moderate_or_above))}")

        gene_disease_supportive_or_below = delimiter.join(diseases_supportive_or_below)
        gene_disease_moderate_or_above = delimiter.join(diseases_moderate_or_above)
        return gene_disease_supportive_or_below, gene_disease_moderate_or_above

    def _create_gene_annotation_version(self, gene_annotation_release, gene_symbols):
        gnomad_gene_constraint = GnomADGeneConstraint.objects.first()

        # Only 1 of each of Gnomad (CachedWebResource - deleted upon reload)
        gav = GeneAnnotationVersion.objects.create(gene_annotation_release=gene_annotation_release,
                                                   last_ontology_import=OntologyImport.objects.order_by("pk").last(),
                                                   gnomad_import_date=gnomad_gene_constraint.cached_web_resource.created)

        # 1st we need to make sure all symbols are matched in a GeneAnnotationRelease (HGNC already done)
        gene_matcher = GeneMatcher(gene_annotation_release)
        gene_matcher.match_unmatched_symbols(gene_symbols)

        missing_genes = Counter()
        # The various annotations are for different genes, so group kwargs by gene
        annotation_by_gene = defaultdict(dict)
        for hgnc_ot in OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC):
            service_terms = {}
            gene_symbol = hgnc_ot.name
            for ontology_service in [OntologyService.OMIM, OntologyService.HPO, OntologyService.MONDO]:
                snake = OntologySnake.terms_for_gene_symbol(gene_symbol, ontology_service, max_depth=0)
                service_terms[ontology_service] = self.TERM_JOIN_STRING.join((str(lt) for lt in snake.leafs()))

            gene_disease_supportive_or_below, gene_disease_moderate_or_above = self._get_gene_disease(gene_symbol, Command.TERM_JOIN_STRING)
            if not (any(service_terms.values()) or
                    gene_disease_supportive_or_below or gene_disease_moderate_or_above):
                continue  # Skip who cares

            uc_symbol = gene_symbol.upper()
            genes_qs = gene_annotation_release.genes_for_symbol(uc_symbol)
            if genes_qs.exists():
                for gene in genes_qs:
                    for ontology_service, terms in service_terms.items():
                        column = f"{str(ontology_service).lower()}_terms"
                        annotation_by_gene[gene][column] = terms

                    annotation_by_gene[gene]["gene_disease_supportive_or_below"] = gene_disease_supportive_or_below
                    annotation_by_gene[gene]["gene_disease_moderate_or_above"] = gene_disease_moderate_or_above
            else:
                print(f"Warning: {gene_annotation_release} has no match for '{uc_symbol}' ({service_terms})")
                missing_genes["ontology"] += 1

        for gene_symbol_id, oe_lof in GnomADGeneConstraint.objects.all().values_list("gene_symbol_id", "oe_lof"):
            genes_qs = gene_annotation_release.genes_for_symbol(gene_symbol_id)
            if genes_qs.exists():
                for gene in genes_qs:
                    annotation_by_gene[gene]["gnomad_oe_lof"] = oe_lof
            else:
                print(f"Warning: {gene_annotation_release} has no match for '{gene_symbol_id}' - GnomADGeneConstraint")
                missing_genes["gnomad_gene_constraints"] += 1

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

        delimiter = '\t'  # Was having trouble with quoting CSVs
        self.stdout.write("Creating CSV to insert into database...\n")
        csv_filename = get_import_processing_filename(gene_annotation_version.pk, "gene_annotation.csv",
                                                      prefix='gene_annotation')
        if os.path.exists(csv_filename):
            os.remove(csv_filename)
        write_sql_copy_csv(gene_annotation_records, csv_filename, delimiter=delimiter)
        partition_table = gene_annotation_version.get_partition_table()
        self.stdout.write(f"Inserting file '{csv_filename}' into partition {partition_table}\n")
        sql_copy_csv(csv_filename, partition_table, self.GENE_ANNOTATION_HEADER, delimiter=delimiter)
        self.stdout.write("Done!\n")
