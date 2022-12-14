import os
from collections import defaultdict, Counter

from django.conf import settings
from django.core.management import BaseCommand
from django.utils import timezone

from annotation.models import GeneAnnotationVersion, OntologyTerm, GenomeBuild, AnnotationVersion, \
    InvalidAnnotationVersionError, GeneAnnotation, DBNSFPGeneAnnotationVersion
from genes.gene_matching import ReleaseGeneMatcher
from genes.models import GeneAnnotationRelease, GnomADGeneConstraint, ReleaseGeneSymbolGene, Gene
from library.django_utils.django_file_utils import get_import_processing_filename
from ontology.models import OntologyService, GeneDiseaseClassification, OntologyTermRelation, \
    OntologyVersion
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv


class Command(BaseCommand):
    GENE_ANNOTATION_HEADER = ["version_id", "gene_id", "dbnsfp_gene_id", "hpo_terms", "omim_terms", "mondo_terms",
                              "gene_disease_moderate_or_above", "gene_disease_supportive_or_below", "gnomad_oe_lof"]
    TERM_JOIN_STRING = " | "

    def add_arguments(self, parser):
        gar_ov_dbsnfp_help = "gene-annotation-release, ontology-version and dbnsfp-gene-version are optional but must" \
                             "be specified together"

        parser.add_argument('--force', action="store_true", help="Force create new GeneAnnotation for same release")
        parser.add_argument('--ontology-version', type=int, help=gar_ov_dbsnfp_help)
        parser.add_argument('--dbnsfp-gene-version', help=gar_ov_dbsnfp_help)

        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--gene-annotation-release', type=int, help=gar_ov_dbsnfp_help)
        group.add_argument('--missing', action="store_true",
                           help="Automatically create for latest AnnotationVersions for each build if missing")
        group.add_argument('--add-new-to-existing', action="store_true",
                           help="Add new columns gene/disease and MONDO terms to existing gene annotation")
        group.add_argument('--add-dbnsfp-gene', action="store_true",
                           help="Add new dbNSFP to existing gene annotation")
        group.add_argument('--fix-bad-gene-annotation', action="store_true",
                           help="Fix for bad gene IDs, dbNSFP not linked")

    def handle(self, *args, **options):
        if not settings.ANNOTATION_GENE_ANNOTATION_VERSION_ENABLED:
            raise ValueError("settings.ANNOTATION_GENE_ANNOTATION_VERSION is disabled.")
        print(f"Started: {timezone.now()}")

        force = options["force"]
        gar_id = options["gene_annotation_release"]
        ov_id = options["ontology_version"]
        dbnsfp_gene_version_id = options["dbnsfp_gene_version"]

        missing = options["missing"]
        self._validate_has_required_data()

        if options["add_new_to_existing"]:
            self._add_new_columns_to_existing()
            return

        dbnsfp_gene_version = DBNSFPGeneAnnotationVersion.objects.all().order_by("created").last()
        if dbnsfp_gene_version is None:
            raise InvalidAnnotationVersionError("You need to import DBNSFPGeneAnnotationVersion")

        if options["add_dbnsfp_gene"]:
            self._add_dbnsfp_gene(dbnsfp_gene_version)
            return

        gene_symbols = set(GnomADGeneConstraint.objects.all().values_list("gene_symbol_id", flat=True))

        if options["fix_bad_gene_annotation"]:
            self._fix_bad_gene_annotation(gene_symbols)
            return

        if gar_id:
            if not ov_id:
                raise ValueError("You must specify ontology-version when gene-annotation-release is specified")
            ontology_version = OntologyVersion.objects.get(pk=ov_id)
            if dbnsfp_gene_version_id:
                dbnsfp_gene_version = DBNSFPGeneAnnotationVersion.objects.get(pk=dbnsfp_gene_version_id)
            gene_annotation_release = GeneAnnotationRelease.objects.get(pk=gar_id)
            if not force and GeneAnnotationVersion.objects.filter(gene_annotation_release=gene_annotation_release).exists():
                raise ValueError("Existing GeneAnnotationVersion for gene_annotation_release={} exists! Use --force?")
            gav = self._create_gene_annotation_version(gene_annotation_release, ontology_version,
                                                       dbnsfp_gene_version)
            self._populate_gene_annotation_version(gav, gene_symbols)
        elif missing:
            if ov_id is not None:
                raise ValueError("Only specify ontology-version when gene-annotation-release also specified")

            for genome_build in GenomeBuild.builds_with_annotation():
                av = AnnotationVersion.latest(genome_build, validate=False, active=False)
                if not av:
                    raise InvalidAnnotationVersionError(f"No AnnotationVersion for {genome_build}")

                if not av.variant_annotation_version:
                    raise InvalidAnnotationVersionError(f"AnnotationVersion {av} has no VariantAnnotationVersion set")

                if not av.ontology_version:
                    raise InvalidAnnotationVersionError(f"AnnotationVersion {av} has no OntologyVersion set")

                if av.gene_annotation_version:
                    if av.ontology_version == av.gene_annotation_version.ontology_version:
                        print(f"Skipping {av} - already has GeneAnnotation")
                        continue
                    else:
                        print("AV ontology version != existing gene annotation ontology version")

                gar = av.variant_annotation_version.gene_annotation_release
                if not gar:
                    raise InvalidAnnotationVersionError(f"VariantAnnotationVersion {av.variant_annotation_version} needs to be assigned an GeneAnnotationRelease")

                gav = self._create_gene_annotation_version(gar, av.ontology_version, dbnsfp_gene_version)
                self._populate_gene_annotation_version(gav, gene_symbols)

    @staticmethod
    def _validate_has_required_data():
        for ontology_service in [OntologyService.OMIM, OntologyService.HPO,
                                 OntologyService.MONDO, OntologyService.HGNC]:
            if not OntologyTerm.objects.filter(ontology_service=ontology_service).exists():
                raise ValueError(f"No {ontology_service.label} records - please import first")

        if not GnomADGeneConstraint.objects.exists():
            raise ValueError("You need to import gnomAD Gene Constraints (see annotation page)")

        otr_qs = OntologyTermRelation.objects.filter(extra__strongest_classification__isnull=False)
        if not otr_qs.exists():
            raise ValueError("You need to import GenCC gene/disease curation (see annotation page)")

    def _add_dbnsfp_gene(self, dbnsfp_gene_version):
        for ga_version in GeneAnnotationVersion.objects.filter(dbnsfp_gene_version__isnull=True):
            print(f"Add dbNSFP gene annotation to {ga_version}")
            ga_by_gene_id = {}
            for ga in ga_version.geneannotation_set.all():
                ga_by_gene_id[ga.gene_id] = ga

            rgsg_qs = ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=ga_version.gene_annotation_release)
            symbol_to_gene_ids = defaultdict(set)
            for gene_symbol, gene_id in rgsg_qs.values_list("release_gene_symbol__gene_symbol", "gene_id"):
                symbol_to_gene_ids[gene_symbol].add(gene_id)

            # Go through and match all the
            gene_annotation = []
            dbnsfp_qs = dbnsfp_gene_version.dbnsfpgeneannotation_set.all()
            for dbnsfp_gene_id, gene_symbol_id in dbnsfp_qs.values_list("pk", "gene_symbol_id"):
                for gene_id in symbol_to_gene_ids.get(gene_symbol_id, []):
                    if ga := ga_by_gene_id.pop(gene_id, None):
                        ga.dbnsfp_gene_id = dbnsfp_gene_id
                        gene_annotation.append(ga)

            if gene_annotation:
                print(f"Updating {len(gene_annotation)} records....")
                GeneAnnotation.objects.bulk_update(gene_annotation, ["dbnsfp_gene_id"], batch_size=2000)

    def _add_new_columns_to_existing(self):
        """ As we only added not changed columns, can just populate existing annotation """
        UPDATE_COLUMNS = ["mondo_terms", "gene_disease_moderate_or_above", "gene_disease_supportive_or_below"]

        print("Updating existing gene annotation...")
        print("Loading HGNC data...")

        # Then go through
        for ga_version in GeneAnnotationVersion.objects.all():
            print(f"Updating GeneAnnotationVersion: {ga_version}")
            ga_by_gene_id = {}
            for ga in ga_version.geneannotation_set.all():
                ga_by_gene_id[ga.gene_id] = ga
            print("Loaded existing gene annotation, matching to HGNC")

            hgnc_data = defaultdict(dict)
            for hgnc_ot in OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC):
                gene_symbol = hgnc_ot.name
                snake = ga_version.ontology_version.terms_for_gene_symbol(hgnc_ot.name, OntologyService.MONDO,
                                                                          max_depth=0)
                uc_symbol = gene_symbol.upper()
                hgnc_data[uc_symbol]["mondo_terms"] = self.TERM_JOIN_STRING.join((str(lt) for lt in snake.leafs()))
                hgnc_data[uc_symbol]["gene_disease"] = self._get_gene_disease(gene_symbol, Command.TERM_JOIN_STRING)

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
    def _bad_gene_annotation():
        bad_gene_annotation = []
        for gav in GeneAnnotationVersion.objects.all():
            if gav.dbnsfp_gene_version and not gav.geneannotation_set.filter(dbnsfp_gene__isnull=False).exists():
                # print(f"{gav} - dbNSFP set but not linked correctly!")
                bad_gene_annotation.append(gav)
                continue

            try:
                if ga := gav.geneannotation_set.first():
                    _ = ga.gene
            except Gene.DoesNotExist:
                bad_gene_annotation.append(gav)
        return bad_gene_annotation

    def _fix_bad_gene_annotation(self, gene_symbols):
        """ dbNSFP wasn't linked from new - due to not including headers
            gene was inserted as str(gene) which caused it to be "ENSG (Ensembl)"
        """

        fixed = []
        for gav in self._bad_gene_annotation():
            print(f"{gav} - needs fixing. Deleting old...")
            gav.geneannotation_set.all().delete()
            print(f"{gav} - regenerating...")
            self._populate_gene_annotation_version(gav, gene_symbols)
            fixed.append(gav)

        print("Fixed gene annotation versions:")
        for gav in fixed:
            print(f"{gav.pk}: {gav}")

    @staticmethod
    def _get_gene_disease(ontology_version, gene_symbol, delimiter: str):
        moderate_or_above = GeneDiseaseClassification.get_above_min(GeneDiseaseClassification.MODERATE)
        supportive_or_below = [gdc.label for gdc in reversed(GeneDiseaseClassification)
                               if gdc.label not in moderate_or_above]

        diseases_supportive_or_below = []
        diseases_moderate_or_above = []
        for otr in ontology_version.gene_disease_relations(gene_symbol,
                                                           min_classification=GeneDiseaseClassification.LIMITED):
            disease = otr.source_term.name
            moi_classifications = otr.get_gene_disease_moi_classifications()
            moi_supportive_or_below = otr.get_moi_summary(moi_classifications, supportive_or_below)
            moi_moderate_or_above = otr.get_moi_summary(moi_classifications, moderate_or_above)

            if moi_supportive_or_below:
                diseases_supportive_or_below.append(f"{disease}={', '.join(sorted(moi_supportive_or_below))}")
            if moi_moderate_or_above:
                diseases_moderate_or_above.append(f"{disease}={', '.join(sorted(moi_moderate_or_above))}")

        gene_disease_supportive_or_below = delimiter.join(diseases_supportive_or_below)
        gene_disease_moderate_or_above = delimiter.join(diseases_moderate_or_above)
        return gene_disease_supportive_or_below, gene_disease_moderate_or_above

    def _create_gene_annotation_version(self, gene_annotation_release, ontology_version,
                                        dbnsfp_gene_version) -> GeneAnnotationVersion:
        gnomad_gene_constraint = GnomADGeneConstraint.objects.first()  # Only ever 1

        # Only 1 of each of Gnomad (CachedWebResource - deleted upon reload)
        # When you create GeneAnnotationVersion (sub version) it automatically creates/bumps a new annotation version
        return GeneAnnotationVersion.objects.create(gene_annotation_release=gene_annotation_release,
                                                    ontology_version=ontology_version,
                                                    dbnsfp_gene_version=dbnsfp_gene_version,
                                                    gnomad_import_date=gnomad_gene_constraint.cached_web_resource.created)

    def _populate_gene_annotation_version(self, gav: GeneAnnotationVersion, gene_symbols):

        # 1st we need to make sure all symbols are matched in a GeneAnnotationRelease (HGNC already done)
        gene_matcher = ReleaseGeneMatcher(gav.gene_annotation_release)
        gene_matcher.match_unmatched_symbols(gene_symbols)

        missing_genes = Counter()
        # The various annotations are for different genes, so group kwargs by gene
        annotation_by_gene = defaultdict(dict)
        for hgnc_ot in OntologyTerm.objects.filter(ontology_service=OntologyService.HGNC):
            service_terms = {}
            gene_symbol = hgnc_ot.name
            for ontology_service in [OntologyService.OMIM, OntologyService.HPO, OntologyService.MONDO]:
                snake = gav.ontology_version.terms_for_gene_symbol(gene_symbol, ontology_service, max_depth=0)
                service_terms[ontology_service] = self.TERM_JOIN_STRING.join((str(lt) for lt in snake.leafs()))

            gene_disease_supportive_or_below, gene_disease_moderate_or_above = self._get_gene_disease(gav.ontology_version,
                                                                                                      gene_symbol,
                                                                                                      Command.TERM_JOIN_STRING)
            if not (any(service_terms.values()) or
                    gene_disease_supportive_or_below or gene_disease_moderate_or_above):
                continue  # Skip who cares

            uc_symbol = gene_symbol.upper()
            genes_qs = gav.gene_annotation_release.genes_for_symbol(uc_symbol)
            if genes_qs.exists():
                for gene_id in genes_qs.values_list("identifier", flat=True):
                    for ontology_service, terms in service_terms.items():
                        column = f"{str(ontology_service).lower()}_terms"
                        annotation_by_gene[gene_id][column] = terms

                    annotation_by_gene[gene_id]["gene_disease_supportive_or_below"] = gene_disease_supportive_or_below
                    annotation_by_gene[gene_id]["gene_disease_moderate_or_above"] = gene_disease_moderate_or_above
            else:
                print(f"Warning: {gav.gene_annotation_release} has no match for '{uc_symbol}' ({service_terms})")
                missing_genes["ontology"] += 1

        rgsg_qs = ReleaseGeneSymbolGene.objects.filter(release_gene_symbol__release=gav.gene_annotation_release)
        symbol_to_gene_ids = defaultdict(set)
        for gene_symbol, gene_id in rgsg_qs.values_list("release_gene_symbol__gene_symbol", "gene_id"):
            symbol_to_gene_ids[gene_symbol].add(gene_id)

        # Go through and match all the
        dbnsfp_qs = gav.dbnsfp_gene_version.dbnsfpgeneannotation_set.all()
        if not dbnsfp_qs.exists():
            raise ValueError(f"{gav.dbnsfp_gene_version} has no entries!")
        for dbnsfp_gene_id, gene_symbol_id in dbnsfp_qs.values_list("pk", "gene_symbol_id"):
            if gene_ids_for_symbol := symbol_to_gene_ids.get(gene_symbol_id):
                for gene_id in gene_ids_for_symbol:
                    annotation_by_gene[gene_id]["dbnsfp_gene_id"] = dbnsfp_gene_id
            else:
                print(f"Warning: {gav.gene_annotation_release} has no match for '{gene_symbol_id}' - dbNSFP Gene Annotation")
                missing_genes["dbnsfp_gene_annotation"] += 1

        for gene_symbol_id, oe_lof in GnomADGeneConstraint.objects.all().values_list("gene_symbol_id", "oe_lof"):
            if gene_ids_for_symbol := symbol_to_gene_ids.get(gene_symbol_id):
                for gene_id in gene_ids_for_symbol:
                    annotation_by_gene[gene_id]["gnomad_oe_lof"] = oe_lof
            else:
                print(f"Warning: {gav.gene_annotation_release} has no match for '{gene_symbol_id}' - GnomADGeneConstraint")
                missing_genes["gnomad_gene_constraints"] += 1

        gene_annotation_records = []
        for gene_id, ga_data in annotation_by_gene.items():
            ga_data["gene_id"] = gene_id
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
