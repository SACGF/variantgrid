import logging
import operator
import shutil
from collections import defaultdict
from typing import Tuple, Optional

from django.conf import settings
from lazy import lazy

from annotation.models.damage_enums import SIFTPrediction, FATHMMPrediction, \
    MutationAssessorPrediction, MutationTasterPrediction, Polyphen2Prediction, \
    PathogenicityImpact, ALoFTPrediction
from annotation.models.models import ColumnVEPField, VariantAnnotation, \
    VariantTranscriptAnnotation, VariantAnnotationVersion
from annotation.models.models_enums import VariantClass
from annotation.vep_annotation import VEPConfig
from genes.models import TranscriptVersion, GeneVersion
from genes.models_enums import AnnotationConsortium
from library.django_utils import get_model_fields
from library.django_utils.django_file_utils import get_import_processing_filename, get_import_processing_dir
from library.log_utils import log_traceback
from library.utils import invert_dict
from snpdb.models import GenomeBuild
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv

VEP_SEPARATOR = '&'
EMPTY_VALUES = {'', '.'}
DELIMITER = '\t'
EXTENSIONS = {",": "csv",
              "\t": "tsv"}


class VEPMissingColumnsError(RuntimeError):
    pass


class VEPColumns:
    HGVS_C = "HGVSc"  # Used to extract transcript version
    PICK = "PICK"  # Representative annotation (most-severe)
    GENE = "Gene"
    FEATURE = "Feature"
    FEATURE_TYPE = "Feature_type"


class BulkVEPVCFAnnotationInserter:
    """ Reads VEP VCF, pulling data from CSQ INFO field.
        @see https://asia.ensembl.org/info/docs/tools/vep/vep_formats.html#vcfout

        Annotation is provided for each transcript, with the representative one being
        flagged with PICK=1.
        All transcript records are copied into VariantTranscriptAnnotation
        Representative records are also copied into VariantAnnotation (used on analysis grid)

        VEP Fields are where they are copied are defined in ColumnVEPField """

    PREFIX = "annotation_run"
    DB_FIXED_COLUMNS = [
        "version_id",
        "annotation_run_id",
    ]
    DB_MANUALLY_POPULATED_COLUMNS = [
        "variant_id",
        "gene_id",
        "transcript_id",
        "transcript_version_id",
        "maxentscan_percent_diff_ref",
    ]
    DB_MANUALLY_POPULATED_VARIANT_ONLY_COLUMNS = [
        "predictions_num_pathogenic",
        "predictions_num_benign",
        "overlapping_symbols",
    ]
    DB_IGNORED_COLUMNS = ["id", "transcript"]
    VEP_NOT_COPIED_FIELDS = [
        "Allele",
        "BIOTYPE",
        "cDNA_position",
        "CDS_position",
        "CLIN_SIG",
        "COSMIC",
        "Feature",
        "Feature_type",
        "Gene",
        "HGVS_OFFSET",
        "PHENO",
        VEPColumns.PICK,
        "SOURCE",  # only populated when using --merged (which we don't use)
        "SpliceAI_pred_SYMBOL",
        "STRAND",
        "SYMBOL_SOURCE",
    ]
    VEP_NOT_COPIED_REFSEQ_ONLY = [
        "REFSEQ_MATCH",
        "REFSEQ_OFFSET",
    ]

    def __init__(self, annotation_run, infos=None, insert_variants=True, validate_columns=True):
        self.annotation_run = annotation_run
        self.rows_processed = 0
        self.variant_transcript_annotation_list = []
        self.variant_annotation_list = []
        self.variant_gene_overlap_list = []
        self.batch_id = 0
        self.insert_variants = insert_variants
        if not insert_variants:
            logging.warning("BulkVEPVCFAnnotationInserter: Not actually inserting variants")

        self.vep_columns = self._get_vep_columns_from_csq(infos)
        logging.info("CSQ: %s", self.vep_columns)
        self._setup_vep_fields_and_db_columns(validate_columns)

    @staticmethod
    def _get_vep_columns_from_csq(infos):
        """
        VEP packs annotations into a single INFO field, separated by '|', eg:
        CSQ=T|intron_variant&non_coding_transcript_variant|MODIFIER|RP5-857K21.4|ENSG00000230021

        The VCF Header has info about what CSQ means, and looks like:
        Description=Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene
        """
        columns = []
        if infos is not None:
            description = infos["CSQ"]["Description"]
            description = description.replace('"', '')  # Strip double quotes

            match = "Format: "
            columns_str = description[description.rfind(match) + len(match):]
            columns = columns_str.split("|")
        return columns

    def _add_vep_field_handlers(self):
        # TOPMED and 1k genomes can return multiple values - take highest
        format_pick_highest_float = get_clean_and_pick_single_value_func(max, float)
        format_pick_highest_int = get_clean_and_pick_single_value_func(max, int)
        remove_empty_multiples = get_clean_and_pick_single_value_func(join_uniq)
        # COSMIC v90 (5/9/2019) switched to COSV (build independent identifiers)
        extract_cosmic = get_extract_existing_variation("COSV")
        extract_dbsnp = get_extract_existing_variation("rs")

        # Some annotations return multiple results eg 2 frequencies eg "0.6764&0.2433"
        # Need to work out what to do (eg pick max)
        self.field_formatters = {
            "af_1kg": format_pick_highest_float,
            "af_uk10k": format_pick_highest_float,
            "aloft_pred": get_choice_formatter_func(ALoFTPrediction.choices),
            "aloft_high_confidence": format_aloft_high_confidence,
            "cosmic_count": format_pick_highest_int,
            "cosmic_id": extract_cosmic,
            "cosmic_legacy_id": remove_empty_multiples,
            "dbsnp_rs_id": extract_dbsnp,
            "fathmm_pred_most_damaging": get_most_damaging_func(FATHMMPrediction),
            "gnomad2_liftover_af": format_pick_highest_float,
            "gnomad_popmax": str.upper,  # nfe -> NFE
            "hgnc_id": format_hgnc_id,
            "impact": get_choice_formatter_func(PathogenicityImpact.CHOICES),
            "interpro_domain": remove_empty_multiples,
            "mastermind_count_1_cdna": get_clean_and_pick_single_value_func(operator.itemgetter(0), int),
            "mastermind_count_2_cdna_prot": get_clean_and_pick_single_value_func(operator.itemgetter(1), int),
            "mastermind_count_3_aa_change": get_clean_and_pick_single_value_func(operator.itemgetter(2), int),
            "mutation_assessor_pred_most_damaging": get_most_damaging_func(MutationAssessorPrediction),
            "mutation_taster_pred_most_damaging": get_most_damaging_func(MutationTasterPrediction),
            "nmd_escaping_variant": format_nmd_escaping_variant,
            # conservation fields are from BigWig, which can return multiple entries
            # for deletions. Higher = more conserved, so for rare disease filtering taking max makes sense
            "phastcons_100_way_vertebrate": format_pick_highest_float,
            "phastcons_30_way_mammalian": format_pick_highest_float,
            "phastcons_46_way_mammalian": format_pick_highest_float,
            "phylop_100_way_vertebrate": format_pick_highest_float,
            "phylop_30_way_mammalian": format_pick_highest_float,
            "phylop_46_way_mammalian": format_pick_highest_float,
            "polyphen2_hvar_pred_most_damaging": get_most_damaging_func(Polyphen2Prediction),
            "sift": format_vep_sift_to_choice,
            "somatic": format_vep_somatic,
            "topmed_af": format_pick_highest_float,
            "variant_class": get_choice_formatter_func(VariantClass.choices),
        }
        if self.genome_build == GenomeBuild.grch38():
            self.field_formatters["gnomad_filtered"] = gnomad_filtered_func

        self.source_field_to_columns = defaultdict(set)
        self.ignored_vep_fields = self.VEP_NOT_COPIED_FIELDS.copy()

        vc = VEPConfig(self.genome_build)
        q_columns_version = ColumnVEPField.get_columns_version_q(vc.columns_version)
        vep_source_qs = ColumnVEPField.filter_for_build(self.genome_build)

        # Sort to have consistent VCF headers
        for cvf in vep_source_qs.filter(q_columns_version).order_by("source_field"):
            try:
                if cvf.vep_custom:  # May not be configured
                    prefix = cvf.get_vep_custom_display()
                    setting_key = prefix.lower()
                    _ = vc[setting_key]  # May throw exception if not setup
                    if cvf.source_field_has_custom_prefix:
                        self.ignored_vep_fields.append(prefix)

                self.source_field_to_columns[cvf.vep_info_field].add(cvf.variant_grid_column_id)
                # logging.info("Handling column %s => %s", cvf.vep_info_field, cvf.variant_grid_column_id)
            except:
                logging.warning("Skipping custom %s due to missing settings", cvf.vep_info_field)

        vav = self.annotation_run.variant_annotation_version
        self.prediction_pathogenic_values = vav.get_functional_prediction_pathogenic_levels()

    def _setup_vep_fields_and_db_columns(self, validate_columns):
        self._add_vep_field_handlers()

        ignore_columns = set(self.DB_FIXED_COLUMNS +
                             self.DB_MANUALLY_POPULATED_COLUMNS +
                             self.DB_IGNORED_COLUMNS)
        if self.annotation_run.annotation_consortium == AnnotationConsortium.REFSEQ:
            ignore_columns.update(self.VEP_NOT_COPIED_REFSEQ_ONLY)

        # Find the ones that don't apply to this version, and exclude them
        vc = VEPConfig(self.genome_build)
        qs = ColumnVEPField.filter_for_build(self.genome_build)
        q_not_this_version = ~ColumnVEPField.get_columns_version_q(vc.columns_version)
        vep_fields_not_this_version = qs.filter(q_not_this_version).values_list("column", flat=True)
        print("*" * 40)
        print(f"{vep_fields_not_this_version=}")
        ignore_columns.update(vep_fields_not_this_version)

        for c in list(ignore_columns):
            if c.endswith("_id"):
                django_column = c.rsplit("_id", 1)[0]
                ignore_columns.add(django_column)

        self.all_variant_columns = set(get_model_fields(VariantAnnotation)) - ignore_columns
        self.transcript_columns = set(get_model_fields(VariantTranscriptAnnotation)) - ignore_columns
        self.variant_only_columns = self.all_variant_columns - self.transcript_columns

        if validate_columns:
            # Display any unused VEP columns -
            handled_vep_columns = set(self.ignored_vep_fields)
            handled_vep_columns.update(self.source_field_to_columns)

            unused_vep_columns = set(self.vep_columns) - handled_vep_columns
            if unused_vep_columns:
                logging.warning("Unhandled VEP CSQ columns (maybe add to VEP_NOT_COPIED_FIELDS or disable in VEP pipeline?) :")
                logging.warning(", ".join(sorted(unused_vep_columns)))

            missing_expected_vep_columns = handled_vep_columns - set(self.vep_columns)
            if missing_expected_vep_columns:
                missing = ", ".join(sorted(missing_expected_vep_columns))
                msg = f"Fields missing from VEP annotated vcf CSQ columns: {missing}"
                logging.error(msg)
                raise VEPMissingColumnsError(msg)

        # Display any unpopulated annotation columns (perhaps forgot map VEP fields?)
        populated_columns = set()
        for columns_set in self.source_field_to_columns.values():
            populated_columns.update(columns_set)

        unpopulated_columns = self.all_variant_columns - populated_columns
        if unpopulated_columns:
            if validate_columns:
                logging.warning("Unpopulated annotation columns:")
                logging.warning(", ".join(sorted(unpopulated_columns)))

            self.all_variant_columns -= unpopulated_columns
            self.variant_only_columns -= unpopulated_columns
            self.transcript_columns -= unpopulated_columns

        self.constant_data = {
            "version_id": self.annotation_run.variant_annotation_version.pk,
            "annotation_run_id": self.annotation_run.pk,
        }

    def get_sql_csv_header(self, base_table_name):
        if base_table_name == VariantAnnotationVersion.VARIANT_GENE_OVERLAP:
            manual_columns = []
            columns = ["variant_id", "gene_id"]
        elif base_table_name == VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION:
            manual_columns = self.DB_MANUALLY_POPULATED_COLUMNS + self.DB_MANUALLY_POPULATED_VARIANT_ONLY_COLUMNS
            columns = self.all_variant_columns
        elif base_table_name == VariantAnnotationVersion.TRANSCRIPT_ANNOTATION:
            manual_columns = self.DB_MANUALLY_POPULATED_COLUMNS
            columns = self.transcript_columns
        else:
            raise ValueError(f"Unknown VariantAnnotationVersion base_table_name: '{base_table_name}'")
        return self.DB_FIXED_COLUMNS + manual_columns + list(sorted(columns))

    def vep_to_db_dict(self, vep_transcript_data, model_columns):
        # logging.debug("vep_to_db_dict:")
        # logging.debug(vep_transcript_data)

        raw_db_data = {}

        # These can be straight copied
        for source_field, dest_columns in self.source_field_to_columns.items():
            for c in dest_columns:
                if c in model_columns:
                    raw_value = vep_transcript_data[source_field]
                    if raw_value is not None:
                        raw_db_data[c] = raw_value

        db_data = {}
        for c, value in raw_db_data.items():
            formatter = self.field_formatters.get(c)
            if formatter:
                try:
                    value = formatter(value)
                except Exception as e:
                    msg = f"Error formatting field: '{c}' formatter: '{formatter}' value: '{value}'"
                    raise RuntimeError(msg) from e
            db_data[c] = value

        return db_data

    def get_gene_id(self, vep_transcript_data):
        gene_id = None
        vep_gene = vep_transcript_data[VEPColumns.GENE]
        if vep_gene:
            if vep_gene in self.gene_identifiers:
                gene_id = vep_gene
            else:
                logging.warning(f"Could not find gene_id: '{vep_gene}'")
        return gene_id

    def get_transcript_and_version_ids(self, vep_transcript_data) -> Tuple[Optional[str], Optional[str]]:
        transcript_id = None
        transcript_version_id = None
        vep_feature_type = vep_transcript_data[VEPColumns.FEATURE_TYPE]
        if vep_feature_type == "Transcript":
            accession = vep_transcript_data[VEPColumns.FEATURE]
            if accession:
                # We pass in --transcript_version so Ensembl IDs will have version
                t_id, version = TranscriptVersion.get_transcript_id_and_version(accession)
                transcript_versions = self.transcript_versions_by_id.get(t_id)
                if transcript_versions:
                    transcript_id = t_id  # Know it's valid to link
                    transcript_version_id = transcript_versions.get(version)
                    if transcript_version_id is None:
                        logging.warning(f"Have transcript '{transcript_id}' but no version: '{version}'")
                else:
                    logging.warning(f"Could not find transcript: '{t_id}'")
        return transcript_id, transcript_version_id

    def add_calculated_transcript_columns(self, transcript_data):
        self._add_calculated_maxentscan(transcript_data)

    def add_calculated_variant_annotation_columns(self, transcript_data):
        self._add_calculated_num_predictions(transcript_data)

    def _add_calculated_num_predictions(self, transcript_data):
        num_pathogenic = 0
        num_benign = 0
        for field, path_values in self.prediction_pathogenic_values.items():
            value = transcript_data.get(field)
            if value is not None:
                if value in path_values:
                    num_pathogenic += 1
                else:
                    num_benign += 1
        transcript_data["predictions_num_pathogenic"] = num_pathogenic
        transcript_data["predictions_num_benign"] = num_benign

    @staticmethod
    def _add_calculated_maxentscan(transcript_data):
        maxentscan_diff = transcript_data.get("maxentscan_diff")
        maxentscan_ref = transcript_data.get("maxentscan_ref")
        if maxentscan_diff and maxentscan_ref:
            maxentscan_diff = float(maxentscan_diff)
            maxentscan_ref = float(maxentscan_ref)
            if maxentscan_ref:
                transcript_data["maxentscan_percent_diff_ref"] = 100 * maxentscan_diff / maxentscan_ref

    @staticmethod
    def _merge_cosmic_ids(transcript_data, custom_vcf_cosmic_id):
        """ VEP ships w/COSMIC so we always try and pull it out of the existing variation field
            We also support a custom vcf of COSMIC (for cosmic_count) - this is often more recent """
        cosmic_ids = set(custom_vcf_cosmic_id.split(VEP_SEPARATOR))
        if existing_cosmic_id := transcript_data.get("cosmic_id"):
            cosmic_ids.update(existing_cosmic_id.split(VEP_SEPARATOR))
        transcript_data["cosmic_id"] = VEP_SEPARATOR.join(sorted(cosmic_ids))

    def process_entry(self, v):
        if len(self.variant_transcript_annotation_list) >= settings.SQL_BATCH_INSERT_SIZE:
            self.bulk_insert()

        csq = v.INFO.get("CSQ")
        if csq is None:  # Legacy data error, VEP failed to annotate as something like C>C
            logging.warning(f"Skipped {v} as no CSQ")
            return

        try:
            variant_id = v.INFO["variant_id"]
            variant_data = None
            overlapping_symbols = set()
            overlapping_gene_ids = set()
            for vep_transcript_record in csq.split(","):
                vep_transcript_columns = empty_to_none(vep_transcript_record.split("|"))
                num_expected_cols = len(self.vep_columns)
                num_actual_cols = len(vep_transcript_columns)
                if num_expected_cols != num_actual_cols:
                    msg = f"Expected CSQ transcript to have {num_expected_cols} cols but had {num_actual_cols}. CSQ for transcript: '{vep_transcript_record}'"
                    raise ValueError(msg)

                vep_transcript_data = dict(zip(self.vep_columns, vep_transcript_columns))
                gene_id = self.get_gene_id(vep_transcript_data)
                transcript_data = self.vep_to_db_dict(vep_transcript_data, self.transcript_columns)
                transcript_data.update(self.constant_data)
                transcript_data["variant_id"] = variant_id
                transcript_data["gene_id"] = gene_id
                transcript_id, transcript_version_id = self.get_transcript_and_version_ids(vep_transcript_data)
                transcript_data["transcript_id"] = transcript_id
                transcript_data["transcript_version_id"] = transcript_version_id
                if symbol := transcript_data.get("symbol"):
                    overlapping_symbols.add(symbol)
                if gene_id:
                    overlapping_gene_ids.add(gene_id)
                self.add_calculated_transcript_columns(transcript_data)
                self.variant_transcript_annotation_list.append(transcript_data)

                representative_transcript = vep_transcript_data.get(VEPColumns.PICK, False)
                if representative_transcript:
                    variant_data = self.vep_to_db_dict(vep_transcript_data, self.variant_only_columns)
                    variant_data.update(transcript_data)
                    self.add_calculated_variant_annotation_columns(variant_data)
                    # If we're using custom COSMIC vcf, merge with those from VEP existing variation
                    if custom_vcf_cosmic_ids := vep_transcript_data.get("COSMIC"):
                        self._merge_cosmic_ids(variant_data, custom_vcf_cosmic_ids)

            # Handle variant data at the end so we can store overlapping symbols
            if variant_data:
                if overlapping_symbols:
                    variant_data["overlapping_symbols"] = ",".join(sorted(overlapping_symbols))
                self.variant_annotation_list.append(variant_data)

                for gene_id in overlapping_gene_ids:
                    overlapping_gene_data = {
                        "variant_id": variant_id,
                        "gene_id": gene_id,
                    }
                    overlapping_gene_data.update(self.constant_data)
                    self.variant_gene_overlap_list.append(overlapping_gene_data)

            self.rows_processed += 1
        except Exception as e:
            log_traceback()
            logging.error("Problem parsing variant: '%s'", v)
            raise e

    def finish(self):
        self.bulk_insert()

    def bulk_insert(self):
        ANNOTATION_TYPE = {
            VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION: self.variant_annotation_list,
            VariantAnnotationVersion.TRANSCRIPT_ANNOTATION: self.variant_transcript_annotation_list,
            VariantAnnotationVersion.VARIANT_GENE_OVERLAP: self.variant_gene_overlap_list,
        }

        for base_table_name, annotations_list in ANNOTATION_TYPE.items():
            logging.info("bulk_insert")
            extension = EXTENSIONS[DELIMITER]
            base_filename = f"{self.PREFIX}_{base_table_name}_{self.batch_id}.{extension}"
            data_filename = get_import_processing_filename(self.annotation_run.pk, base_filename, prefix=self.PREFIX)
            header = self.get_sql_csv_header(base_table_name)
            row_data = self._annotations_list_to_row_data(header, annotations_list)
            write_sql_copy_csv(row_data, data_filename, delimiter=DELIMITER)

            if self.insert_variants:
                vav = self.annotation_run.variant_annotation_version
                partition_table = vav.get_partition_table(base_table_name=base_table_name)

                logging.info("Inserting file '%s' into partition %s", data_filename, partition_table)
                sql_copy_csv(data_filename, partition_table, header, delimiter=DELIMITER)
                logging.info("Done!")

            annotations_list.clear()

        self.batch_id += 1

    @staticmethod
    def _annotations_list_to_row_data(header: list, annotations_list: list):
        """ returns list of lists (in order of header) """

        for data in annotations_list:
            data_values = [data.get(v) for v in header]
            yield data_values

    @staticmethod
    def _gene_overlap_to_row_data(_header, annotations_list: list):
        pass

    def remove_processing_files(self):
        import_processing_dir = get_import_processing_dir(self.annotation_run.pk, prefix=self.PREFIX)
        logging.info("********* Deleting '%s' *******", import_processing_dir)
        shutil.rmtree(import_processing_dir)

    @property
    def genome_build(self):
        vav = self.annotation_run.variant_annotation_version
        return vav.genome_build

    @lazy
    def gene_identifiers(self):
        """ A set of gene identifiers to check whether in DB or not """
        vav = self.annotation_run.variant_annotation_version
        gv_qs = GeneVersion.objects.filter(gene__annotation_consortium=vav.annotation_consortium,
                                           genome_build=vav.genome_build)
        return set(gv_qs.distinct().values_list("gene", flat=True))

    @lazy
    def transcript_versions_by_id(self):
        vav = self.annotation_run.variant_annotation_version
        return TranscriptVersion.transcript_versions_by_id(vav.genome_build, vav.annotation_consortium)


def empty_to_none(it):
    return [v if v not in EMPTY_VALUES else None
            for v in it]


# Field formatters
def gnomad_filtered_func(raw_value):
    """ We use FILTER in Gnomad3 (GRCh38 only) - need to convert back to bool """
    return raw_value not in (None, "PASS")


def format_hgnc_id(raw_value):
    """ VEP GRCh37 returns 55 while GRCh38 returns "HGNC:55" """
    return int(raw_value.replace("HGNC:", ""))


def format_vep_sift_to_choice(vep_sift):
    """ we ignore the low_confidence calls to make it simpler """
    if vep_sift.startswith("deleterious"):
        return SIFTPrediction.DAMAGING
    elif vep_sift.startswith("tolerated"):
        return SIFTPrediction.TOLERATED
    raise ValueError(f"Unknown SIFT value: '{vep_sift}'")


def get_extract_existing_variation(prefix):
    def format_vep_existing_variation(vep_existing_variation):
        ev_list = vep_existing_variation.split(VEP_SEPARATOR)
        cosmic_ids = [ev for ev in ev_list if ev.startswith(prefix)]
        return VEP_SEPARATOR.join(sorted(set(cosmic_ids)))

    return format_vep_existing_variation


def format_vep_somatic(raw_value):
    return "1" in raw_value


def get_choice_formatter_func(choices):
    lookup = invert_dict(dict(choices))

    def format_choice(raw_value):
        return lookup[raw_value]

    return format_choice


def get_clean_and_pick_single_value_func(pick_single_value_func, cast_func=None):
    """ Returns a function to clean and pick single value.
        casting is performed before calling pick_single_value_func so you can call min/max """

    def _clean_and_pick_single_value_func(raw_value):
        it = (tm for tm in raw_value.split(VEP_SEPARATOR) if tm != '')
        # Handle '.'
        if cast_func:
            values = [cast_func(v) for v in it if v not in EMPTY_VALUES]
        else:
            values = [v for v in it if v not in EMPTY_VALUES]
        value = None
        if values:
            value = pick_single_value_func(values)
        return value

    return _clean_and_pick_single_value_func


def join_uniq(multiple_values):
    return VEP_SEPARATOR.join(set(multiple_values))


# Field Processors
def get_most_damaging_func(klass):
    def get_most_damaging(multiple_values):
        prediction_list = multiple_values.split(VEP_SEPARATOR)
        return klass.get_most_damaging(prediction_list)

    return get_most_damaging


def format_nmd_escaping_variant(value) -> bool:
    return value == "NMD_escaping_variant"


def format_aloft_high_confidence(value) -> bool:
    high_confidence = None
    if value == "High Confidence":
        high_confidence = True
    elif value == "Low Confidence":
        high_confidence = False
    return high_confidence
