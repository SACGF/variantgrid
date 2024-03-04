import logging
import operator
import shutil
from collections import defaultdict
from functools import cached_property, reduce
from typing import Optional

from django.conf import settings
from django.db.models import Q

from annotation.models.damage_enums import SIFTPrediction, FATHMMPrediction, \
    MutationAssessorPrediction, MutationTasterPrediction, Polyphen2Prediction, \
    PathogenicityImpact, ALoFTPrediction, AlphaMissensePrediction
from annotation.models.models import ColumnVEPField, VariantAnnotation, \
    VariantTranscriptAnnotation, VariantAnnotationVersion, VariantGeneOverlap, AnnotationRun
from annotation.models.models_enums import VariantClass, VariantAnnotationPipelineType
from annotation.vep_annotation import VEPConfig
from genes.hgvs import HGVSMatcher, HGVSException
from genes.models import TranscriptVersion, GeneVersion
from genes.models_enums import AnnotationConsortium
from library.django_utils import get_model_fields
from library.django_utils.django_file_utils import get_import_processing_filename, get_import_processing_dir
from library.log_utils import log_traceback
from library.utils import invert_dict
from snpdb.models import GenomeBuild, VariantCoordinate
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
        "gnomad_hemi_count",
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
        # "SpliceAI_pred_SYMBOL",
        "STRAND",
        "SYMBOL_SOURCE",
    ]
    VEP_NOT_COPIED_REFSEQ_ONLY = [
        "REFSEQ_MATCH",
        "REFSEQ_OFFSET",
    ]

    # These are present in columns_version > 2
    ALOFT_COLUMNS = ['aloft_prob_tolerant', 'aloft_prob_recessive', 'aloft_prob_dominant',
                     'aloft_pred', 'aloft_high_confidence', 'aloft_ensembl_transcript']

    def __init__(self, annotation_run: AnnotationRun, infos=None, insert_variants=True, validate_columns=True):
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
        self.aloft_columns = False
        logging.info("CSQ: %s", self.vep_columns)
        self._setup_vep_fields_and_db_columns(validate_columns)
        self.hgvs_matcher = HGVSMatcher(annotation_run.genome_build)

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
        empty_mave_float_values = EMPTY_VALUES | {"NA"}
        format_pick_lowest_float = get_clean_and_pick_single_value_func(min, float,
                                                                        empty_values=empty_mave_float_values)
        format_pick_highest_float = get_clean_and_pick_single_value_func(max, float)
        format_pick_highest_int = get_clean_and_pick_single_value_func(max, int)
        remove_empty_multiples = get_clean_and_pick_single_value_func(join_uniq)
        # COSMIC v90 (5/9/2019) switched to COSV (build independent identifiers)
        extract_cosmic = get_extract_existing_variation("COSV")
        extract_dbsnp = get_extract_existing_variation("rs")
        format_empty_as_none = get_format_empty_as_none(empty_values=EMPTY_VALUES)

        # Some annotations return multiple results e.g. 2 frequencies e.g. "0.6764&0.2433"
        # Need to work out what to do (eg pick max)
        self.field_formatters = {
            "af_1kg": format_pick_highest_float,
            "af_uk10k": format_pick_highest_float,
            # ALoFT comes as multiple values, so "." won't have already been ignored, so need to handle
            "aloft_prob_tolerant": format_empty_as_none,
            "aloft_prob_recessive": format_empty_as_none,
            "aloft_prob_dominant": format_empty_as_none,
            "aloft_pred": get_choice_formatter_func(ALoFTPrediction.choices, empty_values=["."]),
            "aloft_high_confidence": format_aloft_high_confidence,
            "aloft_ensembl_transcript": format_empty_as_none,
            "alphamissense_class": get_format_alphamissense_class_func(),
            "canonical": format_canonical,
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
            "mavedb_score": format_pick_lowest_float,
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

        vc = VEPConfig(self.genome_build)
        # gnomad3 wasn't combined using gnomad_data.py so just uses FILTER
        # while combined exome/genomes use "gnomad_filtered=1" (which should auto-convert bool)
        if self.genome_build == GenomeBuild.grch38() and vc.columns_version <= 2:
            self.field_formatters["gnomad_filtered"] = gnomad_filtered_func

        self.source_field_to_columns = defaultdict(set)
        self.ignored_vep_fields = self.VEP_NOT_COPIED_FIELDS.copy()

        cvf_filters = [ColumnVEPField.get_columns_version_q(vc.columns_version)]
        if self.annotation_run.pipeline_type == VariantAnnotationPipelineType.CNV:
            cvf_filters.extend([
                Q(vep_plugin__isnull=True),
                Q(vep_custom__isnull=True),
            ])
        vep_source_qs = ColumnVEPField.filter_for_build(self.genome_build)
        q_cvf = reduce(operator.and_, cvf_filters)
        # Sort to have consistent VCF headers
        for cvf in vep_source_qs.filter(q_cvf).order_by("source_field"):
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
        self.prediction_pathogenic_funcs = vav.get_pathogenic_prediction_funcs()

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
            ignored_vep_columns = set(self.ignored_vep_fields)
            vep_columns = set(self.vep_columns) - ignored_vep_columns
            handled_vep_columns = set(self.source_field_to_columns)

            unused_vep_columns = vep_columns - handled_vep_columns
            if unused_vep_columns:
                logging.warning("Unhandled VEP CSQ columns (maybe add to VEP_NOT_COPIED_FIELDS or disable in VEP pipeline?) :")
                logging.warning(", ".join(sorted(unused_vep_columns)))

            missing_expected_vep_columns = handled_vep_columns - vep_columns
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
        self.aloft_columns = self._has_all_aloft_columns(self.all_variant_columns)

    @staticmethod
    def _has_all_aloft_columns(columns) -> bool:
        return all([ac in columns for ac in BulkVEPVCFAnnotationInserter.ALOFT_COLUMNS])

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

        if self.aloft_columns and self._has_all_aloft_columns(raw_db_data):
            self._pick_aloft_values(raw_db_data)

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

    @staticmethod
    def _pick_aloft_values(raw_db_data: dict):
        """ ALoFT produces values for multiple transcripts, ie raw values before formatting:
                aloft_prob_tolerant     .&0.0516&.&.&.&
                aloft_prob_recessive    .&0.81255&.&.&.&
                aloft_prob_dominant     .&0.13585&.&.&.&
                aloft_pred              .&Recessive&.&.&.&
                aloft_high_confidence        .&High&.&.&.&
                aloft_ensembl_transcript    ENST00000565905&ENST00000361627&ENST00000567348&ENST00000563864&ENST00000543522

            So when we pull a consistent column out of all fields
        """

        aloft_cols = [raw_db_data[a].split("&") for a in BulkVEPVCFAnnotationInserter.ALOFT_COLUMNS]
        aloft_options = []
        for aloft_option in zip(*aloft_cols):
            aloft_options.append(dict(zip(BulkVEPVCFAnnotationInserter.ALOFT_COLUMNS, aloft_option)))

        # Pick High confidence, then Recessive over Dominant
        PREFERENCES = {
            "aloft_high_confidence": ["High", "Low", "."],
            "aloft_pred": ["Recessive", "Dominant", "Tolerant", "."],
        }

        def aloft_preferences(val):
            return [prefs.index(val[k]) for k, prefs in PREFERENCES.items()]

        aloft_options = sorted(aloft_options, key=aloft_preferences)
        best_aloft = aloft_options[0]
        # ensembl_transcript will always be populated, but we don't really want it if the other ALoFT fields are empty
        if best_aloft["aloft_pred"] == ".":
            best_aloft["aloft_ensembl_transcript"] = "."
        raw_db_data.update(best_aloft)

    def get_gene_id(self, vep_transcript_data):
        gene_id = None
        vep_gene = vep_transcript_data[VEPColumns.GENE]
        if vep_gene:
            if vep_gene in self.gene_identifiers:
                gene_id = vep_gene
            else:
                logging.warning(f"Could not find gene_id: '{vep_gene}'")
        return gene_id

    def _get_transcript_id_and_version(self, vep_transcript_data) -> tuple[Optional[str], Optional[str]]:
        transcript_accession = None
        transcript_id = None
        transcript_version_id = None
        vep_feature_type = vep_transcript_data[VEPColumns.FEATURE_TYPE]
        if vep_feature_type == "Transcript":
            if transcript_accession := vep_transcript_data[VEPColumns.FEATURE]:
                # We pass in --transcript_version so Ensembl IDs will have version
                t_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
                transcript_versions = self.transcript_versions_by_id.get(t_id)
                if transcript_versions:
                    transcript_id = t_id  # Know it's valid to link
                    transcript_version_id = transcript_versions.get(version)
                    if transcript_version_id is None:
                        logging.warning(f"Have transcript '{transcript_id}' but no version: '{version}'")
                else:
                    logging.warning(f"Could not find transcript: '{t_id}'")
        return transcript_id, transcript_version_id

    def add_calculated_transcript_columns(self, variant_coordinate: Optional[VariantCoordinate], transcript_data):
        """ variant_coordinate - will only be set for symbolics """
        self._add_calculated_maxentscan(transcript_data)
        self._add_hgvs(variant_coordinate, transcript_data)

    def add_calculated_variant_annotation_columns(self, transcript_data):
        self._add_calculated_num_predictions(transcript_data)
        self._add_hemi_count(transcript_data)

    def _add_calculated_num_predictions(self, transcript_data):
        num_pathogenic = 0
        num_not_pathogenic = 0
        for field, path_func in self.prediction_pathogenic_funcs.items():
            value = transcript_data.get(field)
            if value is not None:
                if path_func(value):
                    num_pathogenic += 1
                else:
                    num_not_pathogenic += 1
        transcript_data["predictions_num_pathogenic"] = num_pathogenic
        transcript_data["predictions_num_benign"] = num_not_pathogenic

    def _add_hemi_count(self, transcript_data):
        """ gnomad_non_par=True means not on pseudoautosomal region on chrX, so XY count = hemizygous """
        if transcript_data.get("gnomad_non_par"):
            transcript_data["gnomad_hemi_count"] = transcript_data.get("gnomad_xy_ac")

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

    def _add_hgvs(self, variant_coordinate: Optional[VariantCoordinate], transcript_data: dict):
        if variant_coordinate is None:
            return

        if abs(variant_coordinate.svlen) > settings.HGVS_MAX_ALT_LENGTH:
            if not transcript_data.get(VEPColumns.PICK, False):
                logging.warning(f"Skipping calculating HGVS for {variant_coordinate} NON-representative transcript "
                                f"because it exceeds {settings.HGVS_MAX_ALT_LENGTH}")
                return

        transcript_id = transcript_data.get("transcript_id")
        version = transcript_data.get("version_id")
        if transcript_id and version:
            transcript_accession = TranscriptVersion.get_accession(transcript_id, version)
            try:
                hgvs_c = self.hgvs_matcher.variant_coordinate_to_hgvs_variant(variant_coordinate, transcript_accession)
                transcript_data['hgvs_c'] = hgvs_c
                # TODO: Protein?? hgvs_p
            except (ValueError, HGVSException):
                pass

    def process_entry(self, v):
        if len(self.variant_transcript_annotation_list) >= settings.SQL_BATCH_INSERT_SIZE:
            self.bulk_insert()

        csq = v.INFO.get("CSQ")
        if csq is None:  # Legacy data error, VEP failed to annotate as something like C>C
            logging.warning(f"Skipped {v} as no CSQ")
            return

        if svlen := v.INFO.get("SVLEN"):
            variant_coordinate = VariantCoordinate(chrom=v.CHROM, position=v.POS, ref=v.REF, alt=v.ALT[0], svlen=svlen)
            # Do now so we only retrieve sequences once
            variant_coordinate = variant_coordinate.as_external_explicit(self.annotation_run.genome_build)
        else:
            variant_coordinate = None

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
                transcript_id, transcript_version_id = self._get_transcript_id_and_version(vep_transcript_data)
                transcript_data["transcript_id"] = transcript_id
                transcript_data["transcript_version_id"] = transcript_version_id
                if symbol := transcript_data.get("symbol"):
                    overlapping_symbols.add(symbol)
                if gene_id:
                    overlapping_gene_ids.add(gene_id)
                self.add_calculated_transcript_columns(variant_coordinate, transcript_data)
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
            VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION: (VariantAnnotation, self.variant_annotation_list),
            VariantAnnotationVersion.TRANSCRIPT_ANNOTATION: (VariantTranscriptAnnotation, self.variant_transcript_annotation_list),
            VariantAnnotationVersion.VARIANT_GENE_OVERLAP: (VariantGeneOverlap, self.variant_gene_overlap_list),
        }

        for base_table_name, (klass, annotations_list) in ANNOTATION_TYPE.items():
            # In unit testing, cursor.copy_from causes:
            # psycopg2.errors.NotNullViolation: null value in column "id" of relation ... violates not-null constraint
            if settings.UNIT_TEST:
                # Goes into main table not partition
                records = [klass(**kwargs) for kwargs in annotations_list]
                klass.objects.bulk_create(records, batch_size=2000)
            else:
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

    @cached_property
    def gene_identifiers(self):
        """ A set of gene identifiers to check whether in DB or not """
        vav = self.annotation_run.variant_annotation_version
        gv_qs = GeneVersion.objects.filter(gene__annotation_consortium=vav.annotation_consortium,
                                           genome_build=vav.genome_build)
        return set(gv_qs.distinct().values_list("gene", flat=True))

    @cached_property
    def transcript_versions_by_id(self):
        vav = self.annotation_run.variant_annotation_version
        return TranscriptVersion.transcript_versions_by_id(vav.genome_build, vav.annotation_consortium)


def empty_to_none(it):
    return [v if v not in EMPTY_VALUES else None
            for v in it]


# Field formatters
def gnomad_filtered_func(raw_value):
    """ We use FILTER in Gnomad3 (GRCh38 only) - need to convert back to bool
        In the combined exomes/genomes (gnomad2, gnomad4) we use gnomad_filtered=1
        So don't need to format this etc
    """
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

def get_format_alphamissense_class_func():
    """ GRCh37 has 'benign' while GRCh38 has 'likely_benign'
        @see https://github.com/Ensembl/VEP_plugins/issues/668
    """
    cff = get_choice_formatter_func(AlphaMissensePrediction.choices)
    def _format_alphamissense_class(alphamissense_class):
        if alphamissense_class == "benign":
            alphamissense_class = "likely_benign"
        return cff(alphamissense_class)
    return _format_alphamissense_class

def get_extract_existing_variation(prefix):
    def format_vep_existing_variation(vep_existing_variation):
        ev_list = vep_existing_variation.split(VEP_SEPARATOR)
        cosmic_ids = [ev for ev in ev_list if ev.startswith(prefix)]
        return VEP_SEPARATOR.join(sorted(set(cosmic_ids)))

    return format_vep_existing_variation


def format_vep_somatic(raw_value):
    return "1" in raw_value


def get_choice_formatter_func(choices, empty_values=None):
    lookup = invert_dict(dict(choices))

    def format_choice(raw_value):
        if empty_values is not None:
            if raw_value in empty_values:
                return None
        return lookup[raw_value]

    return format_choice


def get_clean_and_pick_single_value_func(pick_single_value_func, cast_func=None, empty_values=None):
    """ Returns a function to clean and pick single value.
        casting is performed before calling pick_single_value_func so you can call min/max """

    if empty_values is None:
        empty_values = EMPTY_VALUES

    def _clean_and_pick_single_value_func(raw_value):
        it = (tm for tm in raw_value.split(VEP_SEPARATOR) if tm != '')
        # Handle '.'
        if cast_func:
            values = [cast_func(v) for v in it if v not in empty_values]
        else:
            values = [v for v in it if v not in empty_values]
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


def get_format_empty_as_none(empty_values: set):
    def format_empty_as_none(val):
        if val in empty_values:
            val = None
        return val
    return format_empty_as_none


def format_nmd_escaping_variant(value) -> bool:
    return value == "NMD_escaping_variant"


def format_aloft_high_confidence(value) -> bool:
    high_confidence = None
    if value == "High":
        high_confidence = True
    elif value == "Low":
        high_confidence = False
    return high_confidence


def format_canonical(value) -> bool:
    return value == "YES"
