import logging
import operator
import shutil
import time
from collections import defaultdict, Counter
from functools import cached_property
from typing import Optional, Iterable, TypeAlias

import intervaltree
from django.conf import settings

from annotation import vep_columns as vep_columns_registry
from annotation.models.models import VariantAnnotation, \
    VariantTranscriptAnnotation, VariantAnnotationVersion, VariantGeneOverlap, AnnotationRun
from annotation.models.models_enums import VariantAnnotationPipelineType, VEPCustom
from annotation.refseq_ensembl_resolver import DBNSFPGeneResolver
from annotation.vcf_files.vcf_types import VCFVariant
from annotation.vep_annotation import VEPConfig
from annotation.vep_columns import VEPColumnDef
from annotation.vep_field_formatters import VEP_SEPARATOR, EMPTY_VALUES
from genes.hgvs import HGVSMatcher
from genes.models import TranscriptVersion, GeneVersion
from genes.models_enums import AnnotationConsortium
from library.django_utils import get_model_fields
from library.django_utils.django_file_utils import get_import_processing_filename, get_import_processing_dir
from library.genomics import overlap_fraction, Range, parse_gnomad_coord
from library.log_utils import log_traceback
from library.utils import split_dict_multi_values
from snpdb.models import VariantCoordinate
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv

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


TranscriptData: TypeAlias = dict


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
        "hgvs_g",
        "spliceai_max_ds",
    ]
    DB_IGNORED_COLUMNS = ["id", "transcript", "MaveDB_nt", "MaveDB_pro"]
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

    # dbNSFP source fields that come as &-separated arrays parallel to Ensembl_transcriptid
    # (columns_version >= 4 only — earlier versions only consumed dbNSFP variant-level fields).
    # _pick_dbnsfp_per_transcript_values rewrites these to the single value matching the picked
    # VEP transcript when the resolver finds a mapping; otherwise their formatters collapse
    # the array to a representative value.
    DBNSFP_PER_TRANSCRIPT_SOURCE_FIELDS = (
        "AlphaMissense_pred", "AlphaMissense_score",
        "MPC_score",
        "MetaRNN_pred", "MetaRNN_score",
        "MutPred2_score", "MutPred2_top5_mechanisms",
        "REVEL_score",
        "VARITY_ER_score", "VARITY_R_score",
        "VEST4_score",
    )

    def __init__(self,
                 annotation_run: AnnotationRun,
                 infos: Optional[dict] = None,
                 insert_variants: bool = True,
                 validate_columns: bool = True):
        self.annotation_run = annotation_run
        self.genome_build = self.annotation_run.variant_annotation_version.genome_build
        self.vep_config = VEPConfig(self.genome_build)
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

        cvf_list = vep_columns_registry.filter_for(
            vep_config=self.vep_config,
            pipeline_type=self.annotation_run.pipeline_type,
        )

        self._setup_vep_fields_and_db_columns(validate_columns, cvf_list)
        self.hgvs_matcher = HGVSMatcher(annotation_run.genome_build,
                                        # We only want exact transcript version for annotation
                                        allow_alternative_transcript_version=False)

        sv_overlap_processor = None
        sv_gene_overlap_resolver = None
        if self.annotation_run.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT:
            sv_overlap_processor = SVOverlapProcessor(cvf_list)
            sv_gene_overlap_resolver = SVGeneOverlapResolver(self.annotation_run.variant_annotation_version)
        self.sv_overlap_processor = sv_overlap_processor
        self.sv_gene_overlap_resolver = sv_gene_overlap_resolver
        self._generated_hgvs_c = Counter()

        # Resolver for picking per-transcript dbNSFP values (RefSeq <-> Ensembl translation).
        # Only relevant for columns_version >= 4 (which introduced per-transcript dbNSFP
        # score fields). Persist the resolver name on the VAV so we have provenance for
        # which strategy populated the scores.
        self.transcript_resolver = None
        if self.vep_config.columns_version >= 4:
            self.transcript_resolver = DBNSFPGeneResolver()
            vav = self.annotation_run.variant_annotation_version
            if vav.transcript_resolver != self.transcript_resolver.name:
                vav.transcript_resolver = self.transcript_resolver.name
                vav.save(update_fields=["transcript_resolver"])

    @property
    def description(self) -> str:
        return f"{self.genome_build}/{self.vep_config.columns_version}/{self.annotation_run.get_pipeline_type_display()}"

    @staticmethod
    def _get_vep_columns_from_csq(infos: Optional[dict]):
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

    def _add_vep_field_handlers(self, cvf_list):
        # Per-column value cleaners now live on the VEPColumnDefs (annotation.vep_columns),
        # so build the destination-column -> formatter map straight from the active defs.
        # cvf_list is already filtered through vep_config (build / columns_version / vep_version /
        # gnomad minor), so the correct formatter for this version falls out automatically -
        # including gnomad_filtered (FILTER-sourced defs only) and the columns_version >= 4
        # dbNSFP score / pred fields.
        self.field_formatters = {
            vgc: cvf.formatter
            for cvf in cvf_list
            for vgc in cvf.variant_grid_columns
            if cvf.formatter is not None
        }

        self.source_field_to_columns = defaultdict(set)
        self.ignored_vep_fields = self.VEP_NOT_COPIED_FIELDS.copy()

        # cvf_list is already filtered through vep_config so unconfigured customs are dropped.
        # Sort to have consistent VCF headers (case-insensitive to match postgres `ORDER BY source_field`)
        for cvf in sorted(cvf_list, key=lambda c: (c.source_field or "").lower()):
            for vgc_id in cvf.variant_grid_columns:
                self.source_field_to_columns[cvf.vep_info_field].add(vgc_id)

        vav = self.annotation_run.variant_annotation_version
        if vav.columns_version >= 4:
            self.prediction_pathogenic_funcs = vav.get_raw_score_pathogenic_prediction_funcs()
        else:
            self.prediction_pathogenic_funcs = vav.get_rankscore_pathogenic_prediction_funcs()

    def _setup_vep_fields_and_db_columns(self, validate_columns: bool, cvf_list: tuple[VEPColumnDef, ...]):
        self._add_vep_field_handlers(cvf_list)

        ignore_columns = set(self.DB_FIXED_COLUMNS +
                             self.DB_MANUALLY_POPULATED_COLUMNS +
                             self.DB_IGNORED_COLUMNS)
        if self.annotation_run.annotation_consortium == AnnotationConsortium.REFSEQ:
            ignore_columns.update(self.VEP_NOT_COPIED_REFSEQ_ONLY)

        # Find the ones that don't apply to this version, and exclude them
        in_scope = {vgc for c in cvf_list for vgc in c.variant_grid_columns}
        all_known = vep_columns_registry.all_variant_grid_column_ids()
        ignore_columns.update(all_known - in_scope)

        for c in list(ignore_columns):
            # FIXME, why do we convery ignore_clumns to a list here? it's a set which would work fine
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
                logging.warning(f"{self.description}: Unhandled VEP CSQ columns (maybe add to VEP_NOT_COPIED_FIELDS or disable in VEP pipeline?) :")
                logging.warning(", ".join(sorted(unused_vep_columns)))

            missing_expected_vep_columns = handled_vep_columns - vep_columns
            if missing_expected_vep_columns:
                missing = ", ".join(sorted(missing_expected_vep_columns))
                msg = f"{self.description}: Fields missing from VEP annotated vcf CSQ columns: {missing}"
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
        # Manual columns are emitted first; remove any overlap so the COPY header has no duplicates.
        columns = set(columns) - set(manual_columns)
        return self.DB_FIXED_COLUMNS + manual_columns + sorted(columns)

    def vep_to_db_dict(self, vep_transcript_data: TranscriptData, model_columns: Iterable[str]) -> dict:
        # logging.debug("vep_to_db_dict:")
        # logging.debug(vep_transcript_data)

        if self.transcript_resolver is not None:
            self._pick_dbnsfp_per_transcript_values(vep_transcript_data)

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

        if self.sv_overlap_processor:
            self.sv_overlap_processor.process(raw_db_data)

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

    def _pick_dbnsfp_per_transcript_values(self, vep_transcript_data: TranscriptData):
        """ dbNSFP returns several score/pred fields as &-separated arrays parallel to
            Ensembl_transcriptid. Look up the Ensembl id matching the current VEP transcript
            via self.transcript_resolver, find its index in the array, and rewrite the
            per-transcript fields to that single value. Leaves arrays untouched if no
            mapping is available — the field formatters will then collapse them to a
            representative value (max float / most-damaging pred).

            Note: dbNSFP's Ensembl_transcriptid column is unversioned (e.g. ENST00000262340),
            so we strip the .N suffix from the VEP Feature on both sides. """
        ensembl_ids_raw = vep_transcript_data.get("Ensembl_transcriptid")
        if not ensembl_ids_raw or VEP_SEPARATOR not in ensembl_ids_raw:
            return  # Single value (or absent) — no array to pick from

        feature = vep_transcript_data.get("Feature")
        if not feature:
            return
        feature_unversioned = feature.split(".", 1)[0]

        if self.annotation_run.annotation_consortium == AnnotationConsortium.REFSEQ:
            symbol = vep_transcript_data.get("SYMBOL")
            if not symbol:
                return
            target_ensembl = self.transcript_resolver.refseq_to_ensembl(symbol, feature_unversioned)
        else:
            target_ensembl = feature_unversioned  # Ensembl pipeline: VEP Feature is already ENST

        if not target_ensembl:
            return

        ensembl_ids = ensembl_ids_raw.split(VEP_SEPARATOR)
        try:
            idx = ensembl_ids.index(target_ensembl)
        except ValueError:
            return  # Target ENST not in the array — fall through to formatter collapse

        for source_field in self.DBNSFP_PER_TRANSCRIPT_SOURCE_FIELDS:
            value = vep_transcript_data.get(source_field)
            if value and VEP_SEPARATOR in value:
                parts = value.split(VEP_SEPARATOR)
                if idx < len(parts):
                    vep_transcript_data[source_field] = parts[idx]

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

    def get_gene_id(self, vep_transcript_data: TranscriptData):
        gene_id = None
        vep_gene = vep_transcript_data[VEPColumns.GENE]
        if vep_gene:
            if vep_gene in self.gene_identifiers:
                gene_id = vep_gene
            else:
                logging.warning(f"Could not find gene_id: '{vep_gene}'")
        return gene_id

    def _get_transcript_accession(self, vep_transcript_data: TranscriptData) -> Optional[str]:
        # We pass in --transcript_version so Ensembl IDs will have version
        vep_feature_type = vep_transcript_data[VEPColumns.FEATURE_TYPE]
        transcript_accession = None
        if vep_feature_type == "Transcript":
            transcript_accession = vep_transcript_data[VEPColumns.FEATURE]
        return transcript_accession

    def _get_transcript_id_and_transcript_version_id(self, transcript_accession: Optional[str]) -> tuple[Optional[str], Optional[str]]:
        """ Returns Transcript.pk and TranscriptVersion.pk for linking records """
        transcript_id: Optional[str] = None
        transcript_version_id: Optional[str] = None
        if transcript_accession:
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

    def _fix_multiple_values(self, transcript_data: TranscriptData):
        # T2T gnomad liftover has dupes, sometimes we get multiple entries in which case all will be dupes
        # We are going to pick the HIGHEST, as if there was a liftover error this will cause false positives after
        # gnomAD filter rather than false negatives (potentially throwing away real rare disease causing variants)
        if gnomad_af := transcript_data.get("gnomad_af"):
            if "&" in gnomad_af:
                gnomad = {k: v for k, v in transcript_data.items() if k.startswith("gnomad")}
                gnomad_list = split_dict_multi_values(gnomad, sep='&')
                highest_af = sorted(gnomad_list, key=operator.itemgetter("gnomad_af"), reverse=True)[0]
                # Copy back overwriting old multi-value data with single values
                for k, v in highest_af.items():
                    if v == '.':
                        v = None
                    transcript_data[k] = v

    def add_calculated_transcript_columns(self, variant_coordinate: Optional[VariantCoordinate],
                                          transcript_accession: Optional[str],
                                          transcript_data: TranscriptData):
        """ variant_coordinate - will only be set for symbolics """

        if self.annotation_run.pipeline_type == VariantAnnotationPipelineType.STANDARD:
            self._add_calculated_maxentscan(transcript_data)
        if self.annotation_run.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT:
            self._add_hgvs_c(variant_coordinate, transcript_accession, transcript_data)

    def add_calculated_variant_annotation_columns(self, variant_coordinate: VariantCoordinate, transcript_data: TranscriptData):
        self._add_hemi_count(transcript_data)
        self._add_hgvs_g(variant_coordinate, transcript_data)
        self._add_calculated_num_predictions(transcript_data)
        self._add_spliceai_max_ds(transcript_data)
        if self.annotation_run.pipeline_type == VariantAnnotationPipelineType.STRUCTURAL_VARIANT:
            self._calculate_gnomad_sv_overlap_percentage(variant_coordinate, transcript_data)

    @staticmethod
    def _add_spliceai_max_ds(transcript_data: TranscriptData):
        ds_values = []
        for key in ("spliceai_pred_ds_ag", "spliceai_pred_ds_al",
                    "spliceai_pred_ds_dg", "spliceai_pred_ds_dl"):
            v = transcript_data.get(key)
            if v is not None:
                ds_values.append(float(v))
        if ds_values:
            transcript_data["spliceai_max_ds"] = max(ds_values)

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

    def _add_hemi_count(self, transcript_data: TranscriptData):
        """ gnomad_non_par=True means not on pseudoautosomal region on chrX, so XY count = hemizygous """
        if transcript_data.get("gnomad_non_par"):
            transcript_data["gnomad_hemi_count"] = transcript_data.get("gnomad_xy_ac")

    def _calculate_gnomad_sv_overlap_percentage(self, variant_coordinate: VariantCoordinate, transcript_data: TranscriptData):
        SV_PC_FIELD = "gnomad_sv_overlap_percent"

        # These are currently coordinates that look like eg: chr17:42853981-43253980&chr17:43048298-43898190
        if orig_gnomad_sv_overlap_percent := transcript_data.get(SV_PC_FIELD):
            gnomad_sv_overlap_percent_list = []

            for coord in orig_gnomad_sv_overlap_percent.split(VEP_SEPARATOR):
                chrom, start, end = parse_gnomad_coord(coord)
                of = overlap_fraction(Range(variant_coordinate.position, variant_coordinate.end),
                                      Range(start, end))
                gnomad_sv_overlap_percent_list.append(int(of * 100))

            transcript_data[SV_PC_FIELD] = VEP_SEPARATOR.join([str(p) for p in gnomad_sv_overlap_percent_list])

    @staticmethod
    def _add_calculated_maxentscan(transcript_data: TranscriptData):
        maxentscan_diff = transcript_data.get("maxentscan_diff")
        maxentscan_ref = transcript_data.get("maxentscan_ref")
        if maxentscan_diff and maxentscan_ref:
            maxentscan_diff = float(maxentscan_diff)
            maxentscan_ref = float(maxentscan_ref)
            if maxentscan_ref:
                transcript_data["maxentscan_percent_diff_ref"] = 100 * maxentscan_diff / maxentscan_ref

    @staticmethod
    def _merge_cosmic_ids(transcript_data: TranscriptData, custom_vcf_cosmic_id: str):
        """ VEP ships w/COSMIC so we always try and pull it out of the existing variation field
            We also support a custom vcf of COSMIC (for cosmic_count) - this is often more recent """
        cosmic_ids = {c for c in custom_vcf_cosmic_id.split(VEP_SEPARATOR) if c}
        if existing_cosmic_id := transcript_data.get("cosmic_id"):
            cosmic_ids.update(existing_cosmic_id.split(VEP_SEPARATOR))
        transcript_data["cosmic_id"] = VEP_SEPARATOR.join(sorted(cosmic_ids))

    def _add_hgvs_c(self, variant_coordinate: Optional[VariantCoordinate], transcript_accession: Optional[str],
                    transcript_data: dict):
        # VEP will have already done for non-symbolics, may do them in future version
        if transcript_data.get('hgvs_c'):
            return

        if transcript_data.get(VEPColumns.PICK):
            max_length = settings.HGVS_MAX_SEQUENCE_LENGTH_REPRESENTATIVE_TRANSCRIPT
        else:
            max_length = settings.HGVS_MAX_SEQUENCE_LENGTH

        # Only calculate very long HGVS for representative transcripts
        if variant_coordinate.max_sequence_length > max_length:
            transcript_data['hgvs_c'] = VariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE
            self._generated_hgvs_c["too_long"] += 1
            return

        if transcript_accession:

            try:
                hgvs_c = self.hgvs_matcher.variant_coordinate_to_hgvs_variant(variant_coordinate, transcript_accession)
                self._generated_hgvs_c["OK"] += 1
                # TODO: Protein?? hgvs_p
            except Exception as e:
                logging.debug("Error calculating c.HGVS for '%s'/'%s': %s",
                              variant_coordinate.format_short(), transcript_accession, e)
                hgvs_c = None  # For c.HGVS - it's ok to be blank
                self._generated_hgvs_c["error"] += 1

            transcript_data['hgvs_c'] = hgvs_c

    def _add_hgvs_g(self, variant_coordinate: Optional[VariantCoordinate], transcript_data: TranscriptData):
        # VEP110 has a bug with --hgvsg but we hope to introduce in VEP111+
        if transcript_data.get('hgvs_g'):
            return

        max_length = settings.HGVS_MAX_SEQUENCE_LENGTH_REPRESENTATIVE_TRANSCRIPT  # VariantAnnotation
        if variant_coordinate.max_sequence_length > max_length:
            hgvs_g = VariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE
        else:
            try:
                hgvs_g = self.hgvs_matcher.variant_coordinate_to_g_hgvs(variant_coordinate)
            except Exception as e:
                logging.warning("Error calculating g.HGVS for '%s': %s", variant_coordinate.format_short(), e)
                hgvs_g = VariantAnnotation.SV_HGVS_ERROR_MESSAGE

        transcript_data['hgvs_g'] = hgvs_g

    def process_entry(self, v: VCFVariant):
        if len(self.variant_transcript_annotation_list) >= settings.SQL_BATCH_INSERT_SIZE:
            self.bulk_insert()

        csq = v.INFO.get("CSQ")
        if csq is None:  # Legacy data error, VEP failed to annotate as something like C>C
            logging.warning(f"Skipped {v} as no CSQ")
            return

        svlen = v.INFO.get("SVLEN")
        variant_coordinate = VariantCoordinate(chrom=v.CHROM, position=v.POS, ref=v.REF, alt=v.ALT[0], svlen=svlen)
        # Do now so we only retrieve sequences once. <CNV>/<INS> can't be expanded to explicit
        # ref/alt - leave them symbolic (downstream HGVS/SV-overlap handle the symbolic form)
        if variant_coordinate.can_be_made_explicit:
            variant_coordinate = variant_coordinate.as_external_explicit(self.annotation_run.genome_build)

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
                transcript_accession = self._get_transcript_accession(vep_transcript_data)
                gene_id = self.get_gene_id(vep_transcript_data)
                transcript_data = self.vep_to_db_dict(vep_transcript_data, self.transcript_columns)
                transcript_data.update(self.constant_data)
                transcript_data["variant_id"] = variant_id
                transcript_data["gene_id"] = gene_id
                transcript_id, transcript_version_id = self._get_transcript_id_and_transcript_version_id(transcript_accession)
                transcript_data["transcript_id"] = transcript_id
                transcript_data["transcript_version_id"] = transcript_version_id
                if symbol := transcript_data.get("symbol"):
                    overlapping_symbols.add(symbol)
                if gene_id:
                    overlapping_gene_ids.add(gene_id)
                self.add_calculated_transcript_columns(variant_coordinate, transcript_accession, transcript_data)
                self.variant_transcript_annotation_list.append(transcript_data)

                representative_transcript = vep_transcript_data.get(VEPColumns.PICK, False)
                if representative_transcript:
                    variant_data = self.vep_to_db_dict(vep_transcript_data, self.variant_only_columns)
                    variant_data.update(transcript_data)
                    self._fix_multiple_values(variant_data)
                    self.add_calculated_variant_annotation_columns(variant_coordinate, variant_data)
                    # If we're using custom COSMIC vcf, merge with those from VEP existing variation
                    if custom_vcf_cosmic_ids := vep_transcript_data.get("COSMIC"):
                        self._merge_cosmic_ids(variant_data, custom_vcf_cosmic_ids)

            # Handle variant data at the end so we can store overlapping symbols
            if variant_data:
                if overlapping_symbols:
                    variant_data["overlapping_symbols"] = ",".join(sorted(overlapping_symbols))

                # logging.debug(f"variant_data:")
                # logging.debug(variant_data)
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

    def add_sv_gene_overlaps(self, variant_id: int, variant_coordinate: VariantCoordinate, variant_data: dict):
        """ For long SVs that VEP skipped (TOO_LONG): resolve gene overlaps locally
            from the gene_annotation_release transcripts and populate
            overlapping_symbols + queue VariantGeneOverlap rows.
            Leaves gene_id/transcript_id on the variant row null (a long SV may overlap many genes). """
        if self.sv_gene_overlap_resolver is None:
            return
        symbols, gene_ids = self.sv_gene_overlap_resolver.get_overlaps(variant_coordinate)
        if symbols:
            variant_data["overlapping_symbols"] = ",".join(sorted(symbols))
        for gene_id in gene_ids:
            overlap = {"variant_id": variant_id, "gene_id": gene_id}
            overlap.update(self.constant_data)
            self.variant_gene_overlap_list.append(overlap)

    def finish(self):
        self.bulk_insert()
        if self._generated_hgvs_c:
            logging.info("Generated HGVS C results: %s" % self._generated_hgvs_c)

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
        # ignore_errors so a missing dir (eg cleaned-up retry) doesn't blow up - we just want it gone
        shutil.rmtree(import_processing_dir, ignore_errors=True)

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


class SVOverlapProcessor:
    """
        We use --custom (twice) rather than StructuralVariantOverlap
        due to some issues: https://github.com/Ensembl/VEP_plugins/issues/710
        This requires a bit of post-processing
    """
    def __init__(self, cvf_list: tuple[VEPColumnDef, ...]):
        sv_customs = {VEPCustom.GNOMAD_SV, VEPCustom.GNOMAD_SV_NAME}
        self.sv_fields = {
            vgc for c in cvf_list if c.vep_custom in sv_customs
            for vgc in c.variant_grid_columns
        }

    @staticmethod
    def _get_required_substrings(variant_class: str) -> list[str]:
        required_substrings = []  # Empty means all will go through
        if settings.ANNOTATION_VEP_SV_OVERLAP_SAME_TYPE:
            # gnomad_sv_overlap_name looks like: 'gnomAD-SV_v2.1_INV_17_731&gnomAD-SV_v2.1_DEL_17_160435'
            EXPECTED_NAME = {
                'deletion': ["DEL"],
                'duplication': ['DUP'],
                'inversion': ['INV'],
                'indel': ["DEL", "INS"]
            }
            required_substrings = EXPECTED_NAME.get(variant_class, [])
        return required_substrings

    def process(self, raw_db_data: dict):
        # This one should always there if any overlaps
        if 'gnomad_sv_overlap_name' not in raw_db_data:
            return

        # The StructuralVariantOverlap fields are joined via '&'
        # Get it into a nice structure then process it
        sv_records = VariantAnnotation.vep_multi_fields_to_list_of_dicts(raw_db_data, self.sv_fields)
        if required_substrings := self._get_required_substrings(raw_db_data["variant_class"]):
            filtered_sv_records = []
            for r in sv_records:
                for substring in required_substrings:
                    if substring in r['gnomad_sv_overlap_name']:
                        filtered_sv_records.append(r)
                        break
        else:
            # No filtering required - all go through
            filtered_sv_records = sv_records

        if filtered_sv_records:
            chosen_record = self._pick_record(filtered_sv_records)
            # Update if not special multi-field
            update_fields = {}
            for k, v in chosen_record.items():
                if k not in VariantAnnotation.GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS:
                    update_fields[k] = v

            # Need to re-build the multi-files from what is filtered
            for f in VariantAnnotation.GNOMAD_SV_OVERLAP_MULTI_VALUE_FIELDS:
                update_fields[f] = VEP_SEPARATOR.join([r[f] for r in filtered_sv_records])
        else:
            # Nothing left after filtering - need to blank out all our values
            update_fields = {f: None for f in self.sv_fields}

        raw_db_data.update(update_fields)

    def _pick_record(self, filtered_sv_records):
        if len(filtered_sv_records) == 1:
            return filtered_sv_records[0]

        chosen_record = None
        if settings.ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD == "greatest_overlap":
            raise NotImplementedError("greatest_overlap")
        elif settings.ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD == "lowest_af":
            for record in filtered_sv_records:
                if chosen_record:
                    if record["gnomad_sv_overlap_af"] > chosen_record["gnomad_sv_overlap_af"]:
                        continue
                chosen_record = record
        elif settings.ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD == "exact_or_lowest_af":
            raise NotImplementedError("exact_or_lowest_af")
        else:
            raise ValueError(f"Unknown value for {settings.ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD=}")

        return chosen_record


class SVGeneOverlapResolver:
    """ Resolves gene overlaps for long SVs that VEP skipped due to TOO_LONG.

        Builds an in-memory per-contig IntervalTree of TranscriptVersions in the
        VariantAnnotationVersion's gene_annotation_release. For each variant, returns
        the set of overlapping (symbol, gene_id) pairs.
    """

    def __init__(self, variant_annotation_version: VariantAnnotationVersion):
        self.variant_annotation_version = variant_annotation_version
        gene_annotation_release = variant_annotation_version.gene_annotation_release
        self._trees: dict[str, intervaltree.IntervalTree] = defaultdict(intervaltree.IntervalTree)

        if gene_annotation_release is None:
            logging.warning("SVGeneOverlapResolver: no gene_annotation_release on %s", variant_annotation_version)
            return

        start_time = time.monotonic()
        tv_qs = TranscriptVersion.objects.filter(
            releasetranscriptversion__release=gene_annotation_release,
        ).select_related("gene_version__gene_symbol", "contig")

        count = 0
        for tv in tv_qs:
            try:
                start = tv.start
                end = tv.end
            except (KeyError, IndexError):
                continue
            if end <= start:
                # intervaltree treats zero-length intervals as empty
                end = start + 1
            symbol = tv.gene_version.gene_symbol_id
            gene_id = tv.gene_version.gene_id
            self._trees[tv.contig.name].addi(start, end, (gene_id, symbol))
            count += 1

        elapsed = time.monotonic() - start_time
        logging.info(
            "SVGeneOverlapResolver: built %d intervals across %d contigs for %s in %.2fs",
            count, len(self._trees), variant_annotation_version, elapsed,
        )

    def get_overlaps(self, variant_coordinate: VariantCoordinate) -> tuple[set[str], set[str]]:
        """ Returns (overlapping_symbols, overlapping_gene_ids) for the given variant. """
        symbols: set[str] = set()
        gene_ids: set[str] = set()
        tree = self._trees.get(variant_coordinate.chrom)
        if tree is None:
            return symbols, gene_ids

        start = variant_coordinate.position
        end = variant_coordinate.end
        if end <= start:
            end = start + 1

        for interval in tree.overlap(start, end):
            gene_id, symbol = interval.data
            if gene_id is not None:
                gene_ids.add(gene_id)
            if symbol is not None:
                symbols.add(symbol)
        return symbols, gene_ids
