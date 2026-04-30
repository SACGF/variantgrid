import logging
import os
import re
import uuid
from shlex import shlex
from typing import Optional

from django.conf import settings

from annotation import vep_columns
from annotation.fake_annotation import get_fake_vep_version
from annotation.models.models_enums import VEPPlugin, VEPCustom, VariantAnnotationPipelineType
from annotation.vep_columns import VEPColumnDef
from annotation.vep_config import VEPConfig, parse_gnomad_version_from_filename
from genes.models_enums import AnnotationConsortium
from library.utils import execute_cmd
from library.utils.file_utils import get_extension_without_gzip, mk_path_for_file, open_handle_gzip
from snpdb.models.models_genome import GenomeBuild


class VEPVersionMismatchError(ValueError):
    pass


def _get_dbnsfp_plugin_command(genome_build: GenomeBuild, vc: VEPConfig):
    """ Build from VEPColumnDef.source_field where vep_plugin = DBNSFP """

    dbnsfp_data_path = vc["dbnsfp"]
    fields = vep_columns.source_fields_for(
        vep_config=vc,
        vep_plugin=VEPPlugin.DBNSFP,
    )
    return f"dbNSFP,{dbnsfp_data_path},{','.join(fields)}"


def _get_custom_params_list(cvf_list: list[VEPColumnDef], prefix, data_path) -> list:
    """ All our deployments are VEP >= 110 so we can use key/value pairs """
    int_vep_version = int(settings.ANNOTATION_VEP_VERSION)
    if int_vep_version < 110:
        raise ValueError(f"{int_vep_version=} - we require min of v110")

    extension = get_extension_without_gzip(data_path)
    fields = [cvf.source_field for cvf in cvf_list if cvf.source_field]  # Strip empty/falsey
    params = {
        "file": data_path,
        "type": "overlap",
        "num_records": "all",  # Display all (defaults to 50 then "...")
    }
    if prefix == "RepeatMasker":
        params["num_records"] = "10000"  # repeat masker can get ridiculous - truncates with "..."

    if extension == 'vcf':
        field_delimiter = "%"

        params["short_name"] = prefix
        params["fields"] = field_delimiter.join(fields)
        params["format"] = 'vcf'

        if prefix in ["gnomAD_SV", "gnomAD_SV_name"]:
            del params["fields"]  # Leave fields blank so that it uses ID
            params["type"] = "overlap"
            overlap_cutoff = str(int(100 * settings.ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION))
            params["overlap_cutoff"] = overlap_cutoff
            params["same_type"] = "1"

            # gnomad_sv_name doesn't have any fields
            if prefix == "gnomAD_SV":
                params["fields"] = field_delimiter.join(fields)
                params["coords"] = "1"
        else:
            params["type"] = "exact"

    else:
        if len(cvf_list) != 1:
            raise ValueError(f"Expected exactly 1 ColumnVEPField source field for VEP custom: {prefix}, {cvf_list=}")

        cvf = cvf_list[0]
        if cvf.summary_stats:
            params["short_name"] = prefix
            params["summary_stats"] = cvf.summary_stats
            params["num_records"] = 0
        else:
            params["short_name"] = cvf.source_field

        if extension == 'bed':
            fmt = "bed"
        elif extension == 'bw':
            fmt = "bigwig"
        else:
            msg = "Don't know how to handle custom data: {data_path}"
            raise ValueError(msg)

        params["format"] = fmt

    command = ",".join([f"{k}={v}" for k, v in params.items()])
    logging.info(command)
    return ["--custom", command]


def get_vep_command(vcf_filename, output_filename, genome_build: GenomeBuild, annotation_consortium,
                    pipeline_type: VariantAnnotationPipelineType):
    vc = VEPConfig(genome_build)
    vep_cmd = os.path.join(settings.ANNOTATION_VEP_CODE_DIR, "vep")

    #
    cmd = [
        vep_cmd,
        "-i", vcf_filename,
        "-o", output_filename,
        "--cache",
        "--cache_version", str(vc.cache_version),
        "--dir_cache", settings.ANNOTATION_VEP_CACHE_DIR,
        "--dir_plugins", settings.ANNOTATION_VEP_PLUGINS_DIR,
        # Need to provide VEP a fasta rather than use the default - https://github.com/Ensembl/VEP_plugins/issues/708
        "--fasta", vc["fasta"],
        "--assembly", genome_build.name,
        "--offline", "--use_given_ref", "--vcf", "--compress_output", "gzip",
        "--force_overwrite", "--flag_pick", "--exclude_predicted", "--no_stats",
        "--check_existing",  # COSMIC ids
        "--no_escape",  # Don't URI escape HGVS strings

        # flags for fields
        "--hgvs",
        "--symbol",
        "--numbers",  # EXON/INTRON numbers
        "--domains",
        # "--regulatory", # don't know how to deal with ENSM00525026610
        "--canonical",
        "--protein",
        "--biotype",
        "--transcript_version",  # Makes Ensembl Transcript IDs in Feature have version (RefSeq ones already do)
        "--af",
        "--pubmed",
        "--variant_class",
    ]

    if vc.use_sift:
        cmd.extend(["--sift", "b"])

    if settings.ANNOTATION_VEP_PICK_ORDER:
        # @see https://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_pick_order
        cmd.extend(["--pick_order", settings.ANNOTATION_VEP_PICK_ORDER])

    if settings.ANNOTATION_VEP_DISTANCE is not None:
        cmd.extend(["--distance", str(settings.ANNOTATION_VEP_DISTANCE)])

    if max_sv_size := settings.ANNOTATION_VEP_SV_MAX_SIZE:
        vep_default_max_sv_size = 10_000_000
        if max_sv_size != vep_default_max_sv_size:
            cmd.extend(["--max_sv_size", str(max_sv_size)])

    if annotation_consortium == AnnotationConsortium.REFSEQ:
        cmd.append("--refseq")

    if settings.ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT:
        cmd.insert(0, settings.ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT)

    if settings.ANNOTATION_VEP_FORK and settings.ANNOTATION_VEP_FORK > 1:
        cmd.extend(["--fork", str(settings.ANNOTATION_VEP_FORK)])

    if buffer_size := settings.ANNOTATION_VEP_BUFFER_SIZE.get(pipeline_type):
        cmd.extend(["--buffer_size", str(buffer_size)])

    if settings.ANNOTATION_VEP_ARGS:
        cmd.extend(settings.ANNOTATION_VEP_ARGS)

    # Restrict the number of BRCA1 transcripts to stop a blowout of VEP RAM usage
    # RefSeq GRCh38 annotation (GCF_000001405.40-RS_2023_10) has 368 transcripts for BRCA1, up from 6 previously
    # see https://github.com/SACGF/variantgrid/issues/1228
    try:
        transcript_blocklist_filename = vc["transcript_blocklist"]
        cmd.extend(["--transcript_filter", f"not stable_id in {transcript_blocklist_filename}"])
    except KeyError:
        pass

    if pipeline_type == VariantAnnotationPipelineType.STANDARD:
        # TODO: At the moment we just skip everything, perhaps we should:
        # a) Add ColumnVEPField.pipeline_type (or maybe have a list of all types) so we can configure it to run on both
        # b) Run another pipeline for CNVs

        cmd.extend([
            # Plugins that don't require data
            "--plugin", "Grantham",
            "--plugin", "SpliceRegion",
        ])

        # Plugins that require data - ok for these to fail when retrieving vep config
        plugin_data_func = {
            VEPPlugin.MASTERMIND: lambda: f"Mastermind,{vc['mastermind']},1",  # 1 to not filter
            VEPPlugin.MAXENTSCAN: lambda: f"MaxEntScan,{vc['maxentscan']}",
            VEPPlugin.DBNSFP: lambda: _get_dbnsfp_plugin_command(genome_build, vc),
            VEPPlugin.DBSCSNV: lambda: f"dbscSNV,{vc['dbscsnv']}",
            VEPPlugin.SPLICEAI: lambda: f"SpliceAI,snv={vc['spliceai_snv']},indel={vc['spliceai_indel']}"
        }

        if vc.columns_version >= 2:
            cmd.extend(["--plugin", "NMD"])

        if vc.columns_version >= 3:
            plugin_data_func.update({
                VEPPlugin.MAVEDB: lambda: f"MaveDB,file={vc['mave']},single_aminoacid_changes=0,transcript_match=0 ",
            })

    else:
        plugin_data_func = {}  # No plugins for SVs

    # Custom - filter_for(vep_config=vc) drops anything whose data file isn't configured,
    # so we don't need to probe vc[prefix_lc] separately.
    for vep_custom in VEPCustom:
        prefix = vep_custom.label
        # Match original ColumnVEPField.get(): distinct source_field, ordered by source_field
        # (postgres default collation = case-insensitive).
        cvf_list = sorted(
            {c.source_field: c for c in vep_columns.filter_for(
                vep_config=vc,
                pipeline_type=pipeline_type,
                vep_custom=vep_custom,
            )}.values(),
            key=lambda c: (c.source_field or "").lower(),
        )
        if cvf_list:
            cfg = vc[prefix.lower()]
            cmd.extend(_get_custom_params_list(cvf_list, prefix, cfg))

    for vep_plugin, plugin_arg_func in plugin_data_func.items():
        try:
            cmd.extend(["--plugin", plugin_arg_func()])
        except Exception as e:
            logging.warning(e)
            logging.warning("No annotation set for plugin: %s", vep_plugin)

    return cmd


def run_vep(vcf_filename, output_filename, genome_build: GenomeBuild, annotation_consortium,
            pipeline_type: VariantAnnotationPipelineType):
    """ executes VEP command. Returns (command_line, code, stdout, stderr) """

    cmd = get_vep_command(vcf_filename, output_filename, genome_build, annotation_consortium, pipeline_type)
    return execute_cmd(cmd, shell=False)


def get_vep_version(genome_build: GenomeBuild, annotation_consortium):
    """ returns dictionary of VEP and database versions """

    vcf_filename = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, "fake.vcf")
    if not os.path.exists(vcf_filename):
        mk_path_for_file(vcf_filename)
        FAKE_VCF_STR = "##fileformat=VCFv4.1\n" \
                       "#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO\n" \
                       "1    123456    .    G    T    .    .        \n"

        with open(vcf_filename, "w", encoding="utf-8") as f:
            f.write(FAKE_VCF_STR)

    # Use unique name so we don't get collisions
    output_basename = f"fake.vep_annotated_{genome_build.name}.{uuid.uuid4()}.vcf.gz"
    output_filename = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, output_basename)
    returncode, std_out, std_err = run_vep(vcf_filename, output_filename, genome_build,
                                           annotation_consortium, VariantAnnotationPipelineType.STANDARD)
    if returncode != 0:
        logging.info(std_out)
        logging.error(std_err)
        raise ValueError(f"VEP returned {returncode}")
    vep_version = get_vep_version_from_vcf(output_filename)
    os.remove(output_filename)
    return vep_version


def _spliceai_label(spliceai_filename: str) -> str:
    flavour = "masked" if "masked" in os.path.basename(spliceai_filename).lower() else "raw"
    version = _spliceai_version_from_vcf_header(spliceai_filename)
    if version:
        return f"{flavour} {version}"
    return flavour


def _spliceai_version_from_vcf_header(spliceai_filename: str) -> Optional[str]:
    try:
        with open_handle_gzip(spliceai_filename, "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    break
                if "ID=SpliceAI" in line:
                    if m := re.search(r"SpliceAIv(\d+(?:\.\d+)*)", line):
                        return m.group(1)
    except OSError:
        logging.warning("Could not read SpliceAI VCF header: %s", spliceai_filename)
    return None


def vep_dict_to_variant_annotation_version_kwargs(vep_config, vep_version_dict: dict) -> dict:
    def _vep_int_version(vep_string_version):
        m = re.match(r"v(\d+)", vep_string_version)
        return int(m.group(1))

    # We strip off the hash ie '110.73b02d8' -> '110' so that we can re-use annotation
    # if VEP do a minor bugfix change to their releases etc
    def _major_version(version_str) -> str:
        return version_str.split(".")[0]

    FIELD_CONVERSION = {
        "vep": _vep_int_version,
        "cosmic": int,
        "dbsnp": int,
        "ensembl": _major_version,
        "ensembl_io": _major_version,
        "ensembl_variation": _major_version,
    }

    FIELD_LOOKUP = {
        "VEP": "vep",
        'ensembl': 'ensembl',
        'ensembl-funcgen': 'ensembl_funcgen',
        'ensembl-variation': 'ensembl_variation',
        'ensembl-io': 'ensembl_io',
        '1000genomes': 'thousand_genomes',
        'COSMIC': 'cosmic',
        'HGMD-PUBLIC': 'hgmd',
        'assembly': 'assembly',
        'dbSNP': 'dbsnp',
        'gencode': 'gencode',
        'genebuild': 'genebuild',
        'gnomAD': 'gnomad',
        'refseq': 'refseq',
        'regbuild': 'regbuild',
        'sift': 'sift',
    }

    kwargs = {}
    for json_field, python_field in FIELD_LOOKUP.items():
        if value := vep_version_dict.get(json_field):
            converter = FIELD_CONVERSION.get(python_field)
            if converter:
                value = converter(value)
            kwargs[python_field] = value

    genome_build = vep_config.genome_build
    kwargs["genome_build"] = genome_build
    kwargs["vep_cache"] = vep_config.cache_version
    kwargs["annotation_consortium"] = vep_config.annotation_consortium
    kwargs["columns_version"] = vep_config.columns_version
    distance = getattr(settings, "ANNOTATION_VEP_DISTANCE", None)
    if distance is None:
        distance = 5000
    kwargs["distance"] = distance

    # Plugins are optional
    try:
        dbnsfp_path = vep_config["dbnsfp"]  # KeyError if not set in settings
        if m := re.match(r".*/dbNSFP_?(.*?)\.(GRCh37|GRCh38|hg19|hg38)", dbnsfp_path, flags=re.IGNORECASE):
            kwargs["dbnsfp"] = m.group(1)
        else:
            msg = f"Couldn't determine dbNSFP version from file: {dbnsfp_path}"
            raise ValueError(msg)
    except KeyError:
        pass

    # we use our own gnomAD custom annotation, not the default VEP one
    candidates = [c for c in vep_columns.for_variant_grid_column('gnomad_af')
                  if genome_build.name in c.genome_builds
                  and c.applies_to(columns_version=vep_config.columns_version)]
    if candidates:
        cvf = candidates[0]
        try:
            # annotation_data/GRCh37/gnomad2.1.1_GRCh37_combined_af.vcf.bgz
            # gnomad3.1_GRCh38_merged.vcf.bgz
            gnomad_filename = vep_config[cvf.vep_custom.label.lower()]
            if os.path.exists(gnomad_filename):
                gnomad_version = parse_gnomad_version_from_filename(gnomad_filename)
                if gnomad_version is not None:
                    kwargs["gnomad"] = gnomad_version
                else:
                    msg = f"Couldn't determine gnomAD version from file: {os.path.basename(gnomad_filename)}"
                    raise ValueError(msg)
        except KeyError:
            pass  # Will just use VEP values

    try:
        # COSMIC isn't in T2T
        cosmic_filename = vep_config["cosmic"]
        if os.path.exists(cosmic_filename):
            cosmic_basename = os.path.basename(cosmic_filename)
            if m := re.match(r"^Cosmic.*_v(\d{2,})_.*.vcf.gz", cosmic_basename):
                kwargs["cosmic"] = int(m.group(1))
    except KeyError:
        pass

    try:
        spliceai_snv_filename = vep_config["spliceai_snv"]
        if spliceai_snv_filename and os.path.exists(spliceai_snv_filename):
            kwargs["spliceai"] = _spliceai_label(spliceai_snv_filename)
    except KeyError:
        pass

    try:
        # denovo-db is GRCh37/38 only
        denovo_db_filename = vep_config["denovo_db"]
        if denovo_db_filename and os.path.exists(denovo_db_filename):
            denovo_db_basename = os.path.basename(denovo_db_filename)
            # e.g. denovo-db.variants.v.1.6.1.GRCh37.vcf.gz
            if m := re.match(r"^denovo-db\.variants\.v\.(?P<version>[\d.]+)\.(GRCh37|GRCh38)\.vcf\.gz$",
                             denovo_db_basename):
                kwargs["denovo_db"] = m.group("version")
            else:
                msg = f"Couldn't determine denovo-db version from file: {denovo_db_basename}"
                raise ValueError(msg)
    except KeyError:
        pass

    return kwargs


def get_vep_variant_annotation_version_kwargs(genome_build: GenomeBuild):
    vep_config = VEPConfig(genome_build)
    if settings.ANNOTATION_VEP_FAKE_VERSION:
        return get_fake_vep_version(genome_build, vep_config.annotation_consortium, vep_config.columns_version)

    vep_version_dict = get_vep_version(genome_build, vep_config.annotation_consortium)
    kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_config, vep_version_dict)
    return kwargs


def get_vep_version_from_vcf(output_filename):
    VEP_VERSIONS_LINE_START = "##VEP="

    num_lines_read = 0
    with open_handle_gzip(output_filename, "rt") as f:
        for line in f:
            num_lines_read += 1
            if line.startswith("#"):
                if line.startswith(VEP_VERSIONS_LINE_START):
                    return vep_parse_version_line(line)
            else:
                break

    raise ValueError(f"{output_filename}: Could not find line in header starting with '{VEP_VERSIONS_LINE_START}' "
                     + f"(read {num_lines_read} lines).")

def vep_parse_version_line(line):
    try:
        line = line[2:].strip()  # Remove hashes
        value_sep = "="

        lexer = shlex(line, posix=True)
        # need "-" for e.g. 'ensembl-funcgen' and '.' as they have unquoted versions like: 97.378db18
        lexer.wordchars += value_sep + "-."
        lexer.whitespace = " "
        vep_version_dict = {"refseq": ""}  # Defaults - refseq won't be populated with Ensembl annotation
        vep_version_dict.update(dict([word.split(value_sep, maxsplit=1) for word in lexer]))
        return vep_version_dict
    except:
        logging.error("vep_parse_version_line - couldn't parse '%s'", line)
        raise


def _vep_check_version_match(variant_annotation_version, vep_version_kwargs: dict, vep_version_desc: str):
    for k, v in vep_version_kwargs.items():
        version_value = getattr(variant_annotation_version, k)
        if version_value != v:
            msg = f"Annotation/VEP out of sync! Version '{variant_annotation_version}' and {vep_version_desc} differ, " \
                  f"KEY: '{k}' annotation: '{version_value}', vcf: '{v}'"
            raise VEPVersionMismatchError(msg)


def vep_check_annotated_file_version_match(variant_annotation_version, filename):
    """ Load VEP VCF, check VEP= line and make sure that values match expected VariantAnnotationVersion """
    vep_config = VEPConfig(variant_annotation_version.genome_build)
    vep_dict = get_vep_version_from_vcf(filename)
    vep_version_kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_config, vep_dict)
    _vep_check_version_match(variant_annotation_version, vep_version_kwargs,
                            f"VEP annotated VCF: '{filename}'")


def vep_check_command_line_version_match(variant_annotation_version):
    """ Check vs what will get executed on command line """
    vep_config = VEPConfig(variant_annotation_version.genome_build)
    vep_version = get_vep_version(variant_annotation_version.genome_build, vep_config.annotation_consortium)
    vep_version_kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_config, vep_version)
    _vep_check_version_match(variant_annotation_version, vep_version_kwargs, "VEP command line")
