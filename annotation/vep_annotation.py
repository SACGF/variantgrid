import logging
import os
import re
from shlex import shlex
from typing import Dict

from django.conf import settings

from annotation.fake_annotation import get_fake_vep_version
from annotation.models.models import ColumnVEPField
from annotation.models.models_enums import VEPPlugin, VEPCustom
from genes.models_enums import AnnotationConsortium
from library.file_utils import get_extension_without_gzip, mk_path_for_file, open_handle_gzip
from library.utils import get_single_element, execute_cmd
from snpdb.models.models_genome import GenomeBuild


class VEPVersionMismatchError(ValueError):
    pass


class VEPConfig:

    def __init__(self, genome_build: GenomeBuild):
        self.annotation_consortium = genome_build.annotation_consortium
        self.genome_build = genome_build
        self.vep_data = genome_build.settings["vep_config"]

    def __getitem__(self, key):
        """ Throws KeyError if missing """
        value = self.vep_data[key]  # All callers need to catch KeyError
        return os.path.join(settings.ANNOTATION_VEP_BASE_DIR, value)


def _get_dbnsfp_plugin_command(genome_build: GenomeBuild, vc: VEPConfig):
    """ Build from ColumnVEPField.source_field where vep_plugin = DBNSFP """

    dbnsfp_data_path = vc["dbnsfp"]
    fields = ColumnVEPField.get_source_fields(genome_build, vep_plugin=VEPPlugin.DBNSFP)
    joined_columns = ",".join(fields)
    return f"dbNSFP,{dbnsfp_data_path},{joined_columns}"


def _get_custom_params_list(fields, prefix, data_path) -> list:
    extension = get_extension_without_gzip(data_path)

    if extension == 'vcf':
        joined_columns = ",".join(fields)
        command = f"{data_path},{prefix},vcf,exact,0,{joined_columns}"
    else:
        try:
            field = get_single_element(fields)
        except Exception as e:
            msg = f"Expected exactly 1 ColumnVEPField source field for VEP custom: {prefix}"
            raise ValueError(msg) from e

        if extension == 'bed':
            command = f"{data_path},{field},bed,overlap"
        elif extension == 'bw':
            command = f"{data_path},{field},bigwig,overlap"
        else:
            msg = "Don't know how to handle custom data: {data_path}"
            raise ValueError(msg)
    return ["--custom", command]


def get_vep_command(vcf_filename, output_filename, genome_build: GenomeBuild, annotation_consortium):
    vc = VEPConfig(genome_build)
    vep_cmd = os.path.join(settings.ANNOTATION_VEP_CODE_DIR, "vep")
    reference_fasta = genome_build.reference_fasta

    cmd = [
        vep_cmd,
        "-i", vcf_filename,
        "-o", output_filename,
        "--cache", "--dir", settings.ANNOTATION_VEP_CACHE_DIR,
        "--fasta", reference_fasta,
        "--assembly", genome_build.name,
        "--offline", "--use_given_ref", "--vcf", "--compress_output", "gzip",
        "--force_overwrite", "--flag_pick", "--exclude_predicted", "--no_stats",
        "--check_existing",  # COSMIC ids
        "--no_escape",  # Don't URI escape HGVS strings

        # flags for fields
        "--sift", "b",
        "--uniprot",
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
        # Plugins that don't require data
        "--plugin", "Grantham",
        "--plugin", "SpliceRegion",
        "--plugin", "LoFtool",
    ]

    if settings.ANNOTATION_VEP_PICK_ORDER:
        # @see https://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_pick_order
        cmd.extend(["--pick_order", settings.ANNOTATION_VEP_PICK_ORDER])

    if settings.ANNOTATION_VEP_DISTANCE is not None:
        cmd.extend(["--distance", str(settings.ANNOTATION_VEP_DISTANCE)])

    # Plugins that require data
    PLUGINS = {VEPPlugin.MASTERMIND: lambda: f"Mastermind,{vc['mastermind']},1",  # 1 to not filter
               VEPPlugin.MAXENTSCAN: lambda: f"MaxEntScan,{vc['maxentscan']}",
               VEPPlugin.DBNSFP: lambda: _get_dbnsfp_plugin_command(genome_build, vc),
               VEPPlugin.DBSCSNV: lambda: f"dbscSNV,{vc['dbscsnv']}",
               VEPPlugin.SPLICEAI: lambda: f"SpliceAI,snv={vc['spliceai_snv']},indel={vc['spliceai_indel']}"}

    for vep_plugin, plugin_arg_func in PLUGINS.items():
        try:
            cmd.extend(["--plugin", plugin_arg_func()])
        except Exception as e:
            logging.warning(e)
            logging.warning("No annotation set for plugin: %s", vep_plugin)

    # Custom
    for vep_custom, prefix in dict(VEPCustom.choices).items():
        try:
            if fields := ColumnVEPField.get_source_fields(genome_build, vep_custom=vep_custom):
                prefix_lc = prefix.lower()
                if cfg := vc[prefix_lc]:  # annotation settings are lower case
                    cmd.extend(_get_custom_params_list(fields, prefix, cfg))
                else:
                    logging.info("Skipping due to settings.ANNOTATION[%s][vep_config][%s] = None",
                                 genome_build.name, prefix_lc)
        except Exception as e:
            logging.warning(e)
            # Not all annotations available for all builds - ok to just warn
            logging.warning("Skipped custom annotation: %s", prefix)

    if annotation_consortium == AnnotationConsortium.REFSEQ:
        cmd.append("--refseq")

    if settings.ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT:
        cmd.insert(0, settings.ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT)

    if settings.ANNOTATION_VEP_FORK and settings.ANNOTATION_VEP_FORK > 1:
        cmd.extend(["--fork", str(settings.ANNOTATION_VEP_FORK)])

    if settings.ANNOTATION_VEP_ARGS:
        cmd.extend(settings.ANNOTATION_VEP_ARGS)

    return cmd


def run_vep(vcf_filename, output_filename, genome_build: GenomeBuild, annotation_consortium):
    """ executes VEP command. Returns (command_line, code, stdout, stderr) """

    cmd = get_vep_command(vcf_filename, output_filename, genome_build, annotation_consortium)
    logging.info("Executing VEP:")
    logging.info(" ".join(cmd))
    return execute_cmd(cmd)


def get_vep_version(genome_build: GenomeBuild, annotation_consortium):
    """ returns dictionary of VEP and database versions """

    vcf_filename = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, "fake.vcf")
    if not os.path.exists(vcf_filename):
        mk_path_for_file(vcf_filename)
        FAKE_VCF_STR = "##fileformat=VCFv4.1\n" \
                       "#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO\n" \
                       "1    123456    .    G    T    .    .        \n"

        with open(vcf_filename, "w") as f:
            f.write(FAKE_VCF_STR)

    output_basename = f"fake.vep_annotated_{genome_build.name}.vcf.gz"
    output_filename = os.path.join(settings.ANNOTATION_VCF_DUMP_DIR, output_basename)
    returncode, std_out, std_err = run_vep(vcf_filename, output_filename, genome_build, annotation_consortium)
    if returncode != 0:
        logging.info(std_out)
        logging.error(std_err)
        raise ValueError(f"VEP returned {returncode}")
    return get_vep_version_from_vcf(output_filename)


def vep_dict_to_variant_annotation_version_kwargs(vep_version_dict: Dict) -> Dict:
    def vep_int_version(vep_string_version):
        m = re.match(r"v(\d+)", vep_string_version)
        return int(m.group(1))

    FIELD_CONVERSION = {
        "vep": vep_int_version,
        "cosmic": int,
        "dbsnp": int,
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
        value = vep_version_dict.get(json_field)
        converter = FIELD_CONVERSION.get(python_field)
        if converter:
            value = converter(value)
        kwargs[python_field] = value
    return kwargs


def get_vep_variant_annotation_version_kwargs(genome_build: GenomeBuild):
    vep_config = VEPConfig(genome_build)
    if settings.ANNOTATION_VEP_FAKE_VERSION:
        return get_fake_vep_version(genome_build, vep_config.annotation_consortium)

    vep_version_dict = get_vep_version(genome_build, vep_config.annotation_consortium)
    kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_version_dict)

    vc = VEPConfig(genome_build)
    kwargs["genome_build"] = genome_build
    kwargs["annotation_consortium"] = vc.annotation_consortium
    distance = getattr(settings, "ANNOTATION_VEP_DISTANCE", None)
    if distance is None:
        distance = 5000
    kwargs["distance"] = distance

    # Plugins are optional
    try:
        dbnsfp_path = vc["dbnsfp"]  # KeyError if not set in settings
        if m := re.match(r".*/dbNSFP_?(.*?)\.(GRCh37|GRCh38|hg19|hg38)", dbnsfp_path, flags=re.IGNORECASE):
            kwargs["dbnsfp"] = m.group(1)
        else:
            msg = f"Couldn't determine dbNSFP version from file: {dbnsfp_path}"
            raise ValueError(msg)
    except KeyError:
        kwargs["dbnsfp"] = 'n/a'

    return kwargs


def get_vep_version_from_vcf(output_filename):
    VEP_VERSIONS_LINE_START = "##VEP"

    with open_handle_gzip(output_filename, "rt") as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith(VEP_VERSIONS_LINE_START):
                    return vep_parse_version_line(line)
            else:
                break

    msg = f"Could not find line in header starting with '{VEP_VERSIONS_LINE_START}'"
    raise ValueError(msg)


def vep_parse_version_line(line):
    try:
        line = line[2:].strip()  # Remove hashes
        value_sep = "="

        lexer = shlex(line, posix=True)
        # need "-" for eg 'ensembl-funcgen' and '.' as they have unquoted versions like: 97.378db18
        lexer.wordchars += value_sep + "-."
        lexer.whitespace = " "
        vep_version_dict = {"refseq": ""}  # Defaults - refseq won't be populated with Ensembl annotation
        vep_version_dict.update(dict([word.split(value_sep, maxsplit=1) for word in lexer]))
        return vep_version_dict
    except:
        logging.error("vep_parse_version_line - couldn't parse '%s'", line)
        raise


def vep_check_version_match(variant_annotation_version, filename: str):
    """ Load VEP VCF, check VEP= line and make sure that values match expected VariantAnnotationVersion """
    vep_dict = get_vep_version_from_vcf(filename)
    kwargs = vep_dict_to_variant_annotation_version_kwargs(vep_dict)
    for k, v in kwargs.items():
        version_value = getattr(variant_annotation_version, k)
        if version_value != v:
            msg = f"Annotation/VCF out of sync! Version '{variant_annotation_version}' and VCF '{filename}' differ, " \
                  f"KEY: '{k}' annotation: '{version_value}', vcf: '{v}'"
            raise VEPVersionMismatchError(msg)
