import logging
import os
import shutil
from pathlib import Path
from subprocess import CalledProcessError

import celery
import pandas as pd
from django.conf import settings
from django.db import IntegrityError

from library.django_utils.django_file_utils import get_import_processing_dir
from library.log_utils import log_traceback, get_traceback
from library.utils import execute_cmd
from patients.models_enums import Zygosity
from snpdb.models import VCF, Cohort, Trio, SomalierVCFExtract, SomalierConfig, SomalierAncestryRun, \
    SomalierCohortRelate, SomalierRelate, SomalierTrioRelate, SomalierAncestry, SuperPopulationCode, \
    SomalierSampleExtract, ProcessingStatus, SomalierAllSamplesRelate, AbstractSomalierModel, SomalierRelatePairs
from snpdb.variants_to_vcf import vcf_export_to_file


@celery.task
def somalier_vcf_id(vcf_id: int):
    """ Extract, Ancestry, and Relate (Cohort) """
    vcf = VCF.objects.get(pk=vcf_id)
    SomalierVCFExtract.objects.filter(vcf=vcf).delete()  # Delete any previous
    vcf_extract = SomalierVCFExtract.objects.create(vcf=vcf)
    if not vcf.has_genotype:
        vcf_extract.status = ProcessingStatus.SKIPPED
        vcf_extract.save()
        return

    processing_dir = get_import_processing_dir(vcf_extract.pk, "somalier_vcf_extract")

    try:
        cfg = SomalierConfig()
        vcf_filename = _write_somalier_vcf(cfg, processing_dir, vcf_extract)

        # Extract
        somalier_bin = cfg.get_annotation("command")
        extract_cmd = [somalier_bin, "extract",
                       "--out-dir", vcf_extract.get_somalier_dir(),
                       "--sites", cfg.get_sites(vcf.genome_build),
                       "--fasta", vcf.genome_build.reference_fasta,
                       vcf_filename]
        vcf_extract.execute(extract_cmd)

        _somalier_ancestry(vcf_extract)

        SomalierCohortRelate.objects.filter(cohort=vcf.cohort).delete()  # Delete any previous versions
        relate = SomalierCohortRelate.objects.create(cohort=vcf.cohort, cohort_version=vcf.cohort.version)
        _somalier_relate(relate)
    except CalledProcessError:
        log_traceback()
        logging.info("Somalier failed - carrying on with rest of pipeline")

    if settings.VCF_IMPORT_DELETE_TEMP_FILES_ON_SUCCESS:
        shutil.rmtree(processing_dir)


def _write_somalier_vcf(cfg: SomalierConfig, processing_dir, vcf_extract: SomalierVCFExtract):
    vcf = vcf_extract.vcf
    sites_vcf = cfg.get_sites_vcf(vcf.genome_build)
    sites_qs = sites_vcf.get_variant_qs()
    exported_vcf_filename = os.path.join(processing_dir, f"vcf_{vcf.pk}.vcf.bgz")
    sample_zygosity_count = vcf_export_to_file(vcf, exported_vcf_filename, original_qs=sites_qs,
                                               sample_name_func=AbstractSomalierModel.sample_name)
    ZYG_LOOKUP = {"ref_count": Zygosity.HOM_REF,
                  "het_count": Zygosity.HET,
                  "hom_count": Zygosity.HOM_ALT,
                  "unk_count": Zygosity.UNKNOWN_ZYGOSITY}
    for sample, zy in sample_zygosity_count.items():
        zyg_kwargs = {k: zy[v] for k, v in ZYG_LOOKUP.items()}
        SomalierSampleExtract.objects.create(vcf_extract=vcf_extract, sample=sample, **zyg_kwargs)

    tabix_command = ["tabix", exported_vcf_filename]
    return_code, stdout, stderr = execute_cmd(tabix_command)
    print(f"return_code: {return_code}, stdout: {stdout}, stderr: {stderr}")
    return exported_vcf_filename


def _somalier_ancestry(vcf_extract: SomalierVCFExtract):
    cfg = SomalierConfig()
    ancestry_run = SomalierAncestryRun.objects.create(vcf_extract=vcf_extract)
    compare_samples = os.path.join(cfg.get_annotation("ancestry_somalier_dir"), "*.somalier")
    ancestry_report_dir = Path(ancestry_run.get_report_dir())
    ancestry_report_dir.mkdir(parents=True, exist_ok=True)
    somalier_bin = cfg.get_annotation("command")
    ancestry_cmd = [somalier_bin, "ancestry", "--labels", cfg.get_annotation("ancestry_labels"),
                    compare_samples, "++"] + ancestry_run.get_sample_somalier_filenames()
    # Force Somalier to only use 1 thread - have had it run very slow with multi-core
    # https://github.com/brentp/somalier/issues/61#issuecomment-750492570
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    ancestry_run.execute(ancestry_cmd, cwd=ancestry_report_dir, env=env)

    # Use TSV to write sample specific files
    df = pd.read_csv(ancestry_report_dir / "somalier-ancestry.somalier-ancestry.tsv", sep='\t', index_col=0)
    for sample_extract in vcf_extract.somaliersampleextract_set.all():
        row = df.loc[AbstractSomalierModel.sample_name(sample_extract.sample)]
        predicted_ancestry = getattr(SuperPopulationCode, row["predicted_ancestry"])  # Get enum ie 'EAS'
        SomalierAncestry.objects.create(ancestry_run=ancestry_run,
                                        sample_extract=sample_extract,
                                        predicted_ancestry=predicted_ancestry,
                                        EAS_prob=row["EAS_prob"],
                                        AFR_prob=row["AFR_prob"],
                                        AMR_prob=row["AMR_prob"],
                                        SAS_prob=row["SAS_prob"],
                                        EUR_prob=row["EUR_prob"])


def _somalier_relate(somalier_relate: SomalierRelate) -> Path:
    """ Returns path of relate output """
    cfg = SomalierConfig()
    somalier_bin = cfg.get_annotation("command")
    processing_dir = get_import_processing_dir(somalier_relate.pk, "somalier_relate")

    command = [somalier_bin, "relate"]
    if not somalier_relate.is_joint_called_vcf():
        # Somalier --unknown    set unknown genotypes to hom-ref.
        # it is often preferable to use this with VCF samples that were not jointly called
        command += ["--unknown"]

    if somalier_relate.has_ped_file():
        ped_filename = os.path.join(processing_dir, "temp.ped")
        somalier_relate.write_ped_file(ped_filename)
        command += ["--ped", ped_filename]

    command += somalier_relate.get_sample_somalier_filenames()

    somalier_relate_dir = Path(cfg.related_dir(somalier_relate.uuid))
    somalier_relate_dir.mkdir(parents=True, exist_ok=True)
    somalier_relate.execute(command, cwd=somalier_relate_dir)

    shutil.rmtree(processing_dir)

    return somalier_relate_dir


@celery.task
def somalier_cohort_relate(cohort_id: int):
    cohort = Cohort.objects.get(pk=cohort_id)
    relate = SomalierCohortRelate.objects.create(cohort=cohort, cohort_version=cohort.version)
    _somalier_relate(relate)


@celery.task
def somalier_trio_relate(trio_id: int):
    trio = Trio.objects.get(pk=trio_id)
    relate = SomalierTrioRelate.objects.create(trio=trio)
    _somalier_relate(relate)


@celery.task
def somalier_all_samples():
    somalier_settings = settings.SOMALIER["relatedness"]
    all_samples = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.PROCESSING)
    try:
        related_dir = _somalier_relate(all_samples)
        pairs_filename = os.path.join(related_dir, "somalier.pairs.tsv")
        df = pd.read_csv(pairs_filename, sep='\t')
        shared_het_mask = df["shared_hets"] >= somalier_settings["min_shared_hets"]
        shared_hom_mask = df["shared_hom_alts"] > somalier_settings["min_shared_hom_alts"]
        relateness_mask = df["relatedness"] > somalier_settings["min_relatedness"]
        df = df[shared_het_mask & shared_hom_mask & relateness_mask]

        for _, row in df.iterrows():
            row_data = dict(row)
            row_data["relate"] = all_samples
            row_data.pop("expected_relatedness")
            # Sample_ids are at the start
            sample_a_id = AbstractSomalierModel.sample_name_to_id(row_data.pop("#sample_a"))
            sample_b_id = AbstractSomalierModel.sample_name_to_id(row_data.pop("sample_b"))

            try:
                SomalierRelatePairs.objects.update_or_create(sample_a_id=sample_a_id, sample_b_id=sample_b_id,
                                                             defaults=row_data)
            except IntegrityError:
                pass  # Sample was deleted - just don't store

        all_samples.status = ProcessingStatus.SUCCESS
    except Exception as e:
        tb = get_traceback()
        logging.error(tb)
        all_samples.error_exception = tb
        all_samples.status = ProcessingStatus.ERROR

    all_samples.save()
