"""
Somalier stages: extract a VCF's calls at the somalier sites, predict ancestry, and relate samples
(within a cohort, a trio, or nightly across every sample on the box).

Entry points: `somalier_vcf_id` (run by the VCF import pipeline), `somalier_cohort_relate`,
`somalier_trio_relate` and the beat-scheduled `somalier_all_samples`. Every stage records its own
ProcessingStatus, so one failing never takes the others - or the import - down with it.
"""
import glob
import logging
import os
import shutil
import time
from contextlib import contextmanager
from datetime import timedelta
from pathlib import Path

import celery
import pandas as pd
from django.conf import settings
from django.db.models import F, Max
from django.utils import timezone

from library.django_utils.django_file_utils import get_import_processing_dir
from library.log_utils import get_traceback
from library.utils import execute_cmd
from patients.models_enums import Zygosity
from snpdb.models import (
    VCF,
    AbstractSomalierModel,
    Cohort,
    ProcessingStatus,
    Sample,
    SomalierAllSamplesRelate,
    SomalierAncestry,
    SomalierAncestryRun,
    SomalierCohortRelate,
    SomalierConfig,
    SomalierRelate,
    SomalierRelatePairs,
    SomalierSampleExtract,
    SomalierTrioRelate,
    SomalierVCFExtract,
    SuperPopulationCode,
    Trio,
)
from snpdb.variants_to_vcf import vcf_export_to_file

# A PROCESSING extract older than this was left behind by a worker that died, not by a live run
STALE_PROCESSING_AGE = timedelta(hours=1)


@contextmanager
def _stage_timing(stage: str, vcf_id: int):
    """ #1147 - so how long each stage takes can be answered from the celery log on prod """
    start = time.time()
    try:
        yield
    finally:
        logging.info("somalier %s vcf=%s took %.1fs", stage, vcf_id, time.time() - start)


def _record_error(somalier_model: AbstractSomalierModel, stage: str):
    """ execute() already wrote the command's return code/stderr - keep that over a traceback """
    tb = get_traceback()
    logging.error("Somalier %s failed: %s", stage, tb)
    if somalier_model.status != ProcessingStatus.ERROR:
        somalier_model.error_exception = tb
        somalier_model.status = ProcessingStatus.ERROR
        somalier_model.save()


def _set_skipped(somalier_model: AbstractSomalierModel, reason: str):
    somalier_model.status = ProcessingStatus.SKIPPED
    somalier_model.error_exception = reason
    somalier_model.save()


def _create_vcf_extract(vcf: VCF) -> SomalierVCFExtract:
    """ None when another run for this VCF is still going """
    if existing := SomalierVCFExtract.objects.filter(vcf=vcf).first():
        if existing.status == ProcessingStatus.PROCESSING and \
                existing.modified > timezone.now() - STALE_PROCESSING_AGE:
            logging.info("Somalier extract %s for VCF %s still processing - skipping", existing.pk, vcf.pk)
            return None
        existing.delete()
    return SomalierVCFExtract.objects.create(vcf=vcf)


def _max_genotyped_sites(vcf_extract: SomalierVCFExtract) -> int:
    """ The best sample's het+hom count - below the setting there's nothing to relate or predict on """
    data = vcf_extract.somaliersampleextract_set.aggregate(m=Max(F("het_count") + F("hom_count")))
    return data["m"] or 0


@celery.shared_task
def somalier_vcf_id(vcf_id: int):
    """ Extract, Ancestry, and Relate (Cohort) - each stage independent so one failure isn't fatal """
    vcf = VCF.objects.get(pk=vcf_id)
    vcf_extract = _create_vcf_extract(vcf)
    if vcf_extract is None:
        return

    if not vcf.has_genotype:
        _set_skipped(vcf_extract, "VCF has no genotype field")
        return

    processing_dir = get_import_processing_dir(vcf_extract.pk, "somalier_vcf_extract")
    try:
        with _stage_timing("extract", vcf_id):
            _somalier_vcf_extract(vcf_extract, processing_dir)
    except Exception:
        _record_error(vcf_extract, f"extract vcf={vcf_id}")
        return
    finally:
        if settings.IMPORT_PROCESSING_DELETE_TEMP_FILES_ON_SUCCESS:
            shutil.rmtree(processing_dir, ignore_errors=True)

    genotyped_sites = _max_genotyped_sites(vcf_extract)
    too_few_sites = genotyped_sites < settings.SOMALIER["min_genotyped_sites"]
    skip_reason = f"Only {genotyped_sites} genotyped sites in the best sample"

    ancestry_run = SomalierAncestryRun.objects.create(vcf_extract=vcf_extract)
    if not settings.SOMALIER["ancestry_enabled"]:
        _set_skipped(ancestry_run, "settings.SOMALIER['ancestry_enabled'] is off")
    elif too_few_sites:
        _set_skipped(ancestry_run, skip_reason)
    else:
        try:
            with _stage_timing("ancestry", vcf_id):
                _somalier_ancestry(ancestry_run)
        except Exception:
            _record_error(ancestry_run, f"ancestry vcf={vcf_id}")

    SomalierCohortRelate.objects.filter(cohort=vcf.cohort).delete()  # Delete any previous versions
    relate = SomalierCohortRelate.objects.create(cohort=vcf.cohort, cohort_version=vcf.cohort.version)
    if too_few_sites:
        _set_skipped(relate, skip_reason)
    else:
        try:
            with _stage_timing("relate", vcf_id):
                _somalier_relate(relate)
        except Exception:
            _record_error(relate, f"cohort relate vcf={vcf_id}")


def _somalier_vcf_extract(vcf_extract: SomalierVCFExtract, processing_dir):
    cfg = SomalierConfig()
    vcf = vcf_extract.vcf
    vcf_filename = _write_somalier_vcf(cfg, processing_dir, vcf_extract)

    somalier_bin = cfg.get_annotation("command")
    extract_cmd = [somalier_bin, "extract",
                   "--out-dir", vcf_extract.get_somalier_dir(),
                   "--sites", cfg.get_sites(vcf.genome_build),
                   "--fasta", vcf.genome_build.reference_fasta,
                   vcf_filename]
    vcf_extract.execute(extract_cmd)


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
    logging.info("return_code: %s, stdout: %s, stderr: %s", return_code, stdout, stderr)
    return exported_vcf_filename


def _somalier_ancestry(ancestry_run: SomalierAncestryRun):
    cfg = SomalierConfig()
    compare_samples = os.path.join(cfg.get_annotation("ancestry_somalier_dir"), "*.somalier")
    ancestry_report_dir = Path(ancestry_run.get_report_dir())
    ancestry_report_dir.mkdir(parents=True, exist_ok=True)
    somalier_bin = cfg.get_annotation("command")
    ancestry_cmd = [somalier_bin, "ancestry", "--labels", cfg.get_annotation("ancestry_labels"),
                    compare_samples, "++", *ancestry_run.get_sample_somalier_filenames()]
    # Force Somalier to only use 1 thread - have had it run very slow with multi-core
    # https://github.com/brentp/somalier/issues/61#issuecomment-750492570
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    ancestry_run.execute(ancestry_cmd, cwd=ancestry_report_dir, env=env)

    # Use TSV to write sample specific files
    df = pd.read_csv(ancestry_report_dir / "somalier-ancestry.somalier-ancestry.tsv", sep='\t', index_col=0)
    for sample_extract in ancestry_run.vcf_extract.somaliersampleextract_set.all():
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
    if not somalier_relate.has_hom_ref_calls:
        # Somalier --unknown    set unknown genotypes to hom-ref.
        # Without 0/0 calls in the VCF an absent site means hom-ref, not unknown
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


@celery.shared_task
def somalier_cohort_relate(cohort_id: int):
    cohort = Cohort.objects.get(pk=cohort_id)
    relate = SomalierCohortRelate.objects.create(cohort=cohort, cohort_version=cohort.version)
    _somalier_relate(relate)


@celery.shared_task
def somalier_trio_relate(trio_id: int):
    trio = Trio.objects.get(pk=trio_id)
    relate = SomalierTrioRelate.objects.create(trio=trio)
    _somalier_relate(relate)


def _delete_orphan_somalier_files() -> int:
    """ A deleted Sample takes its SomalierSampleExtract with it but leaves the .somalier file the
        all-samples wildcard picks up. Returns how many were removed. """
    cfg = SomalierConfig()
    filenames = glob.glob(os.path.join(cfg["vcf_base_dir"], "**", "*.somalier"), recursive=True)
    sample_ids = set(SomalierSampleExtract.objects.values_list("sample_id", flat=True))
    num_deleted = 0
    for filename in filenames:
        sample_name = os.path.basename(filename).rsplit(".", 1)[0]
        sample_id = AbstractSomalierModel.sample_name_to_id(sample_name)
        if not sample_id.isdigit() or int(sample_id) not in sample_ids:
            logging.info("Removing somalier file for deleted sample: %s", filename)
            os.remove(filename)
            num_deleted += 1
    return num_deleted


def _all_samples_relate_needed(orphans_removed: int) -> bool:
    """ Nightly, so only worth running when an extract has changed or a sample has gone """
    last_success = SomalierAllSamplesRelate.objects.filter(status=ProcessingStatus.SUCCESS).order_by("-modified").first()
    if last_success is None:
        return True
    if SomalierVCFExtract.objects.filter(modified__gt=last_success.modified).exists():
        return True
    return bool(orphans_removed)


def _load_somalier_pairs(all_samples: SomalierAllSamplesRelate, pairs_filename: str) -> int:
    """ Replaces every stored pair with the ones from this run. Returns how many were kept """
    somalier_settings = settings.SOMALIER["relatedness"]
    df = pd.read_csv(pairs_filename, sep='\t')
    shared_het_mask = df["shared_hets"] >= somalier_settings["min_shared_hets"]
    shared_hom_mask = df["shared_hom_alts"] > somalier_settings["min_shared_hom_alts"]
    relateness_mask = df["relatedness"] > somalier_settings["min_relatedness"]
    df = df[shared_het_mask & shared_hom_mask & relateness_mask]

    pairs = []
    for _, row in df.iterrows():
        row_data = dict(row)
        row_data.pop("expected_relatedness")
        # Sample_ids are at the start
        sample_a_id = AbstractSomalierModel.sample_name_to_id(row_data.pop("#sample_a"))
        sample_b_id = AbstractSomalierModel.sample_name_to_id(row_data.pop("sample_b"))
        pairs.append(SomalierRelatePairs(relate=all_samples, sample_a_id=sample_a_id,
                                         sample_b_id=sample_b_id, **row_data))

    # A sample deleted since its VCF was extracted is still in the report - resolve them in one query
    sample_ids = {p.sample_a_id for p in pairs} | {p.sample_b_id for p in pairs}
    existing_ids = set(Sample.objects.filter(pk__in=sample_ids).values_list("pk", flat=True))
    pairs = [p for p in pairs if int(p.sample_a_id) in existing_ids and int(p.sample_b_id) in existing_ids]

    SomalierRelatePairs.objects.exclude(relate=all_samples).delete()  # unique on (sample_a, sample_b)
    SomalierRelatePairs.objects.bulk_create(pairs, batch_size=2000)
    return len(pairs)


@celery.shared_task
def somalier_all_samples(force: bool = False):
    """ All-vs-all relate across every successfully imported sample (#393) - beat scheduled nightly """
    orphans_removed = _delete_orphan_somalier_files()
    all_samples = SomalierAllSamplesRelate.objects.create(status=ProcessingStatus.PROCESSING)
    if not (force or _all_samples_relate_needed(orphans_removed)):
        _set_skipped(all_samples, "No VCF extract has changed since the last successful run")
        return
    try:
        start = time.time()
        related_dir = _somalier_relate(all_samples)
        pairs_filename = os.path.join(related_dir, "somalier.pairs.tsv")
        num_pairs = _load_somalier_pairs(all_samples, pairs_filename)
        logging.info("somalier all samples relate: %d pairs took %.1fs", num_pairs, time.time() - start)
        all_samples.status = ProcessingStatus.SUCCESS
        all_samples.save()
        # Their report dirs go with them (@see somalier_relate_pre_delete_handler)
        SomalierAllSamplesRelate.objects.exclude(pk=all_samples.pk).delete()
    except Exception:
        tb = get_traceback()
        logging.error(tb)
        all_samples.error_exception = tb
        all_samples.status = ProcessingStatus.ERROR
        all_samples.save()
