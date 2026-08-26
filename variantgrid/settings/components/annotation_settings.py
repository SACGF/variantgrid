import os

from variantgrid.settings.components.secret_settings import get_secret
from variantgrid.settings.components.settings_paths import (
    ANNOTATION_BASE_DIR,
    PRIVATE_DATA_ROOT,
    VARIANTGRID_REPO_REFERENCE_DIR,
)

# GeneAnnotation is only in analyses, as an optimisation to stpre e.g. per-gene ontology records.
ANNOTATION_GENE_ANNOTATION_VERSION_ENABLED = True

ANNOTATION_VEP_FAKE_VERSION = False  # Overridden in unit tests to not call VEP to get version
ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = None  # os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

# #1710: heap ceiling for one VEP process, applied with prlimit (util-linux). None = no limit.
# This bounds the run rather than the worker pool, so a plugin that runs away on one variant takes out
# its own run and nothing else - a pool-wide cgroup cap would have to be divided by however many VEPs
# are in flight, and could pick a healthy run as the OOM victim. Sized against the workload, not the
# box: a full standard run (50,000 variants, every plugin and custom file) holds flat at ~0.4GB, and
# the worst single variant measured is under 1GB, so this is roughly 10x headroom. Deployments whose
# data needs more can raise it - keep ANNOTATION_VEP_MEMORY_LIMIT_GB x annotation_workers concurrency
# comfortably below installed RAM.
ANNOTATION_VEP_MEMORY_LIMIT_GB = 4

# I've had VEP hang on me when running --fork so by default we run in small batches
# This causes a small amount of overhead obtaining an AnnotationRangeLock
# If you get ERROR: Forked process(es) died: read-through of cross-process communication detected
# You may want to reduce ANNOTATION_VEP_BUFFER_SIZE below
# @see https://github.com/Ensembl/ensembl-vep/issues/150
ANNOTATION_VEP_FORK = 1
_VARIANT_ANNOTATION_PIPELINE_STANDARD = "S"
_VARIANT_ANNOTATION_PIPELINE_STRUCTURAL_VARIANT = "C"

ANNOTATION_VEP_BUFFER_SIZE = {
    # Default VEP is 5k but this has crashed out a 16G machine...
    # VEP holds the transcript cache for every 1Mb region the current buffer spans (it drops the rest
    # in AnnotationSource::clean_cache), so peak RSS tracks regions-per-buffer, not variant count.
    # Runtime was flat across every size measured; output byte-identical. Matters most when
    # annotation_workers runs several VEPs concurrently.
    #
    # 12.5k small variants: 4000 -> 1295MB, 2000 -> 1260MB, 1000 -> 1049MB, 500 -> 871MB, 250 -> 858MB.
    # Flattens at 500 for *dense* variants - 500 consecutive ones sit in a handful of regions, so
    # regions-per-buffer was never the dominant term there. A sparse range inverts that: 500 consecutive
    # variants scattered across a build can touch hundreds of distinct regions, which is how #1710 reached
    # 5.3GB and 7.3GB on a 6,250 variant run. 100 costs nothing to hold it down - AnnotationSource::
    # clean_cache keeps the current buffer's regions, and we dump position-sorted, so a region is loaded
    # once per run whatever the buffer size. Keep well above 2 x ANNOTATION_VEP_FORK: VEP divides the
    # buffer between forks (Runner::_forked_buffer_to_output).
    _VARIANT_ANNOTATION_PIPELINE_STANDARD: 100,
    # SVs span whole regions each, so they load far more cache per variant and keep scaling down:
    # 1000 DELs (median 49kb) -> 1774MB @1000, 1212MB @500, 752MB @250, 505MB @100, matching the
    # worst-case regions in one buffer (247 / 126 / 67 / 29). 250 keeps enough per buffer to avoid
    # re-loading shared regions between consecutive buffers.
    _VARIANT_ANNOTATION_PIPELINE_STRUCTURAL_VARIANT: 250,
}
# get_unannotated_count_min_max does quick queries to try and get VEP batch sizes within a range
# If it gets below min, it does a slower query to get range lock.
# The variant table is usually ~55% alt variants but may be different due to data or if you've deleted records
ANNOTATION_VEP_BATCH_MIN = 5000  # Dont' set too low due to overhead of running pipeline etc
ANNOTATION_VEP_BATCH_MAX = 25_000  # Set to None to do all in 1 job (probably want to set FORK higher)

# Annotation dispatcher (#2667). The dispatcher launches at most as many runs as there are free
# annotation_workers slots, keeping the rest as pending DB state it can merge into bigger batches.
# Worker-slot count is derived live from the annotation_workers pool (annotation.celery_utils);
# this fallback is only used when celery inspection sees no workers (eg headless/cron, no pool up).
ANNOTATION_WORKER_SLOTS_FALLBACK = 2
# The annotation DB-upload phase (import_annotation_run) runs on db_workers, not the throttled
# annotation_workers (VEP) pool - #1649. This is its own dispatch budget, sized to match the
# db_workers pool (--concurrency=4 in the reference config) so imports can use the full pool and
# drain fast. Deployments that raise db_workers --concurrency should raise this to match.
ANNOTATION_UPLOAD_WORKER_SLOTS = 4
# Lease reclaims (dead-worker re-dispatch) allowed before a run is failed to ERROR.
ANNOTATION_MAX_RUN_ATTEMPTS = 3
# Lease window for dead-worker reclaim. A live run renews its own lease on a background heartbeat
# (AnnotationRunLeaseHeartbeat in annotate_variants) every ANNOTATION_RUN_LEASE_HEARTBEAT_SECONDS,
# so the window no longer has to exceed the worst-case run time - a long structural-variant VEP run
# (no SV-count cap, can run many hours) stays leased by heartbeating, not by a generous static window.
# So we keep it short: a genuinely dead/SIGKILLed worker stops heartbeating and its run is reclaimed
# within one window (~15 min) instead of hours. The heartbeat interval is well under the window
# (~7 beats/window) so a transient DB/process stall doesn't falsely reclaim a live run.
ANNOTATION_RUN_LEASE_SECONDS = 900
ANNOTATION_RUN_LEASE_HEARTBEAT_SECONDS = 120
# Lease window for the count lane (#1646): the dispatcher leases a batch of uncounted runs to one
# count_annotation_runs task (~15s of work for a full batch, no heartbeat). Kept well under
# ANNOTATION_RUN_LEASE_SECONDS because a live count lease briefly holds its runs out of the VEP lane -
# a dead count worker should release them quickly. Expiry is harmless: reclaim just clears the lease
# (a count lease is never a run attempt).
ANNOTATION_COUNT_LEASE_SECONDS = 300
ANNOTATION_VEP_ARGS = []
ANNOTATION_VEP_VERSION = "116"
ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "VEP")
ANNOTATION_VEP_VERSION_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_code", ANNOTATION_VEP_VERSION)
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "ensembl-vep")
ANNOTATION_VEP_PLUGINS_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "plugins")
ANNOTATION_VEP_CACHE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_cache")

# @see https://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_pick_order
ANNOTATION_VEP_PICK_ORDER = None
ANNOTATION_VEP_DISTANCE = 5000  # VEP --distance arg (default=5000) - how far up/downstream to assign to a transcript
# Read at VAV creation only and snapshotted onto the VariantAnnotationVersion. Existing VAVs are unaffected by changes here.
# Values: "primary" -> --gencode_primary, "basic" -> --gencode_basic, None -> full Ensembl set. Ignored for RefSeq VAVs.
ANNOTATION_VEP_ENSEMBL_GENCODE = "primary"
ANNOTATION_VEP_SV_OVERLAP_SAME_TYPE = True  # Only 'dup' for dups, false is all SVs overlap
ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD = "lowest_af"  # "greatest_overlap", "lowest_af", "exact_or_lowest_af"
ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION = 0.8
ANNOTATION_VEP_SV_MAX_SIZE = 10_000_000  # VEP default = 10M

# Use pyBigWig as optimisation rather than VEP --custom (see #1657)
ANNOTATION_VEP_SV_CONSERVATION_PYBIGWIG_ENABLED = True

ANNOTATION_MAX_BENIGN_RANKSCORE = 0.15
ANNOTATION_MIN_PATHOGENIC_RANKSCORE = 0.85

# dbNSFP rankscores are legacy (replaced by raw scores at columns_version >= 4). New deployments hide
# them so nobody filters/views by them going forward. Deployments that previously used rankscores set
# this True (see env/vgaws.py, env/vgtest2.py). A value that was already set is always shown/applied
# regardless of this flag.
ANNOTATION_SHOW_LEGACY_RANKSCORES = False

_ANNOTATION_FASTA_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "fasta")

BUILD_GRCH37 = "GRCh37"
BUILD_GRCH38 = "GRCh38"
BUILD_T2TV2 = "T2T-CHM13v2.0"


ANNOTATION = {
    # 'reference_fasta' is an NCBI fasta, with contig accessions as the sequence names - required by cdot,
    # and used for VEP. The 'liftover' fasta has sequences named by chromosome, to match the contig names in
    # the chain files - it's only needed when LIFTOVER_BCFTOOLS_ENABLED. bcftools +liftover has no
    # --rename-chrs equivalent, so using 'reference_fasta' here would mean rewriting the chain contigs (see #1373)

    BUILD_GRCH37: {
        "enabled": True,
        "annotation_consortium": "RefSeq",
        "columns_version": 5,
        "cytoband": os.path.join(VARIANTGRID_REPO_REFERENCE_DIR, "hg19", "cytoband.hg19.txt.gz"),
        "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.25_GRCh37.p13_genomic.fna.gz"),
        "reference_fasta_has_chr": False,
        "liftover": {
            "fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"),
            "chain": {
                BUILD_GRCH38: os.path.join(ANNOTATION_BASE_DIR, "liftover/GRCh37_to_GRCh38.chain.gz"),
                BUILD_T2TV2: os.path.join(ANNOTATION_BASE_DIR, "liftover/hg19ToHs1.over.chain.gz"),
            },
        },

        # VEP paths are relative to ANNOTATION_VEP_BASE_DIR - worked out at runtime
        # so you can change just that variable and have everything else work
        # The names correspond to VEPPlugin or VEPCustom entries (but lower case)
        "vep_config": {
            "sift": True,
            "cosmic": "annotation_data/GRCh37/Cosmic_GenomeScreensMutant_Normal_v101_GRCh37.vcf.gz",
            "dbnsfp": "annotation_data/GRCh37/dbNSFP5.3.1a.grch37.stripped.gz",
            "dbscsnv": "annotation_data/GRCh37/dbscSNV1.1_GRCh37.txt.gz",
            "denovo_db": "annotation_data/GRCh37/denovo-db.variants.v.1.6.1.GRCh37.vcf.gz",
            "gnomad2": "annotation_data/GRCh37/gnomad2.1.1_GRCh37_combined_af.vcf.bgz",
            # We use gnomAD SV VCF with --custom twice
            "gnomad_sv": "annotation_data/GRCh37/gnomad_v2.1_sv.sites.grch37.converted.no_filters.vcf.gz",
            "gnomad_sv_name": "annotation_data/GRCh37/gnomad_v2.1_sv.sites.grch37.converted.no_filters.vcf.gz",
            "mastermind": "annotation_data/GRCh37/mastermind_cited_variants_reference-2023.10.02-grch37.vcf.gz",
            "mave": None,  # n/a for GRCh37
            "maxentscan": "annotation_data/all_builds/maxentscan",
            'phastcons100way': "annotation_data/GRCh37/hg19.100way.phastCons.bw",
            'phastcons46way': "annotation_data/GRCh37/hg19.phastCons46way.placental.bw",
            'phastcons30way': None,  # n/a for GRCh37
            'phylop100way': "annotation_data/GRCh37/hg19.100way.phyloP100way.bw",
            'phylop46way': "annotation_data/GRCh37/hg19.phyloP46way.placental.bw",
            'phylop30way': None,  # n/a for GRCh37
            "repeatmasker": "annotation_data/GRCh37/repeatmasker_hg19.bed.gz",
            "spliceai_snv": "annotation_data/GRCh37/spliceai_scores.masked.snv.hg19.vcf.gz",
            "spliceai_indel": "annotation_data/GRCh37/spliceai_scores.masked.indel.hg19.vcf.gz",
            "topmed": "annotation_data/GRCh37/TOPMED_GRCh37.vcf.gz",
            "transcript_blocklist": None,
            "uk10k": "annotation_data/GRCh37/UK10K_COHORT.20160215.sites.vcf.gz",
            # columns_version 5 plugins (#1638) - see pin_annotation_to_columns_version_4() to disable
            "protvar": "annotation_data/all_builds/ProtVar_data.db",
            "open_targets": None,  # n/a for GRCh37 (Open Targets data is GRCh38)
            "eve": None,  # n/a for GRCh37 (GRCh38 only)
            "popeve": None,  # n/a for GRCh37 (GRCh38 only)
            "promoter_ai": None,  # n/a for GRCh37 (GRCh38 only)
        }
    },
    BUILD_GRCH38: {
        # Only 37 is enabled by default - overwrite "enabled" in your server settings to use following builds
        "enabled": False,
        "annotation_consortium": "RefSeq",
        "columns_version": 5,
        "cytoband": os.path.join(VARIANTGRID_REPO_REFERENCE_DIR, "hg38", "cytoband.hg38.txt.gz"),
        "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.39_GRCh38.p13_genomic.fna.gz"),
        "reference_fasta_has_chr": False,
        "liftover": {
            "fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "Homo_sapiens.GRCh38.dna.toplevel.fa.gz"),
            "chain": {
                BUILD_GRCH37: os.path.join(ANNOTATION_BASE_DIR, "liftover/GRCh38_to_GRCh37.chain.gz"),
                BUILD_T2TV2: os.path.join(ANNOTATION_BASE_DIR, "liftover/hg38ToHs1.over.chain.gz"),
            },
        },

        # VEP paths are relative to ANNOTATION_VEP_BASE_DIR - worked out at runtime
        # so you can change just that variable and have everything else work
        # The names correspond to VEPPlugin or VEPCustom entries (but lower case)
        "vep_config": {
            "sift": True,
            "cosmic": "annotation_data/GRCh38/Cosmic_GenomeScreensMutant_Normal_v101_GRCh38.vcf.gz",
            "dbnsfp": "annotation_data/GRCh38/dbNSFP5.3.1a.grch38.stripped.gz",
            "dbscsnv": "annotation_data/GRCh38/dbscSNV1.1_GRCh38.txt.gz",
            "denovo_db": "annotation_data/GRCh38/denovo-db.variants.v.1.6.1.GRCh38.vcf.gz",
            "gnomad2": "annotation_data/GRCh38/gnomad2.1.1_GRCh38_combined_af.vcf.bgz",
            "gnomad3": "annotation_data/GRCh38/gnomad3.1_GRCh38_merged.vcf.bgz",
            "gnomad4": "annotation_data/GRCh38/gnomad4.1_GRCh38_contigs.vcf.gz",
            # We use gnomAD SV VCF with --custom twice
            "gnomad_sv": "annotation_data/GRCh38/gnomad.v4.0.sv.merged.no_filters.vcf.gz",
            "gnomad_sv_name": "annotation_data/GRCh38/gnomad.v4.0.sv.merged.no_filters.vcf.gz",
            "mastermind": "annotation_data/GRCh38/mastermind_cited_variants_reference-2023.10.02-grch38.vcf.gz",
            "mave": "annotation_data/GRCh38/MaveDB_variants_2026-04-30.stripped.tsv.gz",
            "maxentscan": "annotation_data/all_builds/maxentscan",
            'phastcons100way': "annotation_data/GRCh38/hg38.phastCons100way.bw",
            'phastcons46way': None,  # n/a for GRCh38
            'phastcons30way': "annotation_data/GRCh38/hg38.phastCons30way.bw",
            'phylop100way': "annotation_data/GRCh38/hg38.phyloP100way.bw",
            "phylop46way": None,  # n/a for GRCh38
            'phylop30way': "annotation_data/GRCh38/hg38.phyloP30way.bw",
            "repeatmasker": "annotation_data/GRCh38/repeatmasker_hg38.bed.gz",
            "spliceai_snv": "annotation_data/GRCh38/spliceai_scores.masked.snv.hg38.vcf.gz",
            "spliceai_indel": "annotation_data/GRCh38/spliceai_scores.masked.indel.hg38.vcf.gz",
            "topmed": "annotation_data/GRCh38/TOPMED_GRCh38_20180418.vcf.gz",
            "transcript_blocklist": "annotation_data/GRCh38/blocklist_brca1_new_transcripts.txt",
            "uk10k": "annotation_data/GRCh38/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz",
            # columns_version 5 plugins (#1638) - see pin_annotation_to_columns_version_4() to disable
            "protvar": "annotation_data/all_builds/ProtVar_data.db",
            "open_targets": "annotation_data/GRCh38/open_targets_26.03_vep.tsv.bgz",
            "eve": "annotation_data/GRCh38/eve_merged.vcf.gz",  # VEP >= 116
            "popeve": "annotation_data/GRCh38/grch38_popEVE_ukbb_20250715.vcf.gz",  # VEP >= 116
            "promoter_ai": "annotation_data/GRCh38/promoterAI_tss500.tsv.bgz",  # VEP >= 116
        }
    },
    BUILD_T2TV2: {
        "enabled": False,
        "annotation_consortium": "Ensembl",
        "columns_version": 3,
        "cytoband": os.path.join(VARIANTGRID_REPO_REFERENCE_DIR, "t2tv2", "chm13v2.0_cytobands_allchrs.bed.gz"),
        "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"),
        "reference_fasta_has_chr": False,
        "liftover": {
            "fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "Homo_sapiens-GCA_009914755.4-softmasked.fa.gz"),
            "chain": {
                BUILD_GRCH37: os.path.join(ANNOTATION_BASE_DIR, "liftover/hs1ToHg19.over.chain.gz"),
                BUILD_GRCH38: os.path.join(ANNOTATION_BASE_DIR, "liftover/hs1ToHg38.over.chain.gz"),
            },
        },

        "vep_config": {
            "cache_version": 107,  # Otherwise defaults to VEP_ANNO
            "sift": False,
            "cosmic": None,  # N/A
            "dbnsfp": None,
            "dbscsnv": None,
            "denovo_db": None,  # N/A
            "gnomad4": "annotation_data/T2T-CHM13v2.0/gnomad4.1.t2t_liftover_T2T-CHM13v2.0_combined_af.vcf.bgz",
            "gnomad_sv": "annotation_data/T2T-CHM13v2.0/gnomad.v4.0.sv.merged_t2t.no_filters.vcf.gz",
            "gnomad_sv_name": "annotation_data/T2T-CHM13v2.0/gnomad.v4.0.sv.merged_t2t.no_filters.vcf.gz",
            "mastermind": None,  # N/A
            "mave": None,  # N/A
            "maxentscan": "annotation_data/all_builds/maxentscan",
            'phastcons100way': None,
            'phastcons46way': None,
            'phastcons30way': None,  # n/a for GRCh37
            'phylop100way': None,
            'phylop46way': None,
            'phylop30way': None,  # n/a for GRCh37
            "repeatmasker": "annotation_data/T2T-CHM13v2.0/chm13v2.0_rmsk.bed.gz",
            "spliceai_snv": None,
            "spliceai_indel": None,
            "topmed": None,
            "transcript_blocklist": None,
            "uk10k": None,
            # columns_version 5 plugins (#1638) - inactive here (this build is pinned to columns_version 3)
            "protvar": None,
            "open_targets": None,  # n/a (Open Targets data is GRCh38)
            "eve": None,  # n/a (GRCh38 only)
            "popeve": None,  # n/a (GRCh38 only)
            "promoter_ai": None,  # n/a (GRCh38 only)
        }
    },
}

def _disable_columns_version_5_plugins(annotation):
    """Clear the columns_version 5 VEP plugin data keys (#1638) to None across all builds.

    The package default points these at real data files (ProtVar on all builds; OpenTargets / EVE /
    popEVE / PromoterAI on GRCh38). Deployments pinned below columns_version 5 haven't installed that
    plugin data (or VEP 116), so clear the keys - vep_columns.has_data_files then drops the plugin
    columns and get_vep_command skips the matching --plugin invocations.
    """
    for build_settings in annotation.values():
        vep_config = build_settings["vep_config"]
        for plugin_key in ("protvar", "open_targets", "eve", "popeve", "promoter_ai"):
            if plugin_key in vep_config:
                vep_config[plugin_key] = None


def _use_cosmic_v99(annotation):
    """Restore the COSMIC v99 VCFs the package default shipped before #1673.

    The default now points at the v101 release, whose sample count arrives in a different INFO field
    (see the cosmic_count VEPColumnDefs). Deployments pinned to an older columns_version haven't
    downloaded v101, so put them back on the release they have.
    """
    annotation[BUILD_GRCH37]["vep_config"]["cosmic"] = \
        "annotation_data/GRCh37/Cosmic_GenomeScreensMutant_v99_GRCh37.vcf.gz"
    annotation[BUILD_GRCH38]["vep_config"]["cosmic"] = \
        "annotation_data/GRCh38/Cosmic_GenomeScreensMutant_v99_GRCh38.vcf.gz"


def pin_annotation_to_columns_version_4(annotation):
    """Restore the columns_version 4 annotation config (pre-#1638, no VEP 116 plugins).

    columns_version 5 (#1638) layered VEP 116 + the ProtVar / OpenTargets / EVE / PromoterAI plugins
    on top of the v4 data (dbNSFP 5.x raw scores, masked SpliceAI, gnomAD 4.1, denovo-db, ...).
    Deployments that have the v4 data but haven't installed VEP 116 + the plugin data call this from
    their env settings file to stay on v4. The caller is also responsible for keeping
    ANNOTATION_VEP_VERSION on its historical value.
    """
    annotation[BUILD_GRCH37]["columns_version"] = 4
    annotation[BUILD_GRCH38]["columns_version"] = 4
    _disable_columns_version_5_plugins(annotation)
    _use_cosmic_v99(annotation)


def use_pre_vep112_fasta(annotation):
    """Annotate against the chromosome-named (liftover) fasta rather than the NCBI 'reference_fasta'.

    Before v112, VEP renamed contigs to match an NCBI fasta, which silently dropped plugin/custom
    annotations - see https://github.com/Ensembl/ensembl-vep/issues/1635. Deployments that keep
    ANNOTATION_VEP_VERSION below 112 call this from their env settings file.
    """
    for build_settings in annotation.values():
        if fasta := build_settings["liftover"].get("fasta"):
            build_settings["vep_config"]["fasta"] = fasta


def pin_annotation_to_columns_version_3(annotation):
    """Restore the historical (pre-#1625) columns_version 3 annotation config.

    The package default ANNOTATION now ships the latest annotation data (columns_version 5 - VEP 116
    plugins on top of dbNSFP 5.x raw scores, masked SpliceAI, gnomAD 4.1, denovo-db, ...) so a fresh
    deployment gets latest without hunting for config. Existing deployments that haven't loaded the
    newer annotation data files call this from their env settings file to stay on the previous config.
    Note the caller is also responsible for keeping ANNOTATION_VEP_VERSION on its historical value.
    """
    _disable_columns_version_5_plugins(annotation)
    _use_cosmic_v99(annotation)
    annotation[BUILD_GRCH37]["columns_version"] = 3
    annotation[BUILD_GRCH37]["vep_config"].update({
        "denovo_db": None,
        "dbnsfp": "annotation_data/GRCh37/dbNSFP4.5a.grch37.stripped.gz",
        "spliceai_snv": "annotation_data/GRCh37/spliceai_scores.raw.snv.hg19.vcf.gz",
        "spliceai_indel": "annotation_data/GRCh37/spliceai_scores.raw.indel.hg19.vcf.gz",
    })
    annotation[BUILD_GRCH38]["columns_version"] = 3
    annotation[BUILD_GRCH38]["vep_config"].update({
        "denovo_db": None,
        "dbnsfp": "annotation_data/GRCh38/dbNSFP4.5a.grch38.stripped.gz",
        "gnomad4": "annotation_data/GRCh38/gnomad4.0_GRCh38_combined_af.vcf.bgz",
        "mave": "annotation_data/GRCh38/MaveDB_variants_2023-11-29.tsv.gz",
        "spliceai_snv": "annotation_data/GRCh38/spliceai_scores.raw.snv.hg38.vcf.gz",
        "spliceai_indel": "annotation_data/GRCh38/spliceai_scores.raw.indel.hg38.vcf.gz",
    })


ANNOTATION_VCF_DUMP_DIR = os.path.join(PRIVATE_DATA_ROOT, 'annotation_dump')

# #1670: delete a run's dump/VEP/AnnotSV output once imported (failed runs always keep theirs)
ANNOTATION_DELETE_TEMP_FILES_ON_SUCCESS = True
# #1670: don't dispatch new VEP runs below this free space on ANNOTATION_VCF_DUMP_DIR (None = no limit)
ANNOTATION_MIN_FREE_DISK_GIGS = 1

# Gene-level variants (gene fusions) - @see snpdb.gene_level_variants. Only deployments importing
# fusions (TSO 500) ever create one, and a deployment with none still gets an always-empty Gene Level
# AnnotationRun per range lock, so this is off in settings files where fusions can't arrive.
ANNOTATION_GENE_LEVEL_ENABLED = True

# AnnotSV, its own annotation pipeline type since #720. Strictly opt-in: leaving this False means the
# scheduler creates no ANNOTSV AnnotationRuns at all. Deployments enable it in their per-host settings
# file once the binary + bundle are installed.
#
# Turning it on backfills, once the installed AnnotSV is registered with
# `manage.py create_new_annotation_pipeline_version --pipeline-type=A`. The scheduler creates an ANNOTSV
# run for every existing range lock on the active version, each dependent on that lock's
# STRUCTURAL_VARIANT run; the count lane finishes the (majority) SV-free ones on db_workers without ever
# dispatching them.
ANNOTATION_ANNOTSV_ENABLED = False
ANNOTATION_ANNOTSV_BIN = "/data/annotation/AnnotSV/bin/AnnotSV"
# Passed as -annotationsDir: the directory *containing* Annotations_Human (and Annotations_Exomiser),
# ie <install>/share/AnnotSV - not Annotations_Human itself
ANNOTATION_ANNOTSV_ANNOTATIONS_DIR = "/data/annotation/AnnotSV/share/AnnotSV"
# Map our genome build name -> AnnotSV's -genomeBuild value. A build absent from here is one this
# AnnotSV cannot annotate, so the runs page offers no registration for it and the scheduler skips it.
# CHM13 arrived in AnnotSV 3.5 (Annotations_Human 3.5 bundle) - deployments on an older AnnotSV drop
# that entry in their env settings file.
ANNOTATION_ANNOTSV_GENOME_BUILD = {
    BUILD_GRCH37: "GRCh37",
    BUILD_GRCH38: "GRCh38",
    BUILD_T2TV2: "CHM13",
}
# Annotations bundle has no version stamp file; admin sets this to the bundle release string they
# installed (eg "3.5.8"). Becomes AnnotationPipelineVersion.data_version, so a new bundle is registered
# and promoted the same way a new binary is - @see AnnotationPipelineVersion.
ANNOTATION_ANNOTSV_BUNDLE_VERSION = "3.5"  # Annotations_Human_3.5.tar.gz
ANNOTATION_ANNOTSV_EXTRA_ARGS: list[str] = []
ANNOTATION_ANNOTSV_TIMEOUT_SECONDS = 60 * 60
# Admin email used in PubMed queries to contact before throttling/banning
ANNOTATION_ENTREZ_EMAIL = get_secret("ENTREZ.email")  # Automatically set in in annotation.apps.AnnotationConfig
ANNOTATION_ENTREZ_API_KEY = get_secret("ENTREZ.api_key")

ANNOTATION_PUBMED_GENE_SYMBOL_COUNT_CACHE_DAYS = 30
ANNOTATION_PUBMED_SEARCH_TERMS_ENABLED = False


# These need to be able to be passed to URLS so no slashes
CACHED_WEB_RESOURCE_CLINVAR_CITATIONS = "ClinVarCitations"
CACHED_WEB_RESOURCE_GENCC = "GenCC Gene Disease Relationships"
CACHED_WEB_RESOURCE_GNOMAD_GENE_CONSTRAINT = "GnomADGeneConstraint"
CACHED_WEB_RESOURCE_HGNC = "HGNC"
CACHED_WEB_RESOURCE_LRG_REF_SEQ_GENE = "LRGRefSeqGene"
CACHED_WEB_RESOURCE_MANE = "MANE"
CACHED_WEB_RESOURCE_PANEL_APP_AUSTRALIA_PANELS = "PanelApp Australia Panels"
CACHED_WEB_RESOURCE_PANEL_APP_ENGLAND_PANELS = "Genomics England PanelApp Panels"
CACHED_WEB_RESOURCE_PFAM = "Pfam"
CACHED_WEB_RESOURCE_REFSEQ_GENE_SUMMARY = "RefSeq Gene Summary"
CACHED_WEB_RESOURCE_REFSEQ_GENE_INFO = "RefSeq Gene Info"
CACHED_WEB_RESOURCE_REFSEQ_SEQUENCE_INFO = "RefSeq Sequence Info"
CACHED_WEB_RESOURCE_REFSEQ_GENE_PUBMED_COUNTS = "RefSeq Gene PubMed Counts"
CACHED_WEB_RESOURCE_UNIPROT = "UniProt"

DISABLED_CACHED_WEB_RESOURCES = {
    CACHED_WEB_RESOURCE_PFAM: "https://github.com/SACGF/variantgrid/issues/1554",
}

ANNOTATION_CACHED_WEB_RESOURCES = [
    CACHED_WEB_RESOURCE_GNOMAD_GENE_CONSTRAINT,
    CACHED_WEB_RESOURCE_HGNC,
    CACHED_WEB_RESOURCE_LRG_REF_SEQ_GENE,
    CACHED_WEB_RESOURCE_MANE,
    CACHED_WEB_RESOURCE_PANEL_APP_AUSTRALIA_PANELS,
    CACHED_WEB_RESOURCE_PANEL_APP_ENGLAND_PANELS,
    CACHED_WEB_RESOURCE_PFAM,
    CACHED_WEB_RESOURCE_REFSEQ_GENE_SUMMARY,
    CACHED_WEB_RESOURCE_REFSEQ_GENE_INFO,
    CACHED_WEB_RESOURCE_REFSEQ_SEQUENCE_INFO,
    CACHED_WEB_RESOURCE_REFSEQ_GENE_PUBMED_COUNTS,
    CACHED_WEB_RESOURCE_UNIPROT,
    CACHED_WEB_RESOURCE_GENCC,
    CACHED_WEB_RESOURCE_CLINVAR_CITATIONS,
]
