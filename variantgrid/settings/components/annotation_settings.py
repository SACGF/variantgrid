import os

from variantgrid.settings.components.secret_settings import get_secret
from variantgrid.settings.components.settings_paths import ANNOTATION_BASE_DIR, VARIANTGRID_REPO_REFERENCE_DIR, \
    PRIVATE_DATA_ROOT

# GeneAnnotation is only in analyses, as an optimisation to stpre e.g. per-gene ontology records.
ANNOTATION_GENE_ANNOTATION_VERSION_ENABLED = True

ANNOTATION_VEP_FAKE_VERSION = False  # Overridden in unit tests to not call VEP to get version
ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT = None  # os.path.join(BASE_DIR, "scripts", "perlbrew_runner.sh")

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
    _VARIANT_ANNOTATION_PIPELINE_STANDARD: 2000,
    _VARIANT_ANNOTATION_PIPELINE_STRUCTURAL_VARIANT: 1000,
}
# get_unannotated_count_min_max does quick queries to try and get VEP batch sizes within a range
# If it gets below min, it does a slower query to get range lock.
# The variant table is usually ~55% alt variants but may be different due to data or if you've deleted records
ANNOTATION_VEP_BATCH_MIN = 5000  # Dont' set too low due to overhead of running pipeline etc
ANNOTATION_VEP_BATCH_MAX = 50_000  # Set to None to do all in 1 job (probably want to set FORK higher)
ANNOTATION_VEP_ARGS = []
ANNOTATION_VEP_VERSION = "110"
ANNOTATION_VEP_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "VEP")
ANNOTATION_VEP_VERSION_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_code", ANNOTATION_VEP_VERSION)
ANNOTATION_VEP_CODE_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "ensembl-vep")
ANNOTATION_VEP_PLUGINS_DIR = os.path.join(ANNOTATION_VEP_VERSION_DIR, "plugins")
ANNOTATION_VEP_CACHE_DIR = os.path.join(ANNOTATION_VEP_BASE_DIR, "vep_cache")

# @see https://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_pick_order
ANNOTATION_VEP_PICK_ORDER = None
ANNOTATION_VEP_DISTANCE = 5000  # VEP --distance arg (default=5000) - how far up/downstream to assign to a transcript
ANNOTATION_VEP_COLUMNS_VERSION = 1  # 1 = original version, 2 = May 2022
ANNOTATION_VEP_SV_OVERLAP_SAME_TYPE = True  # Only 'dup' for dups, false is all SVs overlap
ANNOTATION_VEP_SV_OVERLAP_SINGLE_VALUE_METHOD = "lowest_af"  # "greatest_overlap", "lowest_af", "exact_or_lowest_af"
ANNOTATION_VEP_SV_OVERLAP_MIN_FRACTION = 0.8
ANNOTATION_VEP_SV_MAX_SIZE = 10_000_000  # VEP default = 10M

ANNOTATION_MAX_BENIGN_RANKSCORE = 0.15
ANNOTATION_MIN_PATHOGENIC_RANKSCORE = 0.85

_ANNOTATION_FASTA_BASE_DIR = os.path.join(ANNOTATION_BASE_DIR, "fasta")

BUILD_GRCH37 = "GRCh37"
BUILD_GRCH38 = "GRCh38"
BUILD_T2TV2 = "T2T-CHM13v2.0"


ANNOTATION = {
    # We need separate 'reference_fasta' as cdot requires a NCBI fasta with contig_ids as the names
    # While VEP has issues with this, so has 'vep_config.fasta' see https://github.com/Ensembl/ensembl-vep/issues/1635

    BUILD_GRCH37: {
        "enabled": True,
        "annotation_consortium": "Ensembl",
        "columns_version": 3,
        "cytoband": os.path.join(VARIANTGRID_REPO_REFERENCE_DIR, "hg19", "cytoband.hg19.txt.gz"),
        "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.25_GRCh37.p13_genomic.fna.gz"),
        "reference_fasta_has_chr": False,
        "liftover": {
            BUILD_GRCH38: os.path.join(ANNOTATION_BASE_DIR,"liftover/GRCh37_to_GRCh38.chain.gz"),
            BUILD_T2TV2: os.path.join(ANNOTATION_BASE_DIR, "liftover/hg19ToHs1.over.chain.gz"),
        },

        # VEP paths are relative to ANNOTATION_VEP_BASE_DIR - worked out at runtime
        # so you can change just that variable and have everything else work
        # The names correspond to VEPPlugin or VEPCustom entries (but lower case)
        "vep_config": {
            "sift": True,
            "cosmic": "annotation_data/GRCh37/Cosmic_GenomeScreensMutant_v99_GRCh37.vcf.gz",
            "dbnsfp": "annotation_data/GRCh37/dbNSFP4.5a.grch37.stripped.gz",
            "dbscsnv": "annotation_data/GRCh37/dbscSNV1.1_GRCh37.txt.gz",
            "gnomad2": "annotation_data/GRCh37/gnomad2.1.1_GRCh37_combined_af.vcf.bgz",
            # We use gnomAD SV VCF with --custom twice
            "gnomad_sv": "annotation_data/GRCh37/gnomad_v2.1_sv.sites.grch37.converted.no_filters.vcf.gz",
            "gnomad_sv_name": "annotation_data/GRCh37/gnomad_v2.1_sv.sites.grch37.converted.no_filters.vcf.gz",
            # We use a VEP specific fasta due to bugs/workarounds, see https://github.com/Ensembl/ensembl-vep/issues/1635
            "fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"),
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
            "spliceai_snv": "annotation_data/GRCh37/spliceai_scores.raw.snv.hg19.vcf.gz",
            "spliceai_indel": "annotation_data/GRCh37/spliceai_scores.raw.indel.hg19.vcf.gz",
            "topmed": "annotation_data/GRCh37/TOPMED_GRCh37.vcf.gz",
            "transcript_blocklist": None,
            "uk10k": "annotation_data/GRCh37/UK10K_COHORT.20160215.sites.vcf.gz",
        }
    },
    # Only 37 is enabled by default - overwrite "enabled" in your server settings to use following builds
    BUILD_GRCH38: {
        "enabled": False,
        "annotation_consortium": "Ensembl",
        "columns_version": 3,
        "cytoband": os.path.join(VARIANTGRID_REPO_REFERENCE_DIR, "hg38", "cytoband.hg38.txt.gz"),
        "reference_fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "GCF_000001405.39_GRCh38.p13_genomic.fna.gz"),
        "reference_fasta_has_chr": False,
        "liftover": {
            BUILD_GRCH37: os.path.join(ANNOTATION_BASE_DIR, "liftover/GRCh38_to_GRCh37.chain.gz"),
            BUILD_T2TV2: os.path.join(ANNOTATION_BASE_DIR, "liftover/hg38ToHs1.over.chain.gz"),
        },

        # VEP paths are relative to ANNOTATION_VEP_BASE_DIR - worked out at runtime
        # so you can change just that variable and have everything else work
        # The names correspond to VEPPlugin or VEPCustom entries (but lower case)
        "vep_config": {
            "sift": True,
            "cosmic": "annotation_data/GRCh38/Cosmic_GenomeScreensMutant_v99_GRCh38.vcf.gz",
            "dbnsfp": "annotation_data/GRCh38/dbNSFP4.5a.grch38.stripped.gz",
            "dbscsnv": "annotation_data/GRCh38/dbscSNV1.1_GRCh38.txt.gz",
            # We use a VEP specific fasta due to bugs/workarounds, see https://github.com/Ensembl/ensembl-vep/issues/1635
            "fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "Homo_sapiens.GRCh38.dna.toplevel.fa.gz"),
            "gnomad2": "annotation_data/GRCh38/gnomad2.1.1_GRCh38_combined_af.vcf.bgz",
            "gnomad3": "annotation_data/GRCh38/gnomad3.1_GRCh38_merged.vcf.bgz",
            "gnomad4": "annotation_data/GRCh38/gnomad4.0_GRCh38_combined_af.vcf.bgz",
            # We use gnomAD SV VCF with --custom twice
            "gnomad_sv": "annotation_data/GRCh38/gnomad.v4.0.sv.merged.no_filters.vcf.gz",
            "gnomad_sv_name": "annotation_data/GRCh38/gnomad.v4.0.sv.merged.no_filters.vcf.gz",
            "mastermind": "annotation_data/GRCh38/mastermind_cited_variants_reference-2023.10.02-grch38.vcf.gz",
            "mave": "annotation_data/GRCh38/MaveDB_variants.tsv.gz",
            "maxentscan": "annotation_data/all_builds/maxentscan",
            'phastcons100way': "annotation_data/GRCh38/hg38.phastCons100way.bw",
            'phastcons46way': None,  # n/a for GRCh38
            'phastcons30way': "annotation_data/GRCh38/hg38.phastCons30way.bw",
            'phylop100way': "annotation_data/GRCh38/hg38.phyloP100way.bw",
            "phylop46way": None,  # n/a for GRCh38
            'phylop30way': "annotation_data/GRCh38/hg38.phyloP30way.bw",
            "repeatmasker": "annotation_data/GRCh38/repeatmasker_hg38.bed.gz",
            "spliceai_snv": "annotation_data/GRCh38/spliceai_scores.raw.snv.hg38.vcf.gz",
            "spliceai_indel": "annotation_data/GRCh38/spliceai_scores.raw.indel.hg38.vcf.gz",
            "topmed": "annotation_data/GRCh38/TOPMED_GRCh38_20180418.vcf.gz",
            "transcript_blocklist": "annotation_data/GRCh38/blocklist_brca1_new_transcripts.txt",
            "uk10k": "annotation_data/GRCh38/UK10K_COHORT.20160215.sites.GRCh38.vcf.gz",
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
            BUILD_GRCH37: os.path.join(ANNOTATION_BASE_DIR, "liftover/hs1ToHg19.over.chain.gz"),
            BUILD_GRCH38: os.path.join(ANNOTATION_BASE_DIR, "liftover/hs1ToHg38.over.chain.gz"),
        },

        "vep_config": {
            "cache_version": 107,  # Otherwise defaults to VEP_ANNO
            "sift": False,
            "cosmic": None,  # N/A
            "dbnsfp": None,
            "dbscsnv": None,
            "gnomad4": "annotation_data/T2T-CHM13v2.0/gnomad4.1.t2t_liftover_T2T-CHM13v2.0_combined_af.vcf.bgz",
            "gnomad_sv": "annotation_data/T2T-CHM13v2.0/gnomad.v4.0.sv.merged_t2t.no_filters.vcf.gz",
            "gnomad_sv_name": "annotation_data/T2T-CHM13v2.0/gnomad.v4.0.sv.merged_t2t.no_filters.vcf.gz",
            # We use a VEP specific fasta due to bugs/workarounds, see https://github.com/Ensembl/ensembl-vep/issues/1635
            "fasta": os.path.join(_ANNOTATION_FASTA_BASE_DIR, "Homo_sapiens-GCA_009914755.4-softmasked.fa.gz"),
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
        }
    },
}

ANNOTATION_VCF_DUMP_DIR = os.path.join(PRIVATE_DATA_ROOT, 'annotation_dump')
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
CACHED_WEB_RESOURCE_UNIPROT = "UniProt"

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
    CACHED_WEB_RESOURCE_UNIPROT,
    CACHED_WEB_RESOURCE_GENCC,
    CACHED_WEB_RESOURCE_CLINVAR_CITATIONS,
]
