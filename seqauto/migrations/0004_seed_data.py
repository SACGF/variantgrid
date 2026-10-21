from django.db import migrations

from library.django_utils.seed_data import insert_seed_data

QC_TYPES = [
    {'name': 'ExecSummaryQC', 'qc_object_path': 'qcexecsummary'},
    {'name': 'Flagstats', 'qc_object_path': 'alignment_files__flagstats', 'total_field': 'total'},
    {'name': 'IlluminaFlowcellQC', 'qc_object_path': 'sequencing_sample__sample_sheet__illuminaflowcellqc'},
    {'name': 'Read_1_FastQC', 'qc_object_path': 'sequencing_sample__unalignedreads__fastq_r1__fastqc'},
    {'name': 'Read_2_FastQC', 'qc_object_path': 'sequencing_sample__unalignedreads__fastq_r2__fastqc'},
]

QC_COLUMNS = [
    {'qc_type_id': 'IlluminaFlowcellQC', 'name': 'Mean Cluster Density', 'field': 'mean_cluster_density'},
    {'qc_type_id': 'IlluminaFlowcellQC', 'name': 'Mean PF Cluster Density', 'field': 'mean_pf_cluster_density'},
    {'qc_type_id': 'IlluminaFlowcellQC', 'name': 'Total Clusters', 'field': 'total_clusters'},
    {'qc_type_id': 'IlluminaFlowcellQC', 'name': 'Total PF Clusters', 'field': 'total_pf_clusters'},
    {'qc_type_id': 'IlluminaFlowcellQC', 'name': 'Percentage of Clusters PF', 'field': 'percentage_of_clusters_pf'},
    {'qc_type_id': 'IlluminaFlowcellQC', 'name': 'Aligned To PhiX', 'field': 'aligned_to_phix'},
    {'qc_type_id': 'Read_1_FastQC', 'name': 'Total Sequences', 'field': 'total_sequences'},
    {'qc_type_id': 'Read_1_FastQC', 'name': 'Filtered Sequences', 'field': 'filtered_sequences'},
    {'qc_type_id': 'Read_1_FastQC', 'name': 'GC %', 'field': 'gc'},
    {'qc_type_id': 'Read_2_FastQC', 'name': 'Total Sequences', 'field': 'total_sequences'},
    {'qc_type_id': 'Read_2_FastQC', 'name': 'Filtered Sequences', 'field': 'filtered_sequences'},
    {'qc_type_id': 'Read_2_FastQC', 'name': 'GC %', 'field': 'gc'},
    {'qc_type_id': 'Flagstats', 'name': 'Total', 'field': 'total'},
    {'qc_type_id': 'Flagstats', 'name': 'Read1', 'field': 'read1'},
    {'qc_type_id': 'Flagstats', 'name': 'Read2', 'field': 'read2'},
    {'qc_type_id': 'Flagstats', 'name': 'Mapped', 'field': 'mapped'},
    {'qc_type_id': 'Flagstats', 'name': 'Properly Paired', 'field': 'properly_paired'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Deduplicated Reads', 'field': 'deduplicated_reads'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Indels in dbSNP (%)', 'field': 'indels_dbsnp_percent'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Mean coverage across genes', 'field': 'mean_coverage_across_genes'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Mean coverage across kit', 'field': 'mean_coverage_across_kit'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Median insert', 'field': 'median_insert'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Number of indels', 'field': 'number_indels'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Number SNPs', 'field': 'number_snps'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Percent of bases >=10x in GOI', 'field': 'percent_10x_goi'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Percent of bases >=20x in GOI', 'field': 'percent_20x_goi'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Percent of bases >=20x in Kit', 'field': 'percent_20x_kit'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Duplication (%)', 'field': 'percent_duplication'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Error Rate (%)', 'field': 'percent_error_rate'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Map to diff chrom (%)', 'field': 'percent_map_to_diff_chr'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Read enrichment (%)', 'field': 'percent_read_enrichment'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Reads (%)', 'field': 'percent_reads'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Softclip (%)', 'field': 'percent_softclip'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Reads', 'field': 'reads'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'SampleID LOD', 'field': 'sample_id_lod'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'SNPs in dbSNP (%)', 'field': 'snp_dbsnp_percent'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Ts/TV ratio', 'field': 'ts_to_tv_ratio'},
    {'qc_type_id': 'ExecSummaryQC', 'name': 'Uniformity of coverage', 'field': 'uniformity_of_coverage'},
]

SEQUENCER_MODELS = [
    {'model': 'HiSeq 2500', 'manufacturer_id': 1, 'data_naming_convention': 'H'},
    {'model': 'MiSeq', 'manufacturer_id': 1, 'data_naming_convention': 'M'},
    {'model': 'NextSeq', 'manufacturer_id': 1, 'data_naming_convention': 'M'},
]


def _insert_seed_data(apps, schema_editor):
    insert_seed_data(apps, schema_editor, [
        ("seqauto.QCType", QC_TYPES),
        ("seqauto.QCColumn", QC_COLUMNS),
        ("seqauto.SequencerModel", SEQUENCER_MODELS),
    ])


class Migration(migrations.Migration):

    initial = True  # manual.squashed_history records it with the initials on a reset database

    dependencies = [
        ('seqauto', '0003_initial'),
        ('snpdb', '0002_seed_data'),
    ]

    operations = [
        migrations.RunPython(_insert_seed_data),
    ]
