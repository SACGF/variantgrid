from django.db import migrations


def _splicegirl_vcf_source_settings(apps, _schema_editor):
    """ SpliceGirl (Illumina TSO 500 RNA splice calls) reuses AD and DP for splice-specific meanings:
        AD is Number=1 splice-supporting reads and DP is *reference* reads, not total depth. Binding
        them by name gives empty allele depth, no VAF, and a read depth column showing reference counts.

        With ref and alt set explicitly, get_ref_alt_allele_depth_function takes its ref/alt branch and
        VAF is derived as AD / (AD + DP) - which in this caller's files is exactly
        INFO/ALTDEDUP / (ALTDEDUP + REFDEDUP). """

    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.update_or_create(
        source_regex="^SpliceGirl",
        defaults={
            "sample_field_overrides": {
                "allele_depth_field": None,  # AD isn't the packed ref,alt array cyvcf2 expects
                "alt_depth_field": "AD",  # splice-supporting reads
                "ref_depth_field": "DP",  # reference reads
                "read_depth_field": None,  # no total-depth field in this caller
                "allele_frequency_field": None,  # so VAF is derived rather than read
            },
        },
    )


def _remove_splicegirl_vcf_source_settings(apps, _schema_editor):
    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.filter(source_regex="^SpliceGirl").delete()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0203_vcfsourcesettings_sample_field_overrides'),
    ]

    operations = [
        migrations.RunPython(_splicegirl_vcf_source_settings, _remove_splicegirl_vcf_source_settings)
    ]
