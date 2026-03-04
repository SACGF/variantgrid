from django.db import migrations


# SA Path config migrated from ANALYSIS_AUTO_LAUNCH_SAMPLE_ANALYSIS_TEMPLATE_CONFIG setting
# Format: (enrichment_kit_name, sample_regex, template_name)
SAPATH_AUTO_LAUNCH_CONFIG = [
    ('idt_rhampseq_idhplus', r"\d+_MO_IDH_FULL\d+_NDC", "Singleton somatic analysis_Haem1"),
    ('idt_rhampseq_idhplus', r"^\d+_MO_IDH_DPYD1_", "Singleton somatic analysis_DPYD"),
    ('idt_rhampseq_idhplus', r'^\\d+_MO_IDH_(?!(DPYD1_|FULL\\d+_NDC))([^_]+)_', "Singleton somatic analysis_rhAmpIDHplus_GOI"),
    ('idt_rhampseq_ffpe_onco', None, "Singleton somatic analysis_rhAmpFFPE_Onco_3"),
    ('idt_brca_cancer', r"^\d+_MO_BRCAPlus_BRCA12_", "Singleton somatic analysis tBRCA"),
    ('idt_brca_cancer', r"^\d+_MO_BRCAPlus_GLI1_", "Singleton somatic analysis Glioma"),
    ("idt_haem", None, "Singleton somatic analysis_Haem (GOI restricted)"),
    ("tso500", None, "TSO500_combo"),
]


def _sapath_create_auto_launch_analysis_templates(apps, _schema_editor):
    AutoLaunchAnalysisTemplate = apps.get_model("analysis", "AutoLaunchAnalysisTemplate")
    AnalysisTemplate = apps.get_model("analysis", "AnalysisTemplate")
    EnrichmentKit = apps.get_model("seqauto", "EnrichmentKit")

    for enrichment_kit_name, sample_regex, template_name in SAPATH_AUTO_LAUNCH_CONFIG:
        enrichment_kit = EnrichmentKit.objects.filter(name=enrichment_kit_name).first()
        if enrichment_kit is None:
            print(f"Skipping: EnrichmentKit '{enrichment_kit_name}' not found")
            continue

        template = AnalysisTemplate.objects.filter(name=template_name).first()
        if template is None:
            print(f"Skipping: AnalysisTemplate '{template_name}' not found")
            continue

        obj, created = AutoLaunchAnalysisTemplate.objects.get_or_create(
            enrichment_kit=enrichment_kit,
            template=template,
            sample_regex=sample_regex,
        )
        if created:
            print(f"Created AutoLaunchAnalysisTemplate: {enrichment_kit_name} / {template_name}")
        else:
            print(f"Already exists: AutoLaunchAnalysisTemplate: {enrichment_kit_name} / {template_name}")


class Migration(migrations.Migration):

    dependencies = [
        ("analysis", "0096_autolaunchanalysistemplate"),
        ("seqauto", "0034_alter_qcgenecoverage_gene_coverage_collection"),
    ]

    operations = [
        migrations.RunPython(_sapath_create_auto_launch_analysis_templates, migrations.RunPython.noop)
    ]
