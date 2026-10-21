from django.db import migrations

DISCORDANCE_REPORT_QUESTIONS = [
    # (order, key, heading, label) - as configured on Shariant
    (0, "curation_method_different", "Curation Method", "Different curation guideline used for interpretation"),
    (1, "curation_method_same", "Curation Method", "Difference in the implementation of the same curation guideline"),
    (2, "curation_method_terminology", "Curation Method", "Differences in use of classification terminology"),
    (3, "curation_method_haplotype", "Curation Method", "Classification of variant as a haplotype"),
    (4, "curation_method_error", "Curation Method", "Curation amended after internal review"),
    (10, "evidence_public_new", "Evidence", "New and/or different publicly available evidence for curation"),
    (11, "evidence_paywall", "Evidence", "Use of data behind a paywall"),
    (12, "evidence_lab", "Evidence", "Use of internal laboratory data"),
    (13, "evidence_weighting", "Evidence", "Difference in the interpretation and/or weighting of the same curation evidence"),
    (19, "clinical_contexts_condition", "Clinical Contexts", "Interpreted against different conditions"),
    (20, "clinical_contexts", "Clinical Contexts", "Interpreted in different clinical contexts"),
    (40, "other", "Other", "Other"),
]


def _create_discordance_report_topic(apps, _schema_editor):
    """ Overlap and DiscordanceReport reviews hard-code this topic key. A deployment that already created it in the
        admin keeps its own name and questions """
    ReviewTopic = apps.get_model("review", "ReviewTopic")
    ReviewQuestion = apps.get_model("review", "ReviewQuestion")

    topic, created = ReviewTopic.objects.get_or_create(
        key="discordance_report",
        defaults={"name": "Discordance Resolution Minutes",
                  "heading": "What were the differences that caused the discordance"},
    )
    if created or not ReviewQuestion.objects.filter(topic=topic).exists():
        for order, key, heading, label in DISCORDANCE_REPORT_QUESTIONS:
            ReviewQuestion.objects.get_or_create(
                key=key,
                defaults={"topic": topic, "heading": heading, "label": label, "order": order},
            )


class Migration(migrations.Migration):

    dependencies = [
        ('review', '0004_review_is_complete'),
    ]

    operations = [
        migrations.RunPython(_create_discordance_report_topic, migrations.RunPython.noop),
    ]
