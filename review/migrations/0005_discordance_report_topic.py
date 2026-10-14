from django.db import migrations

DISCORDANCE_REPORT_QUESTIONS = [
    # (key, heading, label)
    ("discordance_report_criteria", "Evidence", "Different application or weighting of ACMG criteria"),
    ("discordance_report_population", "Evidence", "Population frequency data or thresholds"),
    ("discordance_report_literature", "Evidence", "Published literature"),
    ("discordance_report_internal_data", "Evidence", "Internal or unpublished case data"),
    ("discordance_report_phenotype", "Clinical", "Patient phenotype or clinical information"),
    ("discordance_report_out_of_date", "Other", "A classification was out of date"),
    ("discordance_report_other", "Other", "Other"),
]


def _create_discordance_report_topic(apps, _schema_editor):
    """ Overlap and DiscordanceReport reviews hard-code this topic key. A deployment that already created it in the
        admin keeps its own name and questions """
    ReviewTopic = apps.get_model("review", "ReviewTopic")
    ReviewQuestion = apps.get_model("review", "ReviewQuestion")

    topic, created = ReviewTopic.objects.get_or_create(
        key="discordance_report",
        defaults={"name": "Discordance", "heading": "Reasons for the difference"},
    )
    if created or not ReviewQuestion.objects.filter(topic=topic).exists():
        for order, (key, heading, label) in enumerate(DISCORDANCE_REPORT_QUESTIONS):
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
