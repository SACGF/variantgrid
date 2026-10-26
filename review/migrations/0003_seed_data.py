from django.db import migrations

from library.django_utils.seed_data import insert_seed_data

REVIEW_TOPICS = [
    {
        'key': 'discordance_report',
        'name': 'Discordance Resolution Minutes',
        'heading': 'What were the differences that caused the discordance',
    },
]

REVIEW_QUESTIONS = [
    {
        'key': 'clinical_contexts',
        'topic_id': 'discordance_report',
        'label': 'Interpreted in different clinical contexts',
        'heading': 'Clinical Contexts',
        'order': 20,
    },
    {
        'key': 'clinical_contexts_condition',
        'topic_id': 'discordance_report',
        'label': 'Interpreted against different conditions',
        'heading': 'Clinical Contexts',
        'order': 19,
    },
    {
        'key': 'curation_method_different',
        'topic_id': 'discordance_report',
        'label': 'Different curation guideline used for interpretation',
        'heading': 'Curation Method',
    },
    {
        'key': 'curation_method_error',
        'topic_id': 'discordance_report',
        'label': 'Curation amended after internal review',
        'heading': 'Curation Method',
        'order': 4,
    },
    {
        'key': 'curation_method_haplotype',
        'topic_id': 'discordance_report',
        'label': 'Classification of variant as a haplotype',
        'heading': 'Curation Method',
        'order': 3,
    },
    {
        'key': 'curation_method_same',
        'topic_id': 'discordance_report',
        'label': 'Difference in the implementation of the same curation guideline',
        'heading': 'Curation Method',
        'order': 1,
    },
    {
        'key': 'curation_method_terminology',
        'topic_id': 'discordance_report',
        'label': 'Differences in use of classification terminology',
        'heading': 'Curation Method',
        'order': 2,
    },
    {
        'key': 'evidence_lab',
        'topic_id': 'discordance_report',
        'label': 'Use of internal laboratory data',
        'heading': 'Evidence',
        'order': 12,
    },
    {
        'key': 'evidence_paywall',
        'topic_id': 'discordance_report',
        'label': 'Use of data behind a paywall',
        'heading': 'Evidence',
        'order': 11,
    },
    {
        'key': 'evidence_public_new',
        'topic_id': 'discordance_report',
        'label': 'New and/or different publicly available evidence for curation',
        'heading': 'Evidence',
        'order': 10,
    },
    {
        'key': 'evidence_weighting',
        'topic_id': 'discordance_report',
        'label': 'Difference in the interpretation and/or weighting of the same curation evidence',
        'heading': 'Evidence',
        'order': 13,
    },
    {'key': 'other', 'topic_id': 'discordance_report', 'label': 'Other', 'heading': 'Other', 'order': 40},
]


def _insert_seed_data(apps, schema_editor):
    insert_seed_data(apps, schema_editor, [
        ("review.ReviewTopic", REVIEW_TOPICS),
        ("review.ReviewQuestion", REVIEW_QUESTIONS),
    ])


class Migration(migrations.Migration):

    initial = True  # manual.squashed_history records it with the initials on a reset database

    dependencies = [
        ('review', '0002_initial'),
    ]

    operations = [
        migrations.RunPython(_insert_seed_data),
    ]
