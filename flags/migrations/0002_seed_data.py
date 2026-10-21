from django.db import migrations

from library.django_utils.seed_data import insert_seed_data

FLAG_RESOLUTIONS = [
    {'id': 'closed', 'label': 'Completed', 'description': 'This is complete', 'status': 'C'},
    {'id': 'in_error', 'label': 'In Error', 'description': 'There are outstanding errors for this record'},
    {'id': 'matching_failed', 'label': 'Matching Failed', 'description': 'Could not match this record to a variant'},
    {
        'id': 'np_dont_share',
        'label': 'Do Not Share w ClinVar',
        'description': 'This record wont be considered for exporting to ClinVar',
    },
    {
        'id': 'np_share',
        'label': 'Allow Sharing w ClinVar',
        'description': 'This record will be shared with ClinVar (if normal ClinVar requirements are met)',
        'status': 'C',
    },
    {'id': 'open', 'label': 'In Progress', 'description': 'This is ongoing'},
    {'id': 'pc_applied', 'label': 'Changes Applied', 'description': 'Pending changes have been applied', 'status': 'C'},
    {
        'id': 'pc_pending_changes',
        'label': 'Pending Changes',
        'description': 'This classification has pending changes not yet reflected here',
    },
    {'id': 'pc_rejected', 'label': 'Changes Rejected', 'description': 'Pending changes were rejected', 'status': 'R'},
    {'id': 'rejected', 'label': 'Rejected', 'description': 'This is rejected', 'status': 'R'},
    {
        'id': 'sc_discordance_discussion',
        'label': 'Discordance Discussion',
        'description': 'Discussion with other labs during discordance highlighted data that changed the clinical significance',
        'status': 'C',
    },
    {
        'id': 'sc_internal_review',
        'label': 'Internal Review',
        'description': 'New data became available during an internal review such as new classifications',
        'status': 'C',
    },
    {
        'id': 'sc_new_data',
        'label': 'Summation of Data',
        'description': 'New data became available (outside of discordance discussion) such as new classifications',
        'status': 'C',
    },
    {'id': 'sc_other', 'label': 'Other', 'description': 'Other, please detail in the comments', 'status': 'C'},
    {
        'id': 'tv_accepted',
        'label': 'Change Accepted',
        'description': 'The change of transcript version is deemed acceptable',
        'status': 'C',
    },
    {
        'id': 'tv_rejected',
        'label': 'Change Rejected',
        'description': 'The change of transcript version is deemed not-acceptable',
    },
    {
        'id': 'unwithdrawn',
        'label': 'Un-Withdrawn',
        'description': 'The record was unwithdrawn and is considered again',
        'status': 'C',
    },
    {
        'id': 'vm_confirmed',
        'label': 'Variant Confirmed',
        'description': 'The variant matching has been reviewed and determined to be correct',
        'status': 'C',
    },
    {
        'id': 'vm_rejected',
        'label': 'Variant Rejected',
        'description': 'The variant matched has been determined not to be the one intended by this record',
    },
    {'id': 'withdrawn', 'label': 'Withdrawn', 'description': 'The record is withdrawn and is no longer considered'},
]

FLAG_TYPE_CONTEXTS = [
    {'id': 'allele', 'label': 'Flags for Alleles'},
    {'id': 'classification', 'label': 'Flags for Classifications'},
    {'id': 'clinical_context', 'label': 'Flags for Clinical Context Groups'},
]

FLAG_TYPES = [
    {
        'id': 'classification_clinical_context_change',
        'context_id': 'classification',
        'label': 'Clinical Grouping Changed',
        'description': 'This classification has been placed into another clinical grouping.',
        'permission': 'A',
        'comments_enabled': False,
        'importance': 2,
    },
    {
        'id': 'classification_discordant',
        'context_id': 'classification',
        'label': 'Discordance',
        'description': 'This classification is being reviewed. Data presented here may change.',
        'help_text': 'Have a look at the latest Discordance Report, from there you can initiate conversations with other labs. Be sure to raise flags for any reviews this classification is involved in.',
        'permission': 'A',
        'only_one': True,
        'comments_enabled': False,
        'importance': 2,
    },
    {
        'id': 'classification_internal_review',
        'context_id': 'classification',
        'label': 'Internal Review',
        'description': 'This classification is being internally reviewed.',
        'raise_permission': 'O',
        'permission': 'O',
        'importance': 2,
    },
    {
        'id': 'classification_not_public',
        'context_id': 'classification',
        'label': "Don't Share with ClinVar",
        'description': 'This classification can be shared to other users, but is not suitable for ClinVar or other 3rd parties',
        'raise_permission': 'O',
        'permission': 'O',
        'only_one': True,
        'importance': 2,
    },
    {
        'id': 'classification_outstanding_edits',
        'context_id': 'classification',
        'label': 'Un-submitted Edits',
        'description': 'This classification has been edited since it was last submitted. Please click on the Tick icon to submit these new changes when appropriate.',
        'permission': 'A',
        'only_one': True,
        'comments_enabled': False,
    },
    {
        'id': 'classification_pending_changes',
        'context_id': 'classification',
        'label': 'Pending Changes',
        'description': 'This classification has pending changes not yet reflected here',
        'permission': 'X',
        'only_one': True,
        'comments_enabled': False,
        'importance': 2,
    },
    {
        'id': 'classification_significance_change',
        'context_id': 'classification',
        'label': 'Assertion Change',
        'description': 'This classification record has had its classification changed since it was last published. Please document here what caused the change.',
        'help_text': 'This classification record has had its classification changed since it was last published.',
        'permission': 'O',
        'importance': 2,
    },
    {
        'id': 'classification_submitted',
        'context_id': 'classification',
        'label': 'Submitted',
        'description': "This classification is currently in review and it's clinical significance should not be considered final",
        'permission': 'A',
        'comments_enabled': False,
        'importance': 1,
    },
    {
        'id': 'classification_suggestion',
        'context_id': 'classification',
        'label': 'Suggestion',
        'description': 'Suggestion that one or more fields for this classification should be changed.',
        'help_text': """Review the suggestion made here.
If the suggestion will result in a change of data, update the record in your curation system and mark the suggestion as Closed.
If after reviewing you intend to make no change, click Reject,""",
        'raise_permission': 'U',
        'permission': 'O',
        'importance': 2,
    },
    {
        'id': 'classification_unshared',
        'context_id': 'classification',
        'label': 'Unshared Classification',
        'description': 'This classification has not yet been shared outside your organisation',
        'permission': 'A',
        'comments_enabled': False,
    },
    {
        'id': 'classification_withdrawn',
        'context_id': 'classification',
        'label': 'Classification Withdrawn',
        'description': 'This classification has been withdrawn, no guarentee is made about the quality of the data.',
        'permission': 'A',
        'only_one': True,
        'comments_enabled': False,
        'importance': 1,
    },
    {
        'id': 'clinical_context_discordance',
        'context_id': 'clinical_context',
        'label': 'Discordance',
        'description': 'This clinical group contains shared classifications that are discordant with each other.',
        'permission': 'A',
        'only_one': True,
        'importance': 2,
    },
    {
        'id': 'condition_resolution',
        'context_id': 'classification',
        'label': 'Condition Resolution',
        'description': 'This classification has been recently modified.',
        'permission': 'A',
        'comments_enabled': False,
        'importance': 1,
    },
]

FLAG_TYPE_RESOLUTIONS = [
    {'flag_type_id': 'classification_discordant', 'resolution_id': 'open'},
    {'flag_type_id': 'classification_discordant', 'resolution_id': 'closed'},
    {'flag_type_id': 'classification_suggestion', 'resolution_id': 'open'},
    {'flag_type_id': 'classification_suggestion', 'resolution_id': 'closed'},
    {'flag_type_id': 'classification_suggestion', 'resolution_id': 'rejected'},
    {'flag_type_id': 'classification_clinical_context_change', 'resolution_id': 'closed'},
    {'flag_type_id': 'classification_outstanding_edits', 'resolution_id': 'open'},
    {'flag_type_id': 'classification_outstanding_edits', 'resolution_id': 'closed'},
    {'flag_type_id': 'clinical_context_discordance', 'resolution_id': 'open'},
    {'flag_type_id': 'clinical_context_discordance', 'resolution_id': 'closed'},
    {'flag_type_id': 'classification_unshared', 'resolution_id': 'open'},
    {'flag_type_id': 'classification_unshared', 'resolution_id': 'closed'},
    {'flag_type_id': 'classification_submitted', 'resolution_id': 'closed'},
    {'flag_type_id': 'classification_internal_review', 'resolution_id': 'open'},
    {'flag_type_id': 'classification_internal_review', 'resolution_id': 'closed'},
    {'flag_type_id': 'classification_significance_change', 'resolution_id': 'open'},
    {'flag_type_id': 'classification_significance_change', 'resolution_id': 'sc_discordance_discussion'},
    {'flag_type_id': 'classification_significance_change', 'resolution_id': 'sc_new_data'},
    {'flag_type_id': 'classification_significance_change', 'resolution_id': 'sc_internal_review'},
    {'flag_type_id': 'classification_outstanding_edits', 'resolution_id': 'in_error'},
    {'flag_type_id': 'classification_withdrawn', 'resolution_id': 'withdrawn'},
    {'flag_type_id': 'classification_withdrawn', 'resolution_id': 'unwithdrawn'},
    {'flag_type_id': 'classification_not_public', 'resolution_id': 'np_dont_share'},
    {'flag_type_id': 'classification_not_public', 'resolution_id': 'np_share'},
    {'flag_type_id': 'classification_pending_changes', 'resolution_id': 'pc_pending_changes'},
    {'flag_type_id': 'classification_pending_changes', 'resolution_id': 'pc_applied'},
    {'flag_type_id': 'classification_pending_changes', 'resolution_id': 'pc_rejected'},
    {'flag_type_id': 'condition_resolution', 'resolution_id': 'closed'},
]


def _insert_seed_data(apps, schema_editor):
    insert_seed_data(apps, schema_editor, [
        ("flags.FlagResolution", FLAG_RESOLUTIONS),
        ("flags.FlagTypeContext", FLAG_TYPE_CONTEXTS),
        ("flags.FlagType", FLAG_TYPES),
        ("flags.FlagTypeResolution", FLAG_TYPE_RESOLUTIONS),
    ])


class Migration(migrations.Migration):

    initial = True  # manual.squashed_history records it with the initials on a reset database

    dependencies = [
        ('flags', '0001_initial'),
    ]

    operations = [
        migrations.RunPython(_insert_seed_data),
    ]
