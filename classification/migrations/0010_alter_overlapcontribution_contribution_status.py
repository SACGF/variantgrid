from django.db import migrations
from django.db.models import Q

import classification.enums.overlaps_enums
import library.django_utils.database_utils
from manual.operations.manual_operations import ManualOperation


def _test_overlaps_need_recalc(apps):
    """ Non-comparable contributions, or triages amending a value (possibly to a non-comparable one) """
    OverlapContribution = apps.get_model("classification", "OverlapContribution")
    return OverlapContribution.objects.filter(Q(contribution_status="Z") | Q(triage_state__status="F")).exists()


def _non_comparable_to_contributing(apps, _schema_editor):
    """ Comparability is now judged by the overlap calculator on the effective value, so a shared
        Drug Response / Risk Factor value simply contributes """
    OverlapContribution = apps.get_model("classification", "OverlapContribution")
    OverlapContribution.objects.filter(contribution_status="Z").update(contribution_status="C")


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0009_splicevault_top_events_help_text'),
        ('manual', '0001_initial'),
    ]

    operations = [
        # Overlaps where a triage amends to or from a non-comparable value change status
        ManualOperation.operation_manage(["classification_overlaps", "--recalc"],
                                         test=_test_overlaps_need_recalc),
        migrations.RunPython(_non_comparable_to_contributing, reverse_code=migrations.RunPython.noop),
        migrations.AlterField(
            model_name='overlapcontribution',
            name='contribution_status',
            field=library.django_utils.database_utils.TextFieldChoices(choices=[('P', 'Pending Calculation'), ('C', 'Contributing'), ('N', 'Not-shared'), ('X', 'No value')], choices_type=classification.enums.overlaps_enums.OverlapContributionStatus),
        ),
    ]
