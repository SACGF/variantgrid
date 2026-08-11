# Specimen.mutation_type -> Specimen.tissue_status - @see variantgrid_private#2447

from django.db import migrations, models

from manual.operations.manual_operations import ManualOperation

_MUTATION_SOMATIC = 'S'
_MUTATION_GERMLINE = 'G'
_TISSUE_STATUS_AFFECTED = 'A'
_TISSUE_STATUS_UNKNOWN = 'U'


def _mutation_type_to_tissue_status(apps, _schema_editor):
    """ Somatic means a tumour block was declared. Germline was the field default, so a 'G' can't be
        told apart from nobody touching it - it becomes Unknown and the ManualOperation below asks a
        human to look at deployments that actually populated the column """
    Specimen = apps.get_model("patients", "Specimen")
    # Everything else is already Unknown from the field default the AddField above applied
    Specimen.objects.filter(mutation_type=_MUTATION_SOMATIC).update(tissue_status=_TISSUE_STATUS_AFFECTED)

    PatientRecord = apps.get_model("patients", "PatientRecord")
    PatientRecord.objects.filter(specimen_mutation_type=_MUTATION_SOMATIC).update(
        specimen_tissue_status=_TISSUE_STATUS_AFFECTED)
    PatientRecord.objects.filter(specimen_mutation_type=_MUTATION_GERMLINE).update(
        specimen_tissue_status=_TISSUE_STATUS_UNKNOWN)


def _has_declared_germline_specimens(apps):
    """ A PatientRecord carrying 'G' is evidence the CSV column was populated rather than defaulted,
        so those specimens are worth a human deciding whether they meant Reference """
    PatientRecord = apps.get_model("patients", "PatientRecord")
    return PatientRecord.objects.filter(specimen_mutation_type=_MUTATION_GERMLINE).exists()


class Migration(migrations.Migration):

    dependencies = [
        ("patients", "0015_extraction_reference_id"),
        ("manual", "0002_deployment"),
    ]

    operations = [
        migrations.AddField(
            model_name="specimen",
            name="tissue_status",
            field=models.CharField(choices=[('R', 'Reference (unaffected)'), ('A', 'Affected / lesional'),
                                            ('U', 'Unknown')], default='U', max_length=1),
        ),
        migrations.AddField(
            model_name="patientrecord",
            name="specimen_tissue_status",
            field=models.CharField(choices=[('R', 'Reference (unaffected)'), ('A', 'Affected / lesional'),
                                            ('U', 'Unknown')], max_length=1, null=True),
        ),
        migrations.RunPython(_mutation_type_to_tissue_status, migrations.RunPython.noop),
        ManualOperation.operation_other(
            args=["Review specimens whose patient CSV declared 'Germline' - they became tissue status "
                  "Unknown, and want setting to Reference where the material really was unaffected"],
            note="Only registered where a PatientRecord carries 'G', ie the column was populated rather "
                 "than left at the old field default (variantgrid_private#2447)",
            test=_has_declared_germline_specimens),
        migrations.RemoveField(
            model_name="specimen",
            name="mutation_type",
        ),
        migrations.RemoveField(
            model_name="patientrecord",
            name="specimen_mutation_type",
        ),
    ]
