""" Copies each description's owner and approval off its PatientTextPhenotype / CohortTextPhenotype link row, and
    deletes the descriptions nothing holds (left by deleted patients) with their sentences (#2135). A data-only
    migration, so its deferred FK checks run at its own commit rather than blocking 0023's ALTER TABLEs. """
from django.db import migrations
from django.db.models import OuterRef, Subquery

OWNER_LINK_MODELS = {"patient": "PatientTextPhenotype", "cohort": "CohortTextPhenotype"}


def _copy_owners(apps, _schema_editor):
    PhenotypeDescription = apps.get_model("patients", "PhenotypeDescription")
    for owner, link_model_name in OWNER_LINK_MODELS.items():
        LinkModel = apps.get_model("patients", link_model_name)
        link_qs = LinkModel.objects.filter(phenotype_description=OuterRef("pk"))
        description_qs = PhenotypeDescription.objects.filter(pk__in=LinkModel.objects.values("phenotype_description"))
        description_qs.update(**{owner: Subquery(link_qs.values(owner)),
                                 "approved_by": Subquery(link_qs.values("approved_by"))})


def _copy_owners_back(apps, _schema_editor):
    PhenotypeDescription = apps.get_model("patients", "PhenotypeDescription")
    for owner, link_model_name in OWNER_LINK_MODELS.items():
        LinkModel = apps.get_model("patients", link_model_name)
        owned_qs = PhenotypeDescription.objects.filter(**{f"{owner}__isnull": False})
        LinkModel.objects.bulk_create(
            LinkModel(phenotype_description_id=pk, approved_by_id=approved_by_id, **{f"{owner}_id": owner_id})
            for pk, owner_id, approved_by_id in owned_qs.values_list("pk", owner, "approved_by"))


def _delete_unowned_descriptions(apps, _schema_editor):
    """ Unowned = no patient, no cohort and nothing else pointing at it - any other app's relation in the
        migration state (an SA Path request's link rows) holds a description too """
    PhenotypeDescription = apps.get_model("patients", "PhenotypeDescription")
    held_by = {f"{relation.name}__isnull": True for relation in PhenotypeDescription._meta.related_objects
               if relation.related_model._meta.model_name != "textphenotypesentence"}
    PhenotypeDescription.objects.filter(patient__isnull=True, cohort__isnull=True, **held_by).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("patients", "0021_phenotype_description_owner_fields"),
    ]

    operations = [
        migrations.RunPython(_copy_owners, _copy_owners_back),
        migrations.RunPython(_delete_unowned_descriptions, migrations.RunPython.noop),
    ]
