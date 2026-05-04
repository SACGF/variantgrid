from collections import defaultdict

from django.db import migrations, models
from django.db.models import Q


def _migrate_active_to_status(apps, schema_editor):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")

    in_progress_per_build = defaultdict(list)
    for vav in VariantAnnotationVersion.objects.filter(active=True):
        vav.status = "ACTIVE"
        vav.save(update_fields=["status"])

    for vav in VariantAnnotationVersion.objects.filter(active=False).order_by("genome_build_id", "pk"):
        in_progress_per_build[vav.genome_build_id].append(vav)

    for genome_build_id, vavs in in_progress_per_build.items():
        has_active = VariantAnnotationVersion.objects.filter(
            genome_build_id=genome_build_id, status="ACTIVE"
        ).exists()

        if has_active:
            # Any inactive rows alongside an ACTIVE one are HISTORICAL
            for vav in vavs:
                vav.status = "HISTORICAL"
                vav.save(update_fields=["status"])
        else:
            if len(vavs) > 1:
                pks = [v.pk for v in vavs]
                raise RuntimeError(
                    f"Cannot migrate VariantAnnotationVersion: genome_build_id={genome_build_id} has "
                    f"multiple in-progress (active=False) rows with no ACTIVE: pks={pks}. "
                    f"Resolve manually before running this migration."
                )
            vavs[0].status = "NEW"
            vavs[0].save(update_fields=["status"])


def _reverse_status_to_active(apps, schema_editor):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    for vav in VariantAnnotationVersion.objects.all():
        vav.active = (vav.status == "ACTIVE")
        vav.save(update_fields=["active"])


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0142_annotsv_more_fields"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotationversion",
            name="status",
            field=models.CharField(
                choices=[("NEW", "New"), ("ACTIVE", "Active"), ("HISTORICAL", "Historical")],
                default="NEW",
                max_length=10,
            ),
        ),
        migrations.RunPython(_migrate_active_to_status, _reverse_status_to_active),
        migrations.RemoveField(
            model_name="variantannotationversion",
            name="active",
        ),
        migrations.AddConstraint(
            model_name="variantannotationversion",
            constraint=models.UniqueConstraint(
                fields=("genome_build",),
                condition=Q(status="ACTIVE"),
                name="one_active_vav_per_build",
            ),
        ),
    ]
