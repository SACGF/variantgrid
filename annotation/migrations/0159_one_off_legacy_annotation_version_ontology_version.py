from django.db import migrations


def _one_off_legacy_annotation_version_ontology_version(apps, _schema_editor):
    """ AnnotationVersions created before ontology versioning have no OntologyVersion, which makes them invalid
        (breaking any analysis using them). 0056_one_off_add_ontology_version only fixed deployments that already
        had an OntologyVersion at that point - the rest were left null. """

    AnnotationVersion = apps.get_model("annotation", "AnnotationVersion")
    OntologyVersion = apps.get_model("ontology", "OntologyVersion")

    legacy_ontology_version = OntologyVersion.objects.order_by("pk").first()
    if legacy_ontology_version is None:
        return  # Ontology has never been imported here, nothing to assign

    records = []
    for av in AnnotationVersion.objects.filter(ontology_version__isnull=True):
        ontology_version_id = legacy_ontology_version.pk
        if gav := av.gene_annotation_version:
            # Keep it consistent with the gene annotation, otherwise validate_gene_annotation fails
            if gav.ontology_version_id:
                ontology_version_id = gav.ontology_version_id
        av.ontology_version_id = ontology_version_id
        records.append(av)

    if records:
        print(f"Assigning OntologyVersion to {len(records)} legacy AnnotationVersions")
        # bulk_update so we don't touch the auto_now annotation_date (which decides which version is latest)
        AnnotationVersion.objects.bulk_update(records, ["ontology_version"], batch_size=200)


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0158_alter_variantannotation_variant_class'),
        ('ontology', '0027_alter_ontologyterm_ontology_service'),
    ]

    operations = [
        migrations.RunPython(_one_off_legacy_annotation_version_ontology_version),
    ]
