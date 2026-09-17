import logging

from django.core.exceptions import FieldError
from django.db import migrations

from annotation.models.models_enums import ClinVarReviewStatus

# FilterNodeItem.field stores a copy of VariantGridColumn.variant_column taken when the rule was
# written. snpdb/migrations/0050_change_all_columns.py rewrote these paths but left the copies,
# so nodes holding them raise FieldError on load.
STALE_FIELDS = {
    "clinvar__clinvar_review_status": "clinvar__review_status",
    "variantannotation__uniprot__tissue_specificity":
        "variantannotation__transcript_version__gene_version__hgnc__uniprot__tissue_specificity",
}


def _log_unknown_review_status_data(FilterNodeItem):
    """ A renamed rule can still point at a review status code the enum has since changed """
    valid_values = set(ClinVarReviewStatus.values)
    for pk, data in FilterNodeItem.objects.filter(field="clinvar__review_status").values_list("pk", "data"):
        for value in data.split(","):
            if (value := value.strip()) and value not in valid_values:
                logging.warning("FilterNodeItem %d: clinvar__review_status data '%s' is not a ClinVarReviewStatus",
                                pk, value)


def _log_unresolvable_fields(FilterNodeItem, Variant):
    """ Catches paths that resolve against VariantGridColumn but point somewhere the model has renamed """
    for field in FilterNodeItem.objects.order_by("field").values_list("field", flat=True).distinct():
        try:
            Variant.objects.filter(**{f"{field}__isnull": True})
        except (FieldError, ValueError) as e:
            num_items = FilterNodeItem.objects.filter(field=field).count()
            logging.warning("FilterNodeItem.field '%s' (%d rules) doesn't resolve against Variant: %s",
                            field, num_items, e)


def _fix_stale_filter_node_item_fields(apps, _schema_editor):
    FilterNodeItem = apps.get_model("analysis", "FilterNodeItem")
    Variant = apps.get_model("snpdb", "Variant")

    for old_field, new_field in STALE_FIELDS.items():
        if num_updated := FilterNodeItem.objects.filter(field=old_field).update(field=new_field):
            logging.info("Renamed FilterNodeItem.field '%s' -> '%s' on %d rules", old_field, new_field, num_updated)

    _log_unknown_review_status_data(FilterNodeItem)
    _log_unresolvable_fields(FilterNodeItem, Variant)


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0117_one_off_hide_internal_analysis_templates'),
    ]

    operations = [
        migrations.RunPython(_fix_stale_filter_node_item_fields, migrations.RunPython.noop),
    ]
