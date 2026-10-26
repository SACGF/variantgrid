"""
Rows a fresh database starts with, loaded by each app's `seed_data` migration (#736).

The seed migrations hold the final state of every table a migrated-empty database had rows in, as
`("app_label.Model", [row, ...])` pairs: a row names only the fields that differ from the model default.
Entry point: `insert_seed_data`, called from `RunPython` with the historical apps.
"""
from django.core.management.color import no_style


def insert_seed_data(apps, schema_editor, seed_data: list[tuple[str, list[dict]]]):
    """ Bulk inserts each model's rows in order, then moves the pk sequences past any explicit pks """
    models = []
    for model_label, rows in seed_data:
        model = apps.get_model(model_label)
        if model._meta.parents:  # bulk_create can't write a multi-table child's parent row
            for row in rows:
                model(**row).save()
        else:
            model.objects.bulk_create([model(**row) for row in rows])
        models.append(model)
        models.extend(model._meta.parents)

    connection = schema_editor.connection
    with connection.cursor() as cursor:
        for sql in connection.ops.sequence_reset_sql(no_style(), models):
            cursor.execute(sql)

