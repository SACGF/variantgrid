from django.db import migrations

# The old columns are dropped in a migration of their own: 0120 rewrites IntersectionNode rows (and the
# GenomicInterval delete cascades SET_NULL back onto them), which leaves deferred FK trigger events pending,
# and PostgreSQL refuses to ALTER TABLE in the same transaction. IF EXISTS keeps this a no-op on databases
# that ran the original version of 0120, which dropped the columns itself.
DROP_COLUMNS_SQL = """
ALTER TABLE analysis_intersectionnode
    DROP COLUMN IF EXISTS genomic_interval_id,
    DROP COLUMN IF EXISTS hgvs_string,
    DROP COLUMN IF EXISTS hgvs_variant_id
"""

RESTORE_COLUMNS_SQL = """
ALTER TABLE analysis_intersectionnode
    ADD COLUMN genomic_interval_id integer NULL UNIQUE REFERENCES snpdb_genomicinterval (id) DEFERRABLE INITIALLY DEFERRED,
    ADD COLUMN hgvs_string text NULL,
    ADD COLUMN hgvs_variant_id integer NULL REFERENCES snpdb_variant (id) DEFERRABLE INITIALLY DEFERRED
"""


class Migration(migrations.Migration):

    dependencies = [
        ("analysis", "0127_clinvarnode_remove_no_call_pills"),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            database_operations=[
                migrations.RunSQL(DROP_COLUMNS_SQL, reverse_sql=RESTORE_COLUMNS_SQL),
            ],
            state_operations=[
                migrations.RemoveField(model_name="intersectionnode", name="genomic_interval"),
                migrations.RemoveField(model_name="intersectionnode", name="hgvs_string"),
                migrations.RemoveField(model_name="intersectionnode", name="hgvs_variant"),
            ],
        ),
    ]
