import logging
from contextlib import contextmanager

from django.conf import settings
from django.db import connection, models, transaction
from django.db.utils import ProgrammingError
from django.utils.text import slugify

from library.log_utils import log_traceback
from library.utils import double_quote, single_quote
from library.utils.database_utils import run_sql


def _get_id_sequence(base_table_name: str) -> str:
    """ Postgres only names the sequence behind an identity/serial column '<table>_id_seq' if that name
        happened to be free when the column was created - a table that was rebuilt or renamed while the
        old sequence was still around gets '<table>_id_seq1' instead, so ask for the real name. """

    with connection.cursor() as cursor:
        cursor.execute("SELECT pg_get_serial_sequence(%s, 'id')", [base_table_name])
        sequence_name = cursor.fetchone()[0]

    if sequence_name is None:
        msg = f"'{base_table_name}.id' has no identity/serial sequence, so partition inserts would have no default"
        raise ValueError(msg)
    return sequence_name


def _clear_cached_cols(meta):
    for field in meta.concrete_fields:
        field.__dict__.pop("cached_col", None)


@contextmanager
def temporary_db_table(model, db_table: str):
    """ Point a model at one of its partition tables for the duration of a query.

        Restoring _meta.db_table is not enough on its own: Field.cached_col is a cached_property
        bound to whatever the table was called the first time that column was rendered, so a query
        run while swapped leaves the model's fields pointing at the partition forever after. Later
        queries then render a mix of both table names ('missing FROM-clause entry for table ...'),
        which outlives the swap for the life of the process. Clear the cache both ways. """

    meta = model._meta
    original_db_table = meta.db_table
    _clear_cached_cols(meta)
    meta.db_table = db_table
    try:
        yield
    finally:
        meta.db_table = original_db_table
        _clear_cached_cols(meta)


class RelatedModelsPartitionModel(models.Model):
    """ Partitions related model records by FK
        This should be inherited by the object that HOLDS the collection (that records point to)

        @see https://github.com/SACGF/variantgrid/wiki/Data-Partitioning """
    RECORDS_BASE_TABLE_NAMES = []
    RECORDS_FK_FIELD_TO_THIS_MODEL = None  # FK that points to this
    PARTITION_LABEL_TEXT = None  # Used to join between base_table_name and pk

    class Meta:
        abstract = True

    def save(self, force_insert=False, force_update=False, using=None, update_fields=None):
        created = not self.pk
        super().save(force_insert=force_insert, force_update=force_update, using=using, update_fields=update_fields)
        if created:
            self.create_partition()

    def create_partition(self):
        for base_table_name in self.RECORDS_BASE_TABLE_NAMES:
            self.create_partition_for_base_table(base_table_name)

    def create_partition_for_base_table(self, base_table_name):
        sql_template = """
    CREATE TABLE "%(table_name)s" (
        LIKE %(base_table_name)s including indexes,
        CHECK (%(records_fk_field)s = %(pk)s)
    ) INHERITS (%(base_table_name)s);

    -- If a column in the parent table is an identity column, that property is not inherited
    -- @see https://www.postgresql.org/docs/current/sql-createtable.html
    
    ALTER TABLE "%(table_name)s"
    ALTER COLUMN id SET DEFAULT nextval('%(id_sequence)s');
    """

        table_name = self.get_partition_table(base_table_name=base_table_name)
        logging.info("Creating Partition '%s'", table_name)
        pk = self.pk
        if isinstance(pk, str):
            pk = single_quote(pk)
        sql = sql_template % {"base_table_name": base_table_name,
                              "table_name": table_name,
                              "records_fk_field": self.RECORDS_FK_FIELD_TO_THIS_MODEL,
                              "pk": pk,
                              "id_sequence": _get_id_sequence(base_table_name)}
        run_sql(sql)

    def get_partition_table(self, base_table_name=None):
        if self.pk is None:
            msg = "Cannot set table as model is not saved"
            raise ValueError(msg)

        if base_table_name is not None:
            if base_table_name not in self.RECORDS_BASE_TABLE_NAMES:
                msg = f"get_partition_table(base_table_name={base_table_name}) not in RECORDS_BASE_TABLE_NAMES={self.RECORDS_BASE_TABLE_NAMES}"
                raise ValueError(msg)
        else:
            # Have to keep backwards compatibility
            if len(self.RECORDS_BASE_TABLE_NAMES) == 1:
                base_table_name = self.RECORDS_BASE_TABLE_NAMES[0]
            else:
                msg = f"get_partition_table() called with no argument, and >1 tables: {self.RECORDS_BASE_TABLE_NAMES}"
                raise ValueError(msg)

        return f"{base_table_name}_{self.PARTITION_LABEL_TEXT}_{slugify(self.pk)}"

    def sql_partition_transformer(self, sql):
        """' Modifies SQL generated by QuerySet
             @see library.django_utils.django_queryset_sql_transformer.get_queryset_with_transformer_hook  """

        for base_table_name in self.RECORDS_BASE_TABLE_NAMES:
            # Quote them, otherwise things can get replaced multiple times
            quoted_base_table_name = double_quote(base_table_name)
            quoted_partition_table_name = double_quote(self.get_partition_table(base_table_name=base_table_name))
            sql = sql.replace(quoted_base_table_name, quoted_partition_table_name)
        return sql

    def delete_related_objects(self):
        self._partition_table_op("drop")

    def truncate_related_objects(self):
        self._partition_table_op("truncate")

    def _partition_table_op(self, op):
        for base_table_name in self.RECORDS_BASE_TABLE_NAMES:
            table_name = self.get_partition_table(base_table_name=base_table_name)
            if op == "drop":
                self._warn_if_no_archive(table_name)
            sql = f'{op} table "{table_name}";'
            try:
                # Savepoint so a failed statement (e.g. a missing partition) only rolls
                # back this op, rather than poisoning an enclosing transaction.atomic()
                # (e.g. the VCF archive in snpdb/archive.py).
                with transaction.atomic():
                    run_sql(sql)
            except ProgrammingError:
                if getattr(settings, "LOG_PARTITION_WARNINGS", True):
                    logging.warning(sql)
                    log_traceback(level=logging.WARNING)

    def _warn_if_no_archive(self, table_name: str):
        """ Soft warning when a partition is being dropped without an
            accompanying COMPLETE PartitionArchive row.

            Legitimate-drop callers (VCF re-import, GCC restore, archived-source
            cleanup) bypass the archive pipeline because the original source
            file is canonical for those models -- the warning is defensive
            telemetry only, not enforcement.
        """
        try:
            from snpdb.models.models_partition_archive import PartitionArchive
        except ImportError:
            return
        meta = self._meta
        try:
            has_archive = PartitionArchive.objects.filter(
                source_app_label=meta.app_label,
                source_model=meta.object_name,
                source_pk=str(self.pk),
                source_table_names__contains=[table_name],
                status=PartitionArchive.Status.COMPLETE,
            ).exists()
        except ProgrammingError:
            return
        if not has_archive:
            logging.warning(
                "Partition drop without archive: table=%s source=%s.%s#%s",
                table_name, meta.app_label, meta.object_name, self.pk,
            )
