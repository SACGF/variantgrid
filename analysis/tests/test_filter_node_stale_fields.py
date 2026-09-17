from importlib import import_module

from django.core.exceptions import FieldError
from django.test import TestCase

from snpdb.models import Variant, VariantGridColumn

STALE_FIELDS = import_module("analysis.migrations.0118_one_off_fix_stale_filter_node_item_fields").STALE_FIELDS


class StaleFilterNodeItemFieldsTest(TestCase):
    """ Migration 0118 renames FilterNodeItem.field paths that snpdb/0050_change_all_columns.py
        rewrote on VariantGridColumn but not on the stored copies """

    def test_stale_paths_do_not_resolve(self):
        for stale_field in STALE_FIELDS:
            with self.assertRaises(FieldError):
                Variant.objects.filter(**{f"{stale_field}__isnull": True})

    def test_current_paths_resolve(self):
        for current_field in STALE_FIELDS.values():
            Variant.objects.filter(**{f"{current_field}__isnull": True})

    def test_current_paths_are_real_columns(self):
        for current_field in STALE_FIELDS.values():
            self.assertTrue(VariantGridColumn.objects.filter(variant_column=current_field).exists(),
                            f"{current_field} is not a VariantGridColumn.variant_column")
