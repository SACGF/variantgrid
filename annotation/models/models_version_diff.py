import operator
from collections import Counter, defaultdict
from functools import reduce

import pandas as pd
from django.db import connection
from django.db import models
from django.db.models.deletion import CASCADE
from django.utils.timesince import timesince
from model_utils.managers import InheritanceManager

from annotation.models.models import VariantAnnotationVersion
from library.utils.database_utils import dictfetchall
from snpdb.models import VariantGridColumn


class VersionDiff(models.Model):
    objects = InheritanceManager()
    num_added = models.IntegerField()
    num_modified = models.IntegerField()
    num_removed = models.IntegerField()
    num_unchanged = models.IntegerField()
    COLS_FROM_TO = []
    COLS_CHANGE_COUNT = []

    @property
    def total(self):
        return self.num_added + self.num_modified + self.num_removed + self.num_unchanged

    @property
    def percent_added(self):
        return 100.0 * self.num_added / self.total

    @property
    def percent_modified(self):
        return 100.0 * self.num_modified / self.total

    @property
    def percent_removed(self):
        return 100.0 * self.num_removed / self.total

    @property
    def percent_unchanged(self):
        return 100.0 * self.num_unchanged / self.total

    @staticmethod
    def description(version_from, version_to):
        from_date = version_from.annotation_date.date()
        to_date = version_to.annotation_date.date()
        time_between = timesince(version_to.annotation_date,
                                 version_from.annotation_date)
        return f"{from_date} to {to_date} ({time_between})"

    def get_diff_sql(self, _select_columns):
        return NotImplementedError()

    def create_version_diffs(self):
        vg_column_prefix = self.VG_COLUMN_PREFIX
        cols_from_to = self.COLS_FROM_TO
        cols_change_count = self.COLS_CHANGE_COUNT
        columns = cols_from_to + cols_change_count
        ALIASES = ["v1", "v2"]

        def get_column_alias(alias, column):
            return f'"{alias}"."{column}"'

        alias_columns = defaultdict(list)
        counters_from_to = []
        counters_change = []
        for c in columns:
            for alias in ALIASES:
                column = get_column_alias(alias, c)
                alias_columns[alias].append(column)

            variant_column = f"{vg_column_prefix}__{c}"
            vgc = VariantGridColumn.objects.get(variant_column=variant_column)
            if c in cols_from_to:
                counters_from_to.append((c, vgc, defaultdict(Counter)))
            elif c in cols_change_count:
                # In a list so that we can alter the object rather than copies
                counters_change.append((c, vgc, [0]))

        select_columns = []
        for alias in ALIASES:  # Always add primary key
            select_columns.append(get_column_alias(alias, "id"))

        for ac in alias_columns.values():
            select_columns.extend(ac)

        sql = self.get_diff_sql(select_columns)
        cursor = connection.cursor()
        cursor.execute(sql)

        self.num_added = 0
        self.num_modified = 0
        self.num_removed = 0
        self.num_unchanged = 0

        for data in dictfetchall(cursor, column_names=select_columns):
            v1_id, v2_id = [data[k] for k in select_columns[:2]]
            # Could be modified
            a1 = ALIASES[0]
            a2 = ALIASES[1]
            v1_columns = alias_columns[a1]
            v2_columns = alias_columns[a2]

            v1 = [data[v] for v in v1_columns]
            v2 = [data[v] for v in v2_columns]
            modified = reduce(operator.or_, [a != b for a, b in zip(v1, v2)])

            if v1_id and not v2_id:
                self.num_removed += 1
            elif v2_id and not v1_id:
                self.num_added += 1
            else:
                if modified:
                    self.num_modified += 1
                else:
                    self.num_unchanged += 1

            for column, vgc, counters in counters_from_to:
                val_1 = data[get_column_alias(a1, column)]
                val_2 = data[get_column_alias(a2, column)]
                counters[val_1][val_2] += 1

            for column, vgc, count_list in counters_change:
                val_1 = data[get_column_alias(a1, column)]
                val_2 = data[get_column_alias(a2, column)]
                if val_1 != val_2:
                    count_list[0] += 1

        self.save()

        for column, vgc, count_list in counters_change:
            vdr = VersionDiffChangeCountResult(version_diff=self,
                                               vg_column=vgc,
                                               count=count_list[0])
            vdr.save()

        for column, vgc, counters in counters_from_to:
            for value_from, counter in counters.items():
                for value_to, count in counter.items():
                    vdr = VersionDiffFromToResult(version_diff=self,
                                                  vg_column=vgc,
                                                  value_from=value_from,
                                                  value_to=value_to,
                                                  count=count)
                    vdr.save()

    def get_vg_column(self, base_column):
        variant_column = f"{self.VG_COLUMN_PREFIX}__{base_column}"
        return VariantGridColumn.objects.get(variant_column=variant_column)

    def get_diff_results(self):
        results = defaultdict(list)
        for c in self.COLS_FROM_TO:
            variant_column = self.get_vg_column(c)
            df = self.get_column_version_diff_as_df(variant_column)
            results["cols_from_to"].append((variant_column, df))

        for c in self.COLS_CHANGE_COUNT:
            variant_column = self.get_vg_column(c)
            count = self.get_column_version_diff_change_count(variant_column)
            results["cols_change_count"].append((variant_column, count))

        return results

    def get_column_version_diff_as_df(self, variant_grid_column):
        qs = self.versiondifffromtoresult_set.filter(vg_column=variant_grid_column)
        data = defaultdict(dict)
        for value_from, value_to, count in qs.values_list('value_from', 'value_to', 'count'):
            data[value_from][value_to] = count
        df = pd.DataFrame.from_dict(data).fillna(0)
        NO_ENTRY_DICT = {None: "No Entry"}
        return df.rename(index=NO_ENTRY_DICT, columns=NO_ENTRY_DICT)

    def get_column_version_diff_change_count(self, variant_grid_column):
        vd = self.versiondiffchangecountresult_set.get(vg_column=variant_grid_column)
        return vd.count


class VersionDiffFromToResult(models.Model):
    version_diff = models.ForeignKey(VersionDiff, on_delete=CASCADE)
    vg_column = models.ForeignKey(VariantGridColumn, on_delete=CASCADE)
    value_from = models.TextField(null=True)
    value_to = models.TextField(null=True)
    count = models.IntegerField()

    def __str__(self):
        name = f"{self.version_diff} {self.vg_column}"
        return f"{name}: {self.value_from} => {self.value_to}: {self.count}"


class VersionDiffChangeCountResult(models.Model):
    version_diff = models.ForeignKey(VersionDiff, on_delete=CASCADE)
    vg_column = models.ForeignKey(VariantGridColumn, on_delete=CASCADE)
    count = models.IntegerField()

    def __str__(self):
        return f"{self.version_diff} {self.vg_column}: changed {self.count}"


class VariantAnnotationVersionDiff(VersionDiff):
    VG_COLUMN_PREFIX = "variantannotation"
    COLS_FROM_TO = []
    COLS_CHANGE_COUNT = []

    version_from = models.ForeignKey(VariantAnnotationVersion, related_name='version_diff_from', on_delete=CASCADE)
    version_to = models.ForeignKey(VariantAnnotationVersion, related_name='version_diff_to', on_delete=CASCADE)

    def __str__(self):
        return self.description(self.version_from, self.version_to)
