import abc

from library.database_utils import run_sql
from library.jqgrid_sql import JqGridSQL, get_overrides


class AbstractGenesGrid(JqGridSQL):
    def __init__(self, user, **kwargs):
        self.fields = ["name"] + self.get_column_names()
        super().__init__(user)  # kwargs are thrown away

        #logging.info("fields: " + str(self.fields))
        column_data = [{"label": l} for l in self.get_labels()]
        overrides = get_overrides(self.fields, column_data, model_field=False)
        self.update_overrides(overrides)
        self.extra_config['sortname'] = 'name'

    def get_count(self):
        # Can't use this as we need to adjust for multiple hits for a symbol
        sql, params, _, _ = self.get_sql_params_and_columns(None)
        _, rowcount = run_sql(sql, params)
        return rowcount

    @abc.abstractmethod
    def get_column_names(self):
        pass

    @abc.abstractmethod
    def get_labels(self):
        pass
