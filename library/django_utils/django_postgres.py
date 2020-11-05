from django.db import models


class PostgresRealField(models.Field):
    """ 32 bit float """
    def db_type(self, connection):
        return 'real'
