from django.db import models
from django.db.models import Transform


class PostgresRealField(models.Field):
    """ 32 bit float """
    def db_type(self, connection):
        return 'real'


class MD5(Transform):
    """ This is from Django 3 """
    function = 'MD5'
    lookup_name = 'md5'