from django.db import models
from django_extensions.db.models import TimeStampedModel
from django.db.models import TextField, DateField


class AssayDetailedRNA(TimeStampedModel):
    lab = TextField()
    date = DateField()