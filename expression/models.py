from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import CASCADE
from django.urls.base import reverse

from genes.models import Gene, Transcript
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsAutoInitialSaveMixin
from snpdb.models import ImportStatus
from snpdb.models.models_enums import AnnotationLevel


# TODO: Create a class that represents expression, and then subclass for EdgeR and CuffDiff
class CuffDiffFile(GuardianPermissionsAutoInitialSaveMixin, models.Model):
    user = models.ForeignKey(User, on_delete=CASCADE)
    name = models.TextField(null=False)
    annotation_level = models.CharField(max_length=1, choices=AnnotationLevel.CHOICES)  # Gene or transcript
    import_status = models.CharField(max_length=1, choices=ImportStatus.CHOICES, default=ImportStatus.CREATED)
    sample_1 = models.TextField(null=False)
    sample_2 = models.TextField(null=False)
    imported_records = models.IntegerField(null=True)
    matched_records = models.IntegerField(null=True)

    def can_write(self, user):
        return user == self.user

    def __str__(self):
        return f"{self.name} ({self.annotation_level}) - {self.import_status}"

    def get_absolute_url(self):
        return reverse('view_expression_file', kwargs={"expression_file_id": self.pk})


class CuffDiffRecord(models.Model):
    cuff_diff_file = models.ForeignKey(CuffDiffFile, on_delete=CASCADE)
    test_id = models.TextField()
    reference_id = models.TextField()  # GeneVersion or transcript ID
    locus = models.TextField()
    status = models.TextField()
    value_1 = models.FloatField()
    value_2 = models.FloatField()
    log2_fold_change = models.FloatField()
    test_stat = models.FloatField()
    p_value = models.FloatField()
    q_value = models.FloatField()
    # Linked to ensembl tables
    gene = models.ForeignKey(Gene, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, on_delete=CASCADE)
